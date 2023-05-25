from molscribe import MolScribe
from PIL import ImageGrab
from PIL import Image
from pathlib import Path
from rdkit import Chem
import chemdraw
import rumps
import torch
import pyperclip
import os
import subprocess
import shutil
import csv
import rdkit
from datetime import datetime

import re
import sys
from api_conversions import smiles_to_cas_api,smiles_to_iupac_api,cas_to_price_api, iupac_to_smiles_api,cas_to_smiles_api,cas_to_iupac_api,iupac_to_cas_api
from utils import create_filename_from_smiles,smiles_to_ms_peaks
import requests

print(requests.certs.where())
sys.setrecursionlimit(5000)


class ClipboardImageSaverApp(rumps.App):
    def __init__(self):

        # Create a rumps app instance
        self.app = rumps.App(
            "Clipboard2Smiles",
            icon='pictograms/carlos_helper_logo.png',
            template=True)

        # Define some app settings
        self.confidence_level = 0.66
        self.max_n_images = 100
        self.length_history = 10
        self.model_path = 'models--yujieq--MolScribe/snapshots/601bb0f491fb9598fa40227ae9759f800661cd74/swin_base_char_aux_1m.pth'
        self.model_device = 'cpu'

        # Define image input/output directories
        self.image_input_dir = Path('image_queue/')
        self.image_output_dir = Path('image_output/')
        self.image_generated_dir = Path('image_generated/')

        # Define history file and initialize empty history
        self.csv_file_path = Path("smiles.csv")
        self.smiles_history = []
        self.purchase_link = ''
        self.internet_connection = False


        # Create a MolScribe model instance
        self.model = MolScribe(
            self.model_path, device=torch.device(self.model_device))

        # Create menu items for saving clipboard images and manipulating SMILES
        self.to_smiles_button = rumps.MenuItem(
            "Convert Clipboard to SMILES", callback=self.clipboard_to_smiles, key='V')
        self.to_image_button = rumps.MenuItem(
            "Convert Clipboard to Image", callback=self.clipboard_to_image)
        self.smiles_to_can_smiles_button = rumps.MenuItem(
            "Canonicalize Clipboard SMILES", callback=self.clipboard_smiles_to_canonical_smiles)

        # Api Related checks
        self.to_cas_button = rumps.MenuItem(
            "Convert Clipboard to CAS", callback=self.clipboard_to_cas)
        
        self.to_iupac_button = rumps.MenuItem(
            "Convert Clipboard to IUPAC", callback=self.clipboard_to_iupac)

        self.to_price_button = rumps.MenuItem(
            "Convert Clipboard to Price", callback=self.clipboard_to_price)
        

        # Create a menu item for automatically monitoring the clipboard for molecules
        self.clipboard_to_smiles_watcher_button = rumps.MenuItem(
            "Watch Clipboard for Molecules", callback=self.start_clipboard_to_smiles_timer)

        # Create a menu item for showing the SMILES history
        self.open_history_file_button = rumps.MenuItem(
            title='Open History File', callback=self.open_history_file)
        self.history_item_smiles = rumps.MenuItem(
            "Copy SMILES from History", callback=self.show_history)
        
        self.history_item_price = rumps.MenuItem(
            "Find price from History", callback=self.show_history)
        self.history_image_item = rumps.MenuItem(
            title='Copy Image from History', callback=self.show_history)
        
        self.history_ms_item = rumps.MenuItem(
            title='Common Ions in MS from History', callback=self.show_history)

        # Initialize the SMILES history and load the previously saved SMILES strings
        self.read_history()
        self.show_history()

        # Set up a continuous monitoring queue
        self.monitor_process = None
        self.description_detection_button = rumps.MenuItem(
            title='Continuous Molecule Detection')
        self.description_image_queue_button = rumps.MenuItem(
            title='Monitor Clipboard to build Image queue')

        self.clipboard_to_image_queue_button = rumps.MenuItem(
            title='Start Clipboard monitoring', callback=self.start_clipboard_image_queue_timer)
        self.open_queue_folder_button = rumps.MenuItem(
            title='Add Images to queue', callback=self.open_queue_folder)
        self.process_queue_folder_button = rumps.MenuItem(
            title='Process Image queue', callback=self.batch_images_to_smiles)

        # Image cache, to only run conversion when a new Image is detected
        self.last_image_cache = None

        # Set up a timer to run the save_image function every 0.5 seconds
        self.clipboard_to_smiles_timer = rumps.Timer(
            self.clipboard_image_to_smiles, 0.5)
        self.clipboard_to_image_queue_timer = rumps.Timer(
            self.get_image_from_clipboard, 0.5)

        # Add the menu items to the app's menu
        self.app.menu = [self.to_smiles_button,
                         self.to_image_button,
                         self.to_cas_button,
                         self.to_iupac_button,
                         self.smiles_to_can_smiles_button,
                         None,
                         self.to_price_button,
                         None,
                         self.description_detection_button,
                         self.clipboard_to_smiles_watcher_button,
                         None,
                         self.open_history_file_button,
                         self.history_item_smiles,
                         self.history_image_item,
                         self.history_item_price,
                         self.history_ms_item,
                         # self.history_price_item,
                         None,
                         self.description_image_queue_button,
                         self.clipboard_to_image_queue_button,
                         self.open_queue_folder_button,
                         self.process_queue_folder_button,
                         None]

        # Setup the notification center,
        @rumps.notifications
        def notification_center(info):
            if info.data == 'stop_collecting':
                self.stop_clipboard_image_queue_timer('')
            elif info.data == 'open_purchase_link':
                os.system("open \"\" "+self.purchase_link)

    def start_clipboard_to_smiles_timer(self, sender):
        self.clipboard_to_smiles_watcher_button.title = '🛑 End Molecule watching'
        self.clipboard_to_smiles_watcher_button.set_callback(
            self.stop_clipboard_to_smiles_timer)
        self.clipboard_to_smiles_timer.start()


    def clipoard_content_indentification(self):
        clip = ImageGrab.grabclipboard()

        if isinstance(clip, Image.Image):
            return 'image'
        
        clip = pyperclip.paste()
        if clip:
            try:
                mol = Chem.MolFromSmiles(clip)
                if mol:
                    return 'smiles'
            except:
                print('Not SMILES')

            if re.match(r'^\d+-\d+-\d+$', clip):
                return 'cas'
            
    def clipboard_to_price(self,sender):
        clipboard_input_type = self.clipoard_content_indentification()
        print('clipboard_input_type',clipboard_input_type)
        cas = False
        if clipboard_input_type == 'image':
            smiles = self.clipboard_image_to_smiles(notifications=False)
            cas = smiles_to_cas_api(smiles)
        elif clipboard_input_type == 'smiles':
            smiles = pyperclip.paste()
            cas = smiles_to_cas_api(smiles)
        elif clipboard_input_type == 'cas':
            cas = pyperclip.paste()
        else:
            print('try iupac conversion')
            cas = iupac_to_cas_api(pyperclip.paste())

        if cas:
            print('cas',cas)
            price_dict = cas_to_price_api(cas,vendor='ChemieBrunschwieg')
            self.purchase_link = price_dict['link']
            rumps.notification(f'💰 Buy {price_dict["ItemName"]}',
                f'{price_dict["Packsize"]} for {price_dict["Preis"]} CHF',f'or {price_dict["price_per"]:2} CHF/g',
            action_button='Open Link',
            data='open_purchase_link', icon='pictograms/carlos_helper_good.png')
        elif clipboard_input_type !='image':
            rumps.notification('🔁 Clipboard content could not be converted into CAS', '',
                'The Clipboard content provided was:'+str(pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_bad.png')

    def clipboard_to_iupac(self,sender):
        clipboard_input_type = self.clipoard_content_indentification()
        iupac = False
        if clipboard_input_type == 'image':
            smiles = self.clipboard_image_to_smiles(sender,notifications=False)
            iupac = smiles_to_iupac_api(smiles)
            if iupac:
                pyperclip.copy(iupac)
                rumps.notification("✅ IUPAC found for Image: ", 'Copied into Clipboard', 
                                    str(iupac), sound=False, icon='pictograms/carlos_helper_good.png')
            else: 
                rumps.notification("🚫 Couldn't convert Clipboard into IUPAC "+str(pyperclip.paste()), 'No IUPAC name found for', str(
                pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')

        elif clipboard_input_type == 'smiles':
            smiles = pyperclip.paste()
            iupac = smiles_to_iupac_api(smiles)
            if iupac:
                pyperclip.copy(iupac)
                self.smiles_to_history(smiles)
                rumps.notification("✅ IUPAC found for SMILES: "+str(smiles), 'Copied into Clipboard', 
                                    str(iupac), sound=False, icon='pictograms/carlos_helper_good.png')
            else: 
                rumps.notification("🚫 Couldn't convert Clipboard into IUPAC "+str(pyperclip.paste()), 'No IUPAC name found for', str(
                pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')

        elif clipboard_input_type == 'cas':
            cas = pyperclip.paste()
            iupac = cas_to_iupac_api(cas)
            if iupac:
                pyperclip.copy(iupac)
                rumps.notification("✅ IUPAC found for CAS: "+str(cas), 'Copied into Clipboard', 
                                    str(iupac), sound=False, icon='pictograms/carlos_helper_good.png')
                
                smiles = cas_to_smiles_api(cas)
                self.smiles_to_history(smiles)
            else: 
                rumps.notification("🚫 Couldn't convert Clipboard into IUPAC "+str(pyperclip.paste()), 'No IUPAC name found for', str(
                pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')
        else:
            rumps.notification("🚫 Couldn't convert Clipboard into IUPAC "+str(pyperclip.paste()), 'No IUPAC name found for', str(
            pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')

    def clipboard_to_image(self,sender):
        smiles = self.clipboard_to_smiles(sender)
        if smiles:
            file_name, timestamp, image_path = self.smiles_to_image(smiles)
            if file_name:
                self.copy_image_to_clipboard(image_path=image_path)

    def clipboard_to_smiles(self,sender,notifications=True):
        clipboard_input_type = self.clipoard_content_indentification()
        smiles = False
        if clipboard_input_type == 'image':
            smiles = self.clipboard_image_to_smiles(sender,notifications=False)
            return smiles
        elif clipboard_input_type == 'smiles':
            smiles =  pyperclip.paste()
            self.smiles_to_history(smiles)
            return smiles
        elif clipboard_input_type == 'cas':
            cas = pyperclip.paste()
            smiles = cas_to_smiles_api(cas)
            if smiles:
                pyperclip.copy(smiles)
                self.smiles_to_history(smiles)
                if notifications:
                    rumps.notification("✅ SMILES found for CAS: "+str(cas), 'Copied into Clipboard', str(
                        smiles), sound=False, icon='pictograms/carlos_helper_good.png')
                return smiles
            else:
                if notifications:
                    rumps.notification("🚫 Couldn't convert Clipboard into SMILES "+str(pyperclip.paste()), 'No SMILES found for', str(
                    pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')
        else:
            print('try iupac conversion')
            smiles = iupac_to_smiles_api(pyperclip.paste())
            if smiles:
                if notifications:
                    rumps.notification("✅ SMILES found for IUPAC: "+str(pyperclip.paste()), 'Copied into Clipboard', str(
                        smiles), sound=False, icon='pictograms/carlos_helper_good.png')
                pyperclip.copy(smiles)
                self.smiles_to_history(smiles)
                return smiles
            else:
                if notifications:
                    rumps.notification("🚫 Couldn't convert Clipboard into SMILES "+str(pyperclip.paste()), 'No SMILES found for', str(
                    pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')

    def clipboard_to_cas(self,sender):
        clipboard_input_type = self.clipoard_content_indentification()
        cas = False
        print('clipboard_input_type',clipboard_input_type, pyperclip.paste(),)
        if clipboard_input_type == 'image':
            smiles = self.clipboard_image_to_smiles(sender,notifications=False)
            cas = smiles_to_cas_api(smiles)
            if cas:
                pyperclip.copy(cas)
                rumps.notification("✅ CAS found for Image: ", 'Copied into Clipboard', 
                                    str(cas), sound=False, icon='pictograms/carlos_helper_good.png')
        elif clipboard_input_type == 'smiles':
            smiles = pyperclip.paste()
            cas = smiles_to_cas_api(smiles)
            if cas:
                pyperclip.copy(cas)
                self.smiles_to_history(smiles)
                rumps.notification("✅ CAS found for SMILES: ", 'Copied into Clipboard', 
                                    str(cas), sound=False, icon='pictograms/carlos_helper_good.png')
            else:
                rumps.notification("🚫 Couldn't convert Clipboard into CAS "+str(pyperclip.paste()), 'No CAS found for', str(
                pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')


        elif clipboard_input_type == 'cas':
            smiles = cas_to_smiles_api(pyperclip.paste())
            if smiles:
                self.smiles_to_history(smiles)
                rumps.notification("✅ CAS already in clipboard "+str(pyperclip.paste()), 'Copied into Clipboard', str(
                    pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')
        else:
            cas = iupac_to_cas_api(pyperclip.paste())
            if cas:
                rumps.notification("✅ CAS found for IUPAC: "+str(pyperclip.paste()), 'Copied into Clipboard', str(
                    cas), sound=False, icon='pictograms/carlos_helper_good.png')
                pyperclip.copy(cas)
                smiles = cas_to_smiles_api(cas)
                self.smiles_to_history(smiles)
            else:
                rumps.notification("🚫 Couldn't convert Clipboard into CAS "+str(pyperclip.paste()), 'No CAS found for', str(
                pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')

    def clipboard_smiles_to_canonical_smiles(self, state):
        clip_content = pyperclip.paste()
        try:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(clip_content))
            pyperclip.copy(smiles)
            rumps.notification("✅ Canonicalized SMILES:", 'Copied into Clipboard', str(
                smiles), sound=False, icon='pictograms/carlos_helper_good.png')
            return smiles 
        except:
            rumps.notification("🚫 Not valid smiles entered", "provided smiles:  "+str(
                clip_content), '', sound=False, icon='pictograms/carlos_helper_bad.png')



    def stop_clipboard_to_smiles_timer(self, sender):
        self.clipboard_to_smiles_watcher_button.title = 'Watch Clipboard for Molecules'
        self.clipboard_to_smiles_watcher_button.set_callback(
            self.start_clipboard_to_smiles_timer)
        self.clipboard_to_smiles_timer.stop()

    def start_clipboard_image_queue_timer(self, sender):
        self.clipboard_to_image_queue_button.title = '🛑 End Molecule watching'
        self.clipboard_to_image_queue_button.set_callback(
            self.stop_clipboard_image_queue_timer)
        rumps.notification(title='Image Collector', subtitle='Collecting images...',
                           message='Press "Stop" to stop collecting images',
                           sound=True, action_button='Stop',
                           data='stop_collecting',
                           icon='pictograms/carlos_helper_logo.png')
        self.clipboard_to_image_queue_timer.start()

    def stop_clipboard_image_queue_timer(self, sender):
        self.clipboard_to_image_queue_button.title = 'Watch Clipboard for Molecules'
        self.clipboard_to_image_queue_button.set_callback(
            self.start_clipboard_image_queue_timer)
        self.clipboard_to_image_queue_timer.stop()
        self.batch_images_to_smiles('_')


    def get_image_from_clipboard(self, state):
        im = ImageGrab.grabclipboard()

        if isinstance(im, Image.Image) and im != self.last_image_cache:
            # Save the image as a PNG file
            self.last_image_cache = im
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            file_name = f'{timestamp}_clipboard_image.png'
            file_path = self.image_input_dir / file_name
            im.save(file_path, "PNG")
            print('clipboard save to ', file_path)
            return str(file_path)
        else:
            return False

    def copy_image_to_clipboard(self, image_path):
        try:
            temp_path = os.path.abspath(image_path)
            subprocess.run(
                ["osascript", "-e", 'set the clipboard to (read (POSIX file "'+str(temp_path)+'") as «class PNGf»)'])
            rumps.notification("✅ Successfully draw molecule", '',
                               "Image in clipboard", sound=False, icon='pictograms/carlos_helper_good.png')
        except:
            rumps.notification("🚫 Failed to load Image into clipboard",
                               '', str(image_path), sound=False,  icon='pictograms/carlos_helper_bad.png')

    def clipboard_image_to_smiles(self, state, notifications=True):
        # Get the image data from the clipboard
        file_path = self.get_image_from_clipboard(state)
        if isinstance(state, rumps.Timer):
            notifications = False

        if file_path:
            return self.image_to_smiles(
                file_path, save_smiles_to_clipboard=True, notifications=notifications)
        else:
            if notifications:
                rumps.notification('🔁 Clipboard content is not an Image', '',
                                   'Add an Image of molecule to the clipboard', sound=False, icon='pictograms/carlos_helper_bad.png')

    def image_to_smiles(self, image_path, save_smiles_to_clipboard=True, notifications=True):
        # generate smiles and rename image
        # image_path = self.image_input_dir / filename
        with torch.no_grad():
            output = self.model.predict_image_file(
                image_path, return_atoms_bonds=False, return_confidence=True)
        smiles = output['smiles']
        confidence = output['confidence']
        if confidence > self.confidence_level:
            new_filename, timestamp, _ = self.smiles_to_image(smiles)
            if new_filename:
                # Save original image in the processes image output folder
                new_image_path = self.image_output_dir / new_filename
                shutil.copy2(image_path, new_image_path)
            # remove original image (either from queue folder or clipboard temp)
            os.remove(image_path)

            # delete old images if more than self.max_n_images:
            self.remove_over_max_images()
            self.write_to_history(smiles, new_filename, timestamp)
            self.update_history_menu()

            if save_smiles_to_clipboard:
                pyperclip.copy(smiles)
                if notifications:
                    rumps.notification(
                        '✅ Image successfully converted', 'SMILES copied into clipboard', '', sound=False, icon='pictograms/carlos_helper_bad.png')
            return smiles
        else:
            os.remove(image_path)
            if notifications:
                rumps.notification(
                    '🔁 Confidence lower than'+str(self.confidence_level), 'Retake screenshot', '', sound=False, icon='pictograms/carlos_helper_bad.png')

    def smiles_to_image(self, smiles, notifications=True):
        try:
            file_name, timestamp = create_filename_from_smiles(smiles)
            drawer = chemdraw.Drawer(smiles)
            image_path = self.image_generated_dir / file_name
            # Write image to file
            drawer.draw_img(image_path)
            return file_name, timestamp, image_path
        except:
            if notifications:
                rumps.notification(
                    "🚫 Failed to convert SMILES to image", '', str(smiles), sound=False, icon='pictograms/carlos_helper_bad.png')
            return False, False, False

    def remove_over_max_images(self,):
        # delete old images if more than self.max_n_images:
        all_files = sorted(self.image_output_dir.glob('*.png'))
        if len(all_files) > self.max_n_images:
            for old_file in all_files[:-self.max_n_images:]:
                os.remove(old_file)
        all_files = sorted(self.image_generated_dir.glob('*.png'))
        if len(all_files) > self.max_n_images:
            for old_file in all_files[:-self.max_n_images:]:
                os.remove(old_file)

    def batch_images_to_smiles(self, _):
        print('Start processing all images in folder')
        for i, filename in enumerate(os.listdir(self.image_input_dir)):
            if filename.endswith('.png'):
                image_path = str(self.image_input_dir / str(filename))
                self.image_to_smiles(
                    image_path, save_smiles_to_clipboard=False, notifications=False)

    def smiles_to_history(self, smiles):
        new_filename, timestamp, _ = self.smiles_to_image(smiles,notifications=False)
        self.write_to_history(smiles, new_filename, timestamp)
        self.update_history_menu()

    def write_to_history(self, smiles, new_filename, timestamp):
        with open(self.csv_file_path, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([smiles, new_filename, timestamp])

    def read_history(self):
        with open(self.csv_file_path, 'r') as csvfile:
            reader = csv.reader(csvfile)

            # Get the header row and find the index of the "smiles" column
            header_row = next(reader)
            smiles_col_index = header_row.index("smiles")
            image_file_col_index = header_row.index(" image_file")

            # Loop through the rows of the CSV file and select all entries from the "smiles" column
            smiles_list = []
            for row in reader:
                smiles_list.append((row[smiles_col_index],row[image_file_col_index]))
        self.smiles_history = smiles_list[-self.length_history:]

    def show_history(self):
        # Display the SMILES history
        def factory(smiles):
            def f(smiles):
                pyperclip.copy(smiles.title)
            return f

        def factory_image(smiles):
            def f(smiles):
                file_name, timestamp, image_path = self.smiles_to_image(
                    smiles.title)
                if file_name:
                    self.copy_image_to_clipboard(image_path=image_path)
            return f

        def factory_ms(value):
            def f(value):
                pyperclip.copy(value.title)
                return
            return f

        def factory_price(value):
            def f(value):
                pyperclip.copy(value.title)
                self.clipboard_to_price('')
                return
            return f

        for hist_smi, hist_img_file in reversed(self.smiles_history):
            build_function_smiles = factory(hist_smi)
            build_function_image = factory_image(hist_smi)
            build_function_price = factory_price(hist_smi)
            ms_dict = smiles_to_ms_peaks(hist_smi)
            if hist_img_file=='False':
                self.history_item_smiles.add(rumps.MenuItem(
                hist_smi, callback=build_function_smiles))
                self.history_image_item.add(rumps.MenuItem(
                    hist_smi, callback=build_function_image,))
                self.history_item_price.add(rumps.MenuItem(
                    hist_smi, callback=build_function_price,))
            else:
                self.history_image_item.add(rumps.MenuItem(
                    hist_smi, callback=build_function_image, icon=str(self.image_generated_dir/hist_img_file),dimensions=[75,75]))
                self.history_item_smiles.add(rumps.MenuItem(
                    hist_smi, callback=build_function_smiles,icon=str(self.image_generated_dir/hist_img_file),dimensions=[75,75]))
                self.history_item_price.add(rumps.MenuItem(
                    hist_smi, callback=build_function_price,icon=str(self.image_generated_dir/hist_img_file),dimensions=[75,75]))


            if ms_dict:
                temp_history = rumps.MenuItem(
                    hist_smi, callback=build_function_smiles)
                for key in ms_dict:
                    left_aligned = f"{key:<{9}}"
                    temp = f'{ms_dict[key]:.3f}'
                    right_aligned = f"{temp:>{5}}"
                    build_function_ms = factory_ms(ms_dict[key])
                    temp_history.add(rumps.MenuItem(
                        left_aligned+right_aligned, callback=build_function_ms))
            self.history_ms_item.add(temp_history)

    def update_history_menu(self):
        # Update the history menu with the last 10 SMILES strings
        self.history_item_smiles.clear()
        self.history_image_item.clear()
        self.history_ms_item.clear()
        self.history_item_price.clear()
        self.read_history()
        self.show_history()

    def open_history_file(self, sender):
        os.system(f"open {self.csv_file_path}")

    def open_queue_folder(self, sender):
        os.system(f"open {self.image_input_dir}")

            
    


if __name__ == '__main__':
    # Create an instance of the ClipboardImageSaverApp class

    app = ClipboardImageSaverApp()
    app.app.run()
