from molscribe import MolScribe
from PIL import ImageGrab
from PIL import Image
from datetime import datetime
from pathlib import Path
from rdkit import Chem
import chemdraw
import rumps
import torch
import pyperclip
import os
import time
import subprocess
import multiprocessing
import shutil
import csv

import sys
sys.setrecursionlimit(5000)


def save_images_from_clipboard(modus, model=None):
    temp = []
    while True:
        im = ImageGrab.grabclipboard()
        if temp != im:
            if isinstance(im, Image.Image):
                if modus == 'queue':
                    im.save('image_queue/'+str(time.time()) +
                            "_clipboard_image.png", "PNG")
                else:
                    file_path = 'image_queue/' + \
                        str(time.time())+"_clipboard_image.png"
                    im.save(file_path, "PNG")
                    with torch.no_grad():
                        output = model.predict_image_file(
                            file_path, return_atoms_bonds=False, return_confidence=True)
                        smiles = output['smiles']
                        confidence = output['confidence']
                    if confidence > 0.8:
                        pyperclip.copy(smiles)
        temp = im
        # Wait for 0.5 seconds before checking the clipboard again
        time.sleep(0.5)


class ClipboardImageSaverApp(rumps.App):
    def __init__(self):
        # Create a rumps app instance
        self.app = rumps.App(
            "Carlos Helper", icon='pictograms/carlos_helper_logo.png', template=True)

        # Create a MolScribe model instance
        self.model = MolScribe(
            'models--yujieq--MolScribe/snapshots/601bb0f491fb9598fa40227ae9759f800661cd74/swin_base_char_aux_1m.pth', device=torch.device('cpu'))

        # Define image input/output directories
        self.image_input_dir = Path('image_queue/')
        self.image_output_dir = Path('image_output/')
        self.image_generated_dir = Path('image_generated/')

        # Define some app settings
        self.confidence_level = 0.66
        self.max_n_images = 100
        self.length_history = 10
        self.csv_file_path = Path("smiles.csv")

        # Initialize the SMILES history
        self.smiles_history = []

        # Create menu items for saving clipboard images and manipulating SMILES
        self.save_image_item = rumps.MenuItem(
            "Convert Clipboard Image to SMILES", callback=self.save_image, key='V')
        self.smiles_to_image_button = rumps.MenuItem(
            "Convert Clipboard SMILES to Image", callback=self.clipboard_smiles_to_image)
        self.smiles_to_can_smiles_button = rumps.MenuItem(
            "Canonicalize Clipboard SMILES", callback=self.clipboard_smiles_canonical_smiles)

        # Create a menu item for automatically monitoring the clipboard for molecules
        self.automatic_clipboard_button = rumps.MenuItem(
            "Watch Clipboard for Molecules", callback=self.start_clipboard_timer, key='V')

        # Set up a timer to run the save_image function every 0.5 seconds
        self.clipboard_timer = rumps.Timer(self.save_image, 0.5)

        # Create a menu item for showing the SMILES history
        self.history_item = rumps.MenuItem(
            "Copy SMILES from History", callback=self.show_history)

        # Create an "Options" menu for selecting which model to use (currently only MolScribe is implemented)
        self.options = rumps.MenuItem("Options")
        self.options.add(rumps.MenuItem(
            "Use MolScribe", callback=self.select_molscribe))
        self.options.add(rumps.MenuItem("Use Decimer", callback=None))
        self.options[self.options.keys()[0]].state = 1
        self.model_choice = 'molscribe'

        # Set up a continuous monitoring queue
        self.monitor_process = None
        self.description_button = rumps.MenuItem(
            title='Continuous Molecule Detection')
        self.description_button2 = rumps.MenuItem(
            title='Monitor Clipboard to build Image queue')
        self.monitor_button = rumps.MenuItem(
            title='Start Clipboard monitoring', callback=self.start_monitoring)
        self.open_history_file_button = rumps.MenuItem(
            title='Open History File', callback=self.open_smiles_csv)
        self.open_queue_folder_button = rumps.MenuItem(
            title='Add Images to queue', callback=self.open_queue_folder)
        self.process_queue_folder_button = rumps.MenuItem(
            title='Process Image queue', callback=self.process_batch)
        self.history_image_item = rumps.MenuItem(
            title='Copy Image from History', callback=self.process_batch)
        # SMILES to image
        # self.smiles_to_image_button = rumps.MenuItem(title='SMILES to Image', callback=self.smiles_to_image_window)

        # Add the menu items to the app's menu
        self.app.menu = [self.save_image_item,
                         self.smiles_to_image_button,
                         self.smiles_to_can_smiles_button,
                         None,
                         self.description_button,
                         self.automatic_clipboard_button,
                         None,
                         self.open_history_file_button,
                         self.history_item,
                         self.history_image_item,
                         None,
                         self.description_button2,
                         self.monitor_button,
                         self.open_queue_folder_button,
                         self.process_queue_folder_button,
                         None,
                         self.options]

        # Initialize the SMILES history and load the previously saved SMILES strings
        self.read_history()
        self.show_history()

        self.last_image_cache = None

        # Setup the notification center,

        @rumps.notifications
        def notification_center(info):
            if info.data == 'stop_collecting':
                print('stop monitoring')
                self.stop_monitoring('')
            elif info.data == 'stop_collecting_direct':
                self.stop_monitoring_direct('')
            print('info.data ', info.data)

    def start_clipboard_timer(self, sender):
        self.automatic_clipboard_button.title = '🛑 End Molecule watching'
        self.automatic_clipboard_button.set_callback(self.stop_clipboard_timer)

        self.clipboard_timer.start()

    def stop_clipboard_timer(self, sender):
        self.automatic_clipboard_button.title = 'Watch Clipboard for Molecules'
        self.automatic_clipboard_button.set_callback(
            self.start_clipboard_timer)
        self.clipboard_timer.stop()

    def clipboard_smiles_canonical_smiles(self, state):
        clip_conten = pyperclip.paste()
        try:
            print(clip_conten)
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(clip_conten))
            print(smiles)
            pyperclip.copy(smiles)
            print(smiles)
            rumps.notification("✅ Canonicalized SMILES:", 'Copied into Clipboard', str(
                smiles), sound=False, icon='pictograms/carlos_helper_good.png')
        except:
            rumps.notification("🚫 Not valid smiles entered", "provided smiles:  "+str(
                clip_conten), '', sound=False, icon='pictograms/carlos_helper_bad.png')

    def set_icon(self, status, title=None):
        # global icon_state, title_state
        if status == 'busy':
            icon_state = 'pictograms/carlos_helper_logo_busy.png'
            title_state = title
        else:
            icon_state = 'pictograms/carlos_helper_logo.png'
            title_state = None

    def select_decimer(self, sender):
        self.options[self.options.keys()[0]].state = 0
        self.options[self.options.keys()[1]].state = 1
        self.model_choice = 'decimer'

    def select_molscribe(self, sender):
        self.options[self.options.keys()[0]].state = 1
        self.options[self.options.keys()[1]].state = 0
        self.model_choice = 'molscribe'

    def open_smiles_csv(self, sender):
        os.system(f"open {self.csv_file_path}")

    def open_queue_folder(self, sender):
        os.system(f"open {self.image_input_dir}")

    def smiles_to_filename(self, smiles):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        smiles_for_filename = "".join(x for x in smiles if x.isalnum())
        new_filename = f'{timestamp}_{smiles_for_filename}.png'
        return new_filename, timestamp

    def clipboard_to_image(self,):
        im = ImageGrab.grabclipboard()
        if isinstance(im, Image.Image):
            # Save the image as a PNG file
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            file_name = f'{timestamp}_clipboard_image.png'
            file_path = self.image_output_dir / file_name
            im.save(file_path, "PNG")
            print('clipboard save to ', str(self.image_output_dir / file_name))
            return str(file_path), im
        else:
            return False, False

    def clipboard_smiles_to_image(self, sender):
        clip_conten = pyperclip.paste()
        file_name, timestamp, image_path = self.smiles_to_image(clip_conten)
        if file_name:
            self.copy_image_to_clipboard(image_path=image_path)

    def smiles_to_image(self, smiles):
        try:
            file_name, timestamp = self.smiles_to_filename(smiles)
            print(file_name)
            drawer = chemdraw.Drawer(smiles)
            print(drawer)
            image_path = self.image_generated_dir / file_name
            # Write image to file
            fig = drawer.draw()
            fig.write_image(image_path)
            return file_name, timestamp, image_path
        except:
            rumps.notification(
                "🚫 Failed to convert SMILES to image", '', str(smiles), sound=False, icon='pictograms/carlos_helper_bad.png')
            return False, False, False

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

            # Loop through the rows of the CSV file and select all entries from the "smiles" column
            smiles_list = []
            for row in reader:
                smiles_list.append(row[smiles_col_index])
        self.smiles_history = smiles_list[-self.length_history:]

    def start_monitoring(self, _):
        # stop all other montioring processes
        if self.clipboard_timer.is_alive():
            self.stop_clipboard_timer()
            rumps.notification(title='Continuous Molecule detection Stopped', subtitle='',
                               message='The continious queue was started', sound=False, icon='pictograms/carlos_helper_logo.png')

        self.monitor_button.title = '🛑 Stop Monitoring'
        self.monitor_button.set_callback(self.stop_monitoring)
        # Start the subprocess
        self.monitor_process = multiprocessing.Process(
            target=save_images_from_clipboard, args=('queue',))
        self.monitor_process.start()
        notification = rumps.notification(title='Image Collector', subtitle='Collecting images...',
                                          message='Press "Stop" to stop collecting images', sound=True, action_button='Stop', data='stop_collecting', icon='pictograms/carlos_helper_logo.png')

    def stop_monitoring(self, _):
        self.monitor_button.title = 'Start Clipboard monitoring (queue)'
        self.monitor_button.set_callback(self.start_monitoring)
        self.monitor_process.terminate()
        rumps.notification(title='Monitoring Stopped', subtitle='',
                           message='Process the collected images', sound=False, icon='pictograms/carlos_helper_logo.png')
        # self.smiles_to_image_process = multiprocessing.Process(target=self.process_batch)
        self.process_batch('_')

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

    def process_batch(self, _):
        print('Start processing all images in folder')
        self.set_icon('busy')
        for i, filename in enumerate(os.listdir(self.image_input_dir)):
            self.set_icon('busy', str(i)+'/' +
                          str(len(os.listdir(self.image_input_dir))))
            if filename.endswith('.png'):
                self.process_image(filename)

        self.set_icon('idle')

    def process_image(self, filename):
        # generate smiles and rename image
        image_path = self.image_input_dir / filename
        smiles, confidence = self.image_to_smiles(str(image_path))
        if confidence > self.confidence_level:
            new_filename, timestamp, _ = self.smiles_to_image(smiles)
            if new_filename:
                # make smiles string filename save
                new_image_path = self.image_output_dir / new_filename
                shutil.copy2(image_path, new_image_path)

            self.write_to_history(smiles, new_filename, timestamp)
            # remove old image
            os.remove(image_path)

            # delete old images if more than self.max_n_images:
            self.remove_over_max_images()
            self.update_history_menu()
        else:
            # remove old image
            os.remove(image_path)

    def save_image(self, state):
        # Get the image data from the clipboard
        self.set_icon('busy')
        file_path, image = self.clipboard_to_image()

        if file_path and image != self.last_image_cache:
            self.last_image_cache = image
            smiles, confidence = self.image_to_smiles(file_path)
            if confidence > self.confidence_level:
                new_filename, timestamp, _ = self.smiles_to_image(smiles)
                # Add the SMILES string to the history and save it to the file
                self.write_to_history(smiles, new_filename, timestamp)
                pyperclip.copy(smiles)
                # Update the history
                self.update_history_menu()
                rumps.notification(
                    '✅ Image successfully converted', 'SMILES copied into clipboard', '', sound=False, icon='pictograms/carlos_helper_bad.png')
            else:
                rumps.notification(
                    '🔁 Confidence lower than'+str(self.confidence_level), 'Retake screenshot', '', sound=False, icon='pictograms/carlos_helper_bad.png')
        else:
            if not isinstance(state, rumps.Timer):
                rumps.notification('🔁 Clipboard content is not an Image', '',
                                   'Add an Image of molecule to the clipboard', sound=False, icon='pictograms/wrong-file-icon.png')
            self.set_icon('idle')
        self.set_icon('idle')

    def image_to_smiles(self, image):
        if self.model_choice == 'molscribe':
            print('start image analysis')
            with torch.no_grad():
                output = self.model.predict_image_file(
                    image, return_atoms_bonds=False, return_confidence=True)
            smiles = output['smiles']
            confidence = output['confidence']
        elif self.model_choice == 'decimer':
            print('using Decimer')
            # smiles = predict_SMILES(image)
            confidence = 1
        return smiles, confidence

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

        for hist_entry in reversed(self.smiles_history):
            build_function = factory(hist_entry)
            build_function_image = factory_image(hist_entry)
            self.history_item.add(rumps.MenuItem(
                hist_entry, callback=build_function))
            self.history_image_item.add(rumps.MenuItem(
                hist_entry, callback=build_function_image))

    def update_history_menu(self):
        # Update the history menu with the last 10 SMILES strings
        self.history_item.clear()
        self.history_image_item.clear()
        self.read_history()
        self.show_history()


if __name__ == '__main__':
    # Create an instance of the ClipboardImageSaverApp class
    app = ClipboardImageSaverApp()
    app.app.run()
