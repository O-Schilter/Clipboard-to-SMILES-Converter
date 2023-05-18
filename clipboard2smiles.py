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
import subprocess
import shutil
import csv
import rdkit
import sys

from rdkit.Chem import Descriptors
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

        # Create a MolScribe model instance
        self.model = MolScribe(
            self.model_path, device=torch.device(self.model_device))

        # Create menu items for saving clipboard images and manipulating SMILES
        self.image_to_smiles_button = rumps.MenuItem(
            "Convert Clipboard Image to SMILES", callback=self.clipboard_image_to_smiles, key='V')
        self.smiles_to_image_button = rumps.MenuItem(
            "Convert Clipboard SMILES to Image", callback=self.clipboard_smiles_to_image)
        self.smiles_to_can_smiles_button = rumps.MenuItem(
            "Canonicalize Clipboard SMILES", callback=self.clipboard_smiles_to_canonical_smiles)

        # Create a menu item for automatically monitoring the clipboard for molecules
        self.clipboard_to_smiles_watcher_button = rumps.MenuItem(
            "Watch Clipboard for Molecules", callback=self.start_clipboard_to_smiles_timer, key='V')

        # Create a menu item for showing the SMILES history
        self.open_history_file_button = rumps.MenuItem(
            title='Open History File', callback=self.open_history_file)
        self.history_item_smiles = rumps.MenuItem(
            "Copy SMILES from History", callback=self.show_history)
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
        self.app.menu = [self.image_to_smiles_button,
                         self.smiles_to_image_button,
                         self.smiles_to_can_smiles_button,
                         None,
                         self.description_detection_button,
                         self.clipboard_to_smiles_watcher_button,
                         None,
                         self.open_history_file_button,
                         self.history_item_smiles,
                         self.history_image_item,
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

    def start_clipboard_to_smiles_timer(self, sender):
        self.clipboard_to_smiles_watcher_button.title = '🛑 End Molecule watching'
        self.clipboard_to_smiles_watcher_button.set_callback(
            self.stop_clipboard_to_smiles_timer)
        self.clipboard_to_smiles_timer.start()

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

    def open_queue_folder(self, sender):
        os.system(f"open {self.image_input_dir}")

    def create_filename_from_smiles(self, smiles):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        smiles_for_filename = "".join(x for x in smiles if x.isalnum())
        new_filename = f'{timestamp}_{smiles_for_filename}.png'
        return new_filename, timestamp

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

    def clipboard_image_to_smiles(self, state):
        # Get the image data from the clipboard
        file_path = self.get_image_from_clipboard()
        print(file_path,)
        if not isinstance(state, rumps.Timer):
            notifications = True
        else:
            notifications = False

        if file_path:
            print(file_path)
            self.image_to_smiles(
                file_path, save_smiles_to_clipboard=True, notifications=notifications)
        else:
            if notifications:
                rumps.notification('🔁 Clipboard content is not an Image', '',
                                   'Add an Image of molecule to the clipboard', sound=False, icon='pictograms/carlos_helper_bad.png')

    def clipboard_smiles_to_image(self, sender):
        clip_content = pyperclip.paste()
        file_name, timestamp, image_path = self.smiles_to_image(clip_content)
        if file_name:
            self.copy_image_to_clipboard(image_path=image_path)

    def clipboard_smiles_to_canonical_smiles(self, state):
        clip_content = pyperclip.paste()
        try:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(clip_content))
            pyperclip.copy(smiles)
            rumps.notification("✅ Canonicalized SMILES:", 'Copied into Clipboard', str(
                smiles), sound=False, icon='pictograms/carlos_helper_good.png')
        except:
            rumps.notification("🚫 Not valid smiles entered", "provided smiles:  "+str(
                clip_content), '', sound=False, icon='pictograms/carlos_helper_bad.png')

    def smiles_to_image(self, smiles):
        try:
            file_name, timestamp = self.create_filename_from_smiles(smiles)
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

    def image_to_smiles(self, image_path, save_smiles_to_clipboard=True, notifications=True):
        # generate smiles and rename image
        # image_path = self.image_input_dir / filename
        print(image_path)
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
                rumps.notification(
                    '✅ Image successfully converted', 'SMILES copied into clipboard', '', sound=False, icon='pictograms/carlos_helper_bad.png')

        else:
            os.remove(image_path)
            if notifications:
                rumps.notification(
                    '🔁 Confidence lower than'+str(self.confidence_level), 'Retake screenshot', '', sound=False, icon='pictograms/carlos_helper_bad.png')

    def smiles_to_ms_peaks(self, smiles):
        try:
            mol_weight = Descriptors.ExactMolWt(
                rdkit.Chem.MolFromSmiles(smiles))
            return {'[M]': mol_weight,
                    '[M+H]': mol_weight+1.007276,
                    '[M+NH4]': mol_weight+18.033823,
                    '[M+Na]': mol_weight+22.989218,
                    '[M+K]': mol_weight+38.963158,
                    '[2M+H]': 2*mol_weight+1.007276,
                    '[2M+NH4]': 2*mol_weight+18.033823,
                    '[2M+Na]': 2*mol_weight+22.989218,
                    '[2M+K]': 2*mol_weight+38.963158}
        except:
            return False

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
                return
            return f

        for hist_entry in reversed(self.smiles_history):
            build_function_smiles = factory(hist_entry)
            build_function_image = factory_image(hist_entry)
            ms_dict = self.smiles_to_ms_peaks(hist_entry)

            self.history_item_smiles.add(rumps.MenuItem(
                hist_entry, callback=build_function_smiles))
            self.history_image_item.add(rumps.MenuItem(
                hist_entry, callback=build_function_image))

            if ms_dict:
                temp_history = rumps.MenuItem(
                    hist_entry, callback=build_function_smiles)
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
        self.read_history()
        self.show_history()

    def open_history_file(self, sender):
        os.system(f"open {self.csv_file_path}")


if __name__ == '__main__':
    # Create an instance of the ClipboardImageSaverApp class

    app = ClipboardImageSaverApp()
    app.app.run()
