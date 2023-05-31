import sys
import csv
import os
import re
import subprocess

import rumps
import torch
import pyperclip

from functools import partial
from pathlib import Path
from datetime import datetime

from PIL import ImageGrab
from PIL import Image

import selfies as sf
from rdkit import Chem
from rdkit.Chem import AllChem

from molscribe import MolScribe

from utils import smiles_to_molecular_properties
from converter import Converter

sys.setrecursionlimit(5000)


class ClipboardImageSaverApp(rumps.App):
    def __init__(self):

        # Create a rumps app instance
        self.app = rumps.App(
            "Clipboard2Smiles",
            icon='pictograms/carlos_helper_logo.png',
            template=True)

        # Define app settings
        self.confidence_level = 0.66
        self.max_n_images = 25
        self.length_history = 8
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

        # Create Converter class
        self.converter = Converter(model=self.model,
                                   image_input_dir=self.image_input_dir, image_generated_dir=self.image_generated_dir)

        # Create menu items for chemical structure conversion
        self.to_smiles_button = rumps.MenuItem(
            "Convert Clipboard to SMILES", callback=partial(self.clipboard_to, output_format='smiles'))
        self.to_image_button = rumps.MenuItem(
            "Convert Clipboard to Image", callback=partial(self.clipboard_to, output_format='image'))
        self.to_cas_button = rumps.MenuItem(
            "Convert Clipboard to CAS", callback=partial(self.clipboard_to, output_format='cas'))
        self.to_iupac_button = rumps.MenuItem(
            "Convert Clipboard to IUPAC", callback=partial(self.clipboard_to, output_format='iupac'))
        self.to_inchi_button = rumps.MenuItem(
            "Convert Clipboard to InChI", callback=partial(self.clipboard_to, output_format='inchi'))
        self.to_mol_button = rumps.MenuItem(
            "Convert Clipboard to RDkit Mol", callback=partial(self.clipboard_to, output_format='mol'))
        self.to_selfies_button = rumps.MenuItem(
            "Convert Clipboard to Selfies", callback=partial(self.clipboard_to, output_format='selfies'))

        # Clipboard SMILES operations
        self.smiles_to_smiles_button = rumps.MenuItem(
            "Clipboard SMILES Operations")
        self.smiles_to_smiles_button.add(rumps.MenuItem(
            "Canonicalize SMILES", callback=partial(self.clipboard_to, output_format='canonicalize')))
        self.smiles_to_smiles_button.add(rumps.MenuItem(
            "Augment SMILES", callback=partial(self.clipboard_to, output_format='augment')))
        # self.smiles_to_smiles_button.add(rumps.MenuItem(
        #     "Standardize Reaction SMILES", callback=partial(self.clipboard_to, output_format='standardize')))
        self.smiles_to_smiles_button.add(rumps.MenuItem(
            "Remove Atom Mapping from SMILES",
            callback=partial(self.clipboard_to, output_format='remove_atom_mapping')))

        # Finde Price for molecules
        self.to_price_button = rumps.MenuItem(
            "Convert Clipboard to Price", callback=partial(self.clipboard_to, output_format='price'))
        self.history_item_price = rumps.MenuItem(
            "Find Price from History", callback=self.update_history_menu)

        # Create a menu item for automatically monitoring the clipboard for molecules
        self.clipboard_to_smiles_watcher_button = rumps.MenuItem(
            "Watch Clipboard for Molecules", callback=self.start_clipboard_to_smiles_timer)

        # Create a menu item for showing the history

        self.history_item_smiles = rumps.MenuItem(
            "Copy SMILES from History", callback=self.update_history_menu)
        self.history_structure_item = rumps.MenuItem(
            title='Copy Structure from History', callback=self.update_history_menu)
        self.history_molecular_properties_item = rumps.MenuItem(
            title='Molecular Properties from History', callback=self.update_history_menu)
        self.open_history_file_button = rumps.MenuItem(
            title='Open History File', callback=self.open_history_file)

        # Initialize the SMILES history and load the previously saved molecules
        self.read_history()
        self.show_history()

        # Set up a continuous monitoring queue
        self.description_detection_button = rumps.MenuItem(
            title='Continuous Molecule Detection')
        self.description_image_queue_button = rumps.MenuItem(
            title='Monitor Clipboard to Build Image Queue')

        self.clipboard_to_image_queue_button = rumps.MenuItem(
            title='Start Clipboard Monitoring', callback=self.start_clipboard_image_queue_timer)
        self.open_queue_folder_button = rumps.MenuItem(
            title='Add Images to Queue', callback=self.open_queue_folder)
        self.process_queue_folder_button = rumps.MenuItem(
            title='Process Image Queue', callback=self.batch_images_to_smiles)

        # Options
        self.options_button = rumps.MenuItem(title='Options')

        # Select Vendor Menu
        temp_option_menu = rumps.MenuItem(title='Select Vendor')
        temp_option_menu.add(rumps.MenuItem(
            title='Enamine', callback=self.select_vendor))
        temp_option_menu.add(rumps.MenuItem(
            title='Chemie Brunschwieg', callback=self.select_vendor))
        temp_option_menu['Enamine'].state = 1
        self.options_button.add(temp_option_menu)

        # Image cache, to only run conversion when a new Image is detected
        self.last_image_cache = None
        self.last_clipboard_cache = None

        # Set up a timer to run the save_image function every 0.5 seconds
        self.clipboard_to_smiles_timer = rumps.Timer(
            partial(self.clipboard_to, output_format='smiles'), 0.5)
        self.clipboard_to_image_queue_timer = rumps.Timer(
            self.save_image_from_clipboard, 0.5)

        # Add the menu items to the app's menu
        self.app.menu = [self.to_smiles_button,
                         self.to_image_button,
                         self.to_cas_button,
                         self.to_iupac_button,
                         self.to_inchi_button,
                         self.to_mol_button,
                         self.to_selfies_button,
                         None,
                         self.smiles_to_smiles_button,
                         None,
                         self.to_price_button,
                         self.history_item_price,
                         None,
                         self.description_detection_button,
                         self.clipboard_to_smiles_watcher_button,
                         None,
                         self.history_item_smiles,
                         self.history_structure_item,
                         self.history_molecular_properties_item,
                         self.open_history_file_button,
                         None,
                         self.description_image_queue_button,
                         self.clipboard_to_image_queue_button,
                         self.open_queue_folder_button,
                         self.process_queue_folder_button,
                         None,
                         self.options_button]

        # Setup the notification center,
        @rumps.notifications
        def notification_center(info):
            if info.data == 'stop_collecting':
                self.stop_clipboard_image_queue_timer('')
            elif info.data == 'open_purchase_link':
                os.system("open \"\" " + self.purchase_link)

    def select_vendor(self, sender):
        self.converter.vendors.select_vendor(sender.title.replace(' ', ''))
        for menu_item in self.options_button['Select Vendor']:
            if menu_item == sender.title:
                self.options_button['Select Vendor'][menu_item].state = 1
            else:
                self.options_button['Select Vendor'][menu_item].state = 0

    def start_clipboard_to_smiles_timer(self, sender):
        self.clipboard_to_smiles_watcher_button.title = '🛑 End Molecule watching'
        self.clipboard_to_smiles_watcher_button.set_callback(
            self.stop_clipboard_to_smiles_timer)
        self.clipboard_to_smiles_timer.start()

    def clipboard_content_identification(self):
        """
        Identifies the content type of the clipboard and returns a dictionary with the format and content.

        Supported formats:
        - Image
        - SMILES
        - InChI
        - MOL
        - Selfies
        - Reaction SMILES
        - CAS number
        - IUPAC name

        Returns:
        dict: A dictionary with the format ('format') and content ('content') of the clipboard.
        """
        clip = ImageGrab.grabclipboard()

        if isinstance(clip, Image.Image):
            return {'format': 'image', 'content': clip}

        clip = pyperclip.paste()

        if clip:
            try:
                mol = Chem.MolFromSmiles(clip)
                if mol:
                    return {'format': 'smiles', 'content': clip}
            except:
                print('Not SMILES')

            if clip.startswith("InChI="):
                try:
                    mol = Chem.MolFromInchi(clip)
                    if mol is not None:
                        return {'format': 'inchi', 'content': clip}
                except:
                    print("Invalid InChI")

            try:
                mol = Chem.MolFromMolBlock(clip)
                if mol is not None:
                    return {'format': 'mol', 'content': clip}
            except:
                print("Invalid MOL file")

            try:
                smi = sf.decoder(clip)
                if smi is not None and smi != '' and len(smi) > 0:
                    print('Selfies converted to SMILES:', smi)
                    return {'format': 'selfies', 'content': clip}
            except:
                print("Invalid Selfies")

            if '>' in clip:
                print('">" detected')
                try:
                    rxn = AllChem.ReactionFromSmarts(clip)
                    print('Reaction:', rxn)
                    if rxn:
                        return {'format': 'smiles', 'content': clip}
                except:
                    print('Not Reaction SMILES')

            if re.match(r'^\d+-\d+-\d+$', clip):
                return {'format': 'cas', 'content': clip}
            else:
                return {'format': 'iupac', 'content': clip}

    def clipboard_to(self, sender, output_format):
        clipboard_input = self.clipboard_content_identification()
        # Only do the check if the last clipboard was the same as previous if the called by a Timer
        if self.last_clipboard_cache == clipboard_input['content'] and isinstance(sender, rumps.Timer):
            print('same cache',
                  clipboard_input['content'], self.last_clipboard_cache)
            return
        else:
            self.last_clipboard_cache = clipboard_input['content']
            output, smiles = self.converter.convert(
                clipboard_input['content'], clipboard_input['format'], output_format)
            print('output_format', output_format,
                  "clipboard_input['format']", clipboard_input['format'], clipboard_input['content'])

            if output:
                if output_format in ['remove_atom_mapping', 'canonicalize', 'standardize', 'augment']:
                    pyperclip.copy(output)
                    rumps.notification(f"✅ {output_format.capitalize().replace('_', ' ')} copied to Clipboard for ",
                                       f"{clipboard_input['content']}",
                                       str(output), action_button=' ')

                elif output_format == 'price':
                    self.purchase_link = output['link']
                    rumps.notification(
                        f'💰 Buy {output["item_name"]}',
                        f'{output["amount"]} for {output["price"]} CHF',
                        f'or {output["price_per"]:.2f} CHF/g',
                        action_button='Open Link',
                        data='open_purchase_link',
                        icon='pictograms/carlos_helper_good.png')
                else:
                    if output_format == 'image':
                        self.copy_image_to_clipboard(output)
                    else:
                        pyperclip.copy(output)
                        self.last_clipboard_cache = output

                    self.smiles_to_history(smiles)
                    rumps.notification(f"✅ {output_format.capitalize()} copied to Clipboard for ",
                                       f"{clipboard_input['format'].capitalize()} : {clipboard_input['content']}",
                                       str(output), action_button=' ')

            else:
                rumps.notification(f"🚫 Couldn't convert Clipboard {output_format.capitalize()}  ",
                                   f"No {output_format.capitalize()} found for", str(
                                       pyperclip.paste()), sound=False, icon='pictograms/carlos_helper_good.png')
            return output

    def stop_clipboard_to_smiles_timer(self, sender):
        """
        Stops the clipboard-to-SMILES timer and updates the watcher button.

        Args:
        - sender: The sender object that triggered the event.
        """
        self.clipboard_to_smiles_watcher_button.title = 'Watch Clipboard for Molecules'
        self.clipboard_to_smiles_watcher_button.set_callback(self.start_clipboard_to_smiles_timer)
        self.clipboard_to_smiles_timer.stop()

    def start_clipboard_image_queue_timer(self, sender):
        self.clipboard_to_image_queue_button.title = '🛑 End Molecule Watching'
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
        image = ImageGrab.grabclipboard()
        if isinstance(image, Image.Image) and image != self.last_image_cache:
            # Save the image as a PNG file
            self.last_image_cache = image
            return image
        else:
            return False

    def save_image_from_clipboard(self, state):
        image = self.get_image_from_clipboard(state)
        if image:
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            file_name = f'{timestamp}_clipboard_image.png'
            file_path = self.image_input_dir / file_name
            image.save(file_path, "PNG")

    def copy_image_to_clipboard(self, image_path):
        try:
            temp_path = os.path.abspath(image_path)
            subprocess.run(
                ["osascript", "-e",
                 'set the clipboard to (read (POSIX file "' + str(temp_path) + '") as «class PNGf»)'])
        except:
            rumps.notification("🚫 Failed to load Image into clipboard",
                               '', str(image_path),
                               sound=False,
                               icon='pictograms/carlos_helper_bad.png')

    def remove_over_max_images(self, ):
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
        for filename in os.listdir(self.image_input_dir):
            if filename.endswith('.png'):
                image_path = str(self.image_input_dir / str(filename))
                smiles = self.converter.image_to_smiles(
                    image_path, image_or_path='path')
                if smiles:
                    self.smiles_to_history(smiles)

    def smiles_to_history(self, smiles):
        try:
            filename = self.converter.smiles_to_image(smiles)
            new_filename = os.path.basename(filename)
            timestamp = ''.join(str(new_filename).split("_")[0])
        except:
            new_filename = 'False'
            timestamp = 'False'

        # Check for duplicates if the duplication option is enabled
        rows_to_write = []
        with open(self.csv_file_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row and row[0] == smiles:
                    continue
                else:
                    rows_to_write.append(row)

        rows_to_write.append([smiles, new_filename, timestamp])

        with open(self.csv_file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(rows_to_write)

        self.update_history_menu()
        self.remove_over_max_images()

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
                smiles_list.append(
                    (row[smiles_col_index], row[image_file_col_index]))
        self.smiles_history = smiles_list[-self.length_history:]

    def show_history(self, sender=None):
        # Display the SMILES history
        def factory_copy_to_clip(value):
            def f(value):
                pyperclip.copy(value.title)

            return f

        def factory_structure(value, input, output_format):
            def f(value):
                if output_format in ['cas', 'selfies', 'inchi', 'iupac', 'mol']:
                    output, _ = self.converter.convert(
                        input, 'smiles', output_format)
                    if output:
                        pyperclip.copy(output)
                        rumps.notification(f"✅ {output_format.capitalize()} copied to Clipboard",
                                           f"",
                                           str(output))
                    else:
                        rumps.notification(f"🚫 Couldn't convert to {output_format.capitalize()}  ",
                                           f"No {output_format.capitalize()} found for", str(
                                               input), sound=False)
                elif output_format == 'smiles':
                    pyperclip.copy(input)
                    rumps.notification(f"✅ {output_format.capitalize()} copied to Clipboard",
                                       '',
                                       str(input))
                elif output_format == 'image':
                    self.copy_image_to_clipboard(image_path=input)
                    rumps.notification(
                        "✅ Image copied to Clipboard", '', f'{input}')

                elif output_format == 'price':
                    pyperclip.copy(input)
                    self.clipboard_to(input, output_format='price')
                else:
                    pyperclip.copy(input)

            return f

        for hist_smi, hist_img_file in reversed(self.smiles_history):
            build_function_smiles = factory_structure('', hist_smi, 'smiles')
            build_function_price = factory_structure('', hist_smi, 'price')
            build_function_inchi = factory_structure('', hist_smi, 'inchi')
            build_function_selfies = factory_structure('', hist_smi, 'selfies')
            build_function_mol = factory_structure('', hist_smi, 'mol')
            build_function_iupac = factory_structure('', hist_smi, 'iupac')
            build_function_cas = factory_structure('', hist_smi, 'cas')

            properties_dict = smiles_to_molecular_properties(hist_smi)
            if hist_img_file == 'False':
                self.history_item_smiles.add(rumps.MenuItem(
                    hist_smi, callback=build_function_smiles))

                self.history_item_price.add(rumps.MenuItem(
                    hist_smi, callback=build_function_price, ))
                temp_menu_structure = rumps.MenuItem(
                    hist_smi, icon=str(self.image_generated_dir / hist_img_file), dimensions=[75, 75])
            else:
                build_function_image = factory_structure(
                    '', str(self.image_generated_dir / hist_img_file), 'image')
                temp_menu_structure = rumps.MenuItem(
                    '', icon=str(self.image_generated_dir / hist_img_file), dimensions=[75, 75])
                temp_menu_structure.add(rumps.MenuItem(
                    "Image", callback=build_function_image))

                self.history_item_smiles.add(rumps.MenuItem(
                    '', callback=build_function_smiles, icon=str(self.image_generated_dir / hist_img_file),
                    dimensions=[75, 75]))
                self.history_item_price.add(rumps.MenuItem(
                    '', callback=build_function_price, icon=str(self.image_generated_dir / hist_img_file),
                    dimensions=[75, 75]))

            temp_menu_structure.add(rumps.MenuItem(
                "SMILES", callback=build_function_smiles))
            temp_menu_structure.add(rumps.MenuItem(
                "CAS", callback=build_function_cas))
            temp_menu_structure.add(rumps.MenuItem(
                "IUPAC", callback=build_function_iupac))
            temp_menu_structure.add(rumps.MenuItem(
                "InChI", callback=build_function_inchi))
            temp_menu_structure.add(rumps.MenuItem(
                "RDkit Mol", callback=build_function_mol))
            temp_menu_structure.add(rumps.MenuItem(
                "Selfie", callback=build_function_selfies))

            self.history_structure_item.add(temp_menu_structure)
            if properties_dict:
                if hist_img_file == 'False':
                    temp_menu_properties = rumps.MenuItem(
                        hist_smi, callback=build_function_smiles)
                else:
                    temp_menu_properties = rumps.MenuItem(
                        '', callback=build_function_smiles,
                        icon=str(self.image_generated_dir / hist_img_file), dimensions=[75, 75])
                for key in properties_dict:
                    left_aligned = f"{key:<{18}}"
                    temp = f'{properties_dict[key]:.3f}'
                    right_aligned = f"{temp:>{5}}"
                    build_function_properties = factory_copy_to_clip(
                        properties_dict[key])
                    temp_menu_properties.add(rumps.MenuItem(
                        left_aligned + right_aligned, callback=build_function_properties))
            self.history_molecular_properties_item.add(temp_menu_properties)

    def update_history_menu(self, sender=None):
        # Update the history menu with the last 10 SMILES strings
        self.history_item_smiles.clear()
        self.history_structure_item.clear()
        self.history_molecular_properties_item.clear()
        self.history_item_price.clear()
        self.read_history()
        self.show_history()

    def open_history_file(self, sender):
        os.system(f"open {self.csv_file_path}")

    def open_queue_folder(self, sender):
        os.system(f"open {self.image_input_dir}")


if __name__ == '__main__':
    app = ClipboardImageSaverApp()
    app.app.run()
