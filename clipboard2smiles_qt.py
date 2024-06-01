import sys
import csv
import os
import re
import subprocess
from functools import partial
from pathlib import Path
from datetime import datetime
import torch
import pyperclip
import requests
from PIL import ImageGrab, Image

import selfies as sf
from rdkit import Chem
from rdkit.Chem import AllChem

from molscribe.interface import MolScribe

from utils import smiles_to_molecular_properties
from converter import Converter
import darkdetect
from PySide6.QtCore import QSize, Qt, QTimer, Signal, QDir
from PySide6.QtGui import QAction, QIcon
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QMenu,
    QSystemTrayIcon,
    QProxyStyle,
    QStyle,
)

basedir = os.path.dirname(__file__)


class MyProxyStyle(QProxyStyle):
    pass

    def pixelMetric(self, QStyle_PixelMetric, option=None, widget=None):
        if QStyle_PixelMetric == QStyle.PM_SmallIconSize:
            return 80
        else:
            return QProxyStyle.pixelMetric(self, QStyle_PixelMetric, option, widget)


class ClipboardImageSaverApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Define app settings
        self.confidence_level = 0.66
        self.max_n_images = 25
        self.length_history = 8

        self.model_path = os.path.join(
            basedir,
            "models--yujieq--MolScribe",
            "snapshots",
            "601bb0f491fb9598fa40227ae9759f800661cd74",
            "swin_base_char_aux_1m.pth",
        )
        self.model_device = "cpu"

        # Define image input/output directories
        self.image_input_dir = os.path.join(
            basedir, "image_queue"
        )  # Path("image_queue/")
        self.image_output_dir = os.path.join(
            basedir, "image_output"
        )  # Path("image_output/")
        self.image_generated_dir = os.path.join(basedir, "image_generated")

        # Define history file and initialize empty history
        self.csv_file_path = os.path.join(basedir, "smiles.csv")
        self.smiles_history = []
        self.purchase_link = ""
        self.internet_connection = False

        # Create a MolScribe model instance
        self.model = MolScribe(self.model_path, device=torch.device(self.model_device))

        # Create Converter class
        self.converter = Converter(
            model=self.model,
            image_input_dir=self.image_input_dir,
            image_generated_dir=self.image_generated_dir,
        )

        # Create menu and actions
        self.create_menu()

        # Create empty Timer
        self.clipboard_to_smiles_timer = QTimer(self)

        # Image cache, to only run conversion when a new Image is detected
        self.last_image_cache = None
        self.last_clipboard_cache = None

        # Initialize the SMILES history and load the previously saved molecules
        self.read_history()
        self.show_history()

        # Setup system tray icon
        self.tray_icon = QSystemTrayIcon(self)
        if darkdetect.isDark():
            self.pictogram_folder = os.path.join(basedir, "pictograms", "dark")
        else:
            self.pictogram_folder = os.path.join(basedir, "pictograms", "light")
            
        self.tray_icon.Information = QIcon(
            os.path.join(self.pictogram_folder, "carlos_helper_logo.png")
        )
        self.tray_icon.Warning = QIcon(
            os.path.join(self.pictogram_folder, "carlos_helper_bad.png")
        )
        self.tray_icon.Critical = QIcon(
            os.path.join(self.pictogram_folder, "carlos_helper_bad.png")
        )

        self.normal_status_icon = QIcon(
            os.path.join(self.pictogram_folder, "carlos_helper_logo.png")
        )
        self.recording_status_icon = QIcon(
            os.path.join(self.pictogram_folder, "carlos_helper_recording.png")
        )

        self.tray_icon.setIcon(
            QIcon(os.path.join(self.pictogram_folder, "carlos_helper_good.png"))
        )
        self.tray_icon.setContextMenu(self.tray_menu)
        self.tray_icon.show()

    def create_menu(self):
        # Create menu and actions
        self.menu = QMenu("test")

        # Chemical structure conversion actions
        self.to_smiles_action = QAction("Convert Clipboard to SMILES", self)
        self.to_smiles_action.triggered.connect(
            partial(self.clipboard_to, output_format="smiles")
        )
        self.menu.addAction(self.to_smiles_action)

        self.to_image_action = QAction("Convert Clipboard to Image", self)
        self.to_image_action.triggered.connect(
            partial(self.clipboard_to, output_format="image")
        )
        self.menu.addAction(self.to_image_action)

        self.to_cas_action = QAction("Convert Clipboard to CAS", self)
        self.to_cas_action.triggered.connect(
            partial(self.clipboard_to, output_format="cas")
        )
        self.menu.addAction(self.to_cas_action)

        self.to_iupac_action = QAction("Convert Clipboard to IUPAC", self)
        self.to_iupac_action.triggered.connect(
            partial(self.clipboard_to, output_format="iupac")
        )
        self.menu.addAction(self.to_iupac_action)

        self.to_inchi_action = QAction("Convert Clipboard to InChI", self)
        self.to_inchi_action.triggered.connect(
            partial(self.clipboard_to, output_format="inchi")
        )
        self.menu.addAction(self.to_inchi_action)

        self.to_mol_action = QAction("Convert Clipboard to RDkit Mol", self)
        self.to_mol_action.triggered.connect(
            partial(self.clipboard_to, output_format="mol")
        )
        self.menu.addAction(self.to_mol_action)

        self.to_selfies_action = QAction("Convert Clipboard to Selfies", self)
        self.to_selfies_action.triggered.connect(
            partial(self.clipboard_to, output_format="selfies")
        )
        self.menu.addAction(self.to_selfies_action)

        self.menu.addSeparator()

        self.smiles_to_smiles_menu = QMenu("Clipboard SMILES Operations", self)

        canonicalize_action = QAction("Canonicalize SMILES", self)
        canonicalize_action.triggered.connect(
            partial(self.clipboard_to, output_format="canonicalize")
        )
        self.smiles_to_smiles_menu.addAction(canonicalize_action)

        augment_action = QAction("Augment SMILES", self)
        augment_action.triggered.connect(partial(print, "augment"))
        self.smiles_to_smiles_menu.addAction(augment_action)

        remove_atom_mapping_action = QAction("Remove Atom Mapping from SMILES", self)
        remove_atom_mapping_action.triggered.connect(
            partial(self.clipboard_to, output_format="remove_atom_mapping")
        )
        self.smiles_to_smiles_menu.addAction(remove_atom_mapping_action)

        self.menu.addMenu(self.smiles_to_smiles_menu)
        self.menu.addSeparator()

        # Find Price for molecules
        self.to_price_action = QAction("Convert Clipboard to Price", self)
        self.to_price_action.triggered.connect(
            partial(self.clipboard_to, output_format="price")
        )
        self.menu.addAction(self.to_price_action)

        self.history_item_price_menu = QMenu("Find Price from History", self)
        self.menu.addMenu(self.history_item_price_menu)

        self.menu.addSeparator()

        # Clipboard monitoring
        self.description_detection_action = QAction(
            "Continuous Molecule Detection", self
        )
        self.description_detection_action.setEnabled(False)
        self.menu.addAction(self.description_detection_action)
        self.clipboard_to_smiles_watcher_action = QAction(
            "Watch Clipboard for Molecules", self
        )
        self.clipboard_to_smiles_watcher_action.triggered.connect(
            self.start_clipboard_to_smiles_timer
        )
        self.menu.addAction(self.clipboard_to_smiles_watcher_action)

        self.menu.addSeparator()

        # History menu
        self.history_item_smiles_menu = QMenu("Copy SMILES from History", self)
        self.menu.addMenu(self.history_item_smiles_menu)

        self.history_structure_menu = QMenu("Copy Structure from History", self)
        self.menu.addMenu(self.history_structure_menu)

        self.history_molecular_properties_menu = QMenu(
            "Molecular Properties from History", self
        )
        self.menu.addMenu(self.history_molecular_properties_menu)

        self.open_history_file_action = QAction("Open History File", self)
        self.open_history_file_action.triggered.connect(partial(self.open_history_file))
        self.menu.addAction(self.open_history_file_action)

        self.menu.addSeparator()
        # Continuous monitoring queue
        self.description_queue_action = QAction(
            "Monitor Clipboard to Build Image Queue", self
        )
        self.description_queue_action.setEnabled(False)
        self.menu.addAction(self.description_queue_action)
        self.menu.addAction(self.description_queue_action)

        self.clipboard_to_image_queue_action = QAction(
            "Start Clipboard Monitoring", self
        )
        self.clipboard_to_image_queue_action.triggered.connect(
            self.start_clipboard_image_queue_timer
        )
        self.menu.addAction(self.clipboard_to_image_queue_action)

        self.open_queue_folder_action = QAction("Add Images to Queue", self)
        self.open_queue_folder_action.triggered.connect(self.open_queue_folder)
        self.menu.addAction(self.open_queue_folder_action)

        self.process_queue_folder_action = QAction("Process Image Queue", self)
        self.process_queue_folder_action.triggered.connect(self.batch_images_to_smiles)
        self.menu.addAction(self.process_queue_folder_action)

        self.menu.addSeparator()

        # Options
        self.options_menu = QMenu("Options", self)
        self.menu.addMenu(self.options_menu)

        # Select Vendor Menu
        self.select_vendor_menu = QMenu("Select Vendor", self)
        self.options_menu.addMenu(self.select_vendor_menu)

        self.enamine_action = QAction("Enamine", self, checkable=True)
        self.enamine_action.triggered.connect(partial(self.select_vendor, "Enamine"))
        self.select_vendor_menu.addAction(self.enamine_action)
        self.enamine_action.setChecked(False)

        self.chemie_brunschwieg_action = QAction(
            "Chemie Brunschwieg", self, checkable=True
        )
        self.chemie_brunschwieg_action.triggered.connect(
            partial(self.select_vendor, "ChemieBrunschwieg")
        )
        self.chemie_brunschwieg_action.setChecked(True)
        self.select_vendor_menu.addAction(self.chemie_brunschwieg_action)

        # Credits
        share_message = "I Just discovered an amazing app that lets you convert images to SMILES and other chemical identifiers directly in the clipboard 🤯. Say goodbye to manual drawing molecules and hello to seamless chemistry workflows. Check it out now https://github.com/O-Schilter/Clipboard-to-SMILES-Converter Credits to @OSchilter"
        encoded_share_message = requests.utils.quote(share_message, safe="")

        share_url = f"https://twitter.com/intent/tweet?text={encoded_share_message}"
        self.share_action = QAction("Share the App", self)
        self.share_action.triggered.connect(lambda: self.open_url(share_url))
        self.options_menu.addAction(self.share_action)

        credits_url = "https://github.com/O-Schilter/Clipboard-to-SMILES-Converter"
        self.credits_action = QAction("Credits & Updates", self)
        self.credits_action.triggered.connect(lambda: self.open_url(credits_url))
        self.options_menu.addAction(self.credits_action)

        self.quit_action = QAction("Quit", self)
        self.quit_action.triggered.connect(QApplication.instance().quit)
        self.menu.addAction(self.quit_action)

        # Set up the system tray menu
        self.tray_menu = self.menu

    def read_history(self):
        with open(self.csv_file_path, "r") as csvfile:
            reader = csv.reader(csvfile)
            # Get the header row and find the index of the "smiles" column
            header_row = next(reader)
            smiles_col_index = header_row.index("smiles")
            image_file_col_index = header_row.index(" image_file")
            # Loop through the rows of the CSV file and select all entries from the "smiles" column
            smiles_list = []
            for row in reader:
                smiles_list.append((row[smiles_col_index], row[image_file_col_index]))
        self.smiles_history = smiles_list[-self.length_history :]

    def show_history(self):
        def factory_copy_to_clip(value):
            def copy_to_clipboard():
                pyperclip.copy(value)

            return copy_to_clipboard

        def factory_structure(value, input_format, output_format):
            def convert_and_copy():
                if output_format in ["cas", "selfies", "inchi", "iupac", "mol"]:
                    output, _ = self.converter.convert(
                        input_format, "smiles", output_format
                    )
                    if output:
                        pyperclip.copy(output)
                        self.tray_icon.showMessage(
                            f"✅ {output_format.capitalize()} copied to Clipboard",
                            str(output),
                            self.tray_icon.Information,
                            5000,
                        )
                    else:
                        self.tray_icon.showMessage(
                            f"🚫 Couldn't convert to {output_format.capitalize()}  ",
                            f"No {output_format.capitalize()} found for \n "
                            + str(input),
                            self.tray_icon.Critical,
                            5000,
                        )

                elif output_format == "smiles":
                    pyperclip.copy(input_format)
                    self.tray_icon.showMessage(
                        f"✅ {output_format.capitalize()} copied to Clipboard",
                        str(input_format),
                        self.tray_icon.Information,
                        5000,
                    )

                elif output_format == "image":
                    output, _ = self.converter.convert(
                        input_format, "smiles", output_format
                    )
                    self.copy_image_to_clipboard(image_path=output)
                    self.tray_icon.showMessage(
                        f"✅ {output_format.capitalize()} copied to Clipboard",
                        str(output),
                        self.tray_icon.Information,
                        5000,
                    )
                elif output_format == "price":
                    output, _ = self.converter.convert(
                        input_format, "smiles", output_format
                    )
                    self.tray_icon.showMessage(
                        f"💰 Buy {output['item_name']}".format(output=output),
                        f"{output['amount']} for {output['price']} CHF\nor {output['price_per']:.2f} CHF/g".format(
                            output=output
                        ),
                        self.tray_icon.Information,
                        10000,
                    )
                    self.tray_icon.messageClicked.connect(self.open_url(output["link"]))

                else:
                    pyperclip.copy(input_format)

            return convert_and_copy

        # Clear existing menus
        self.history_item_smiles_menu.clear()
        self.history_structure_menu.clear()
        self.history_molecular_properties_menu.clear()
        self.history_item_price_menu.clear()

        for hist_smi, hist_img_file in reversed(self.smiles_history):
            build_function_smiles = factory_structure(hist_smi, hist_smi, "smiles")
            build_function_price = factory_structure(hist_smi, hist_smi, "price")
            build_function_inchi = factory_structure(hist_smi, hist_smi, "inchi")
            build_function_selfies = factory_structure(hist_smi, hist_smi, "selfies")
            build_function_mol = factory_structure(hist_smi, hist_smi, "mol")
            build_function_iupac = factory_structure(hist_smi, hist_smi, "iupac")
            build_function_cas = factory_structure(hist_smi, hist_smi, "cas")

            properties_dict = smiles_to_molecular_properties(hist_smi)
            if hist_img_file == "False":
                self.history_item_smiles_menu.addAction(
                    QAction(hist_smi, self, triggered=build_function_smiles)
                )
                self.history_item_price_menu.addAction(
                    QAction(hist_smi, self, triggered=build_function_price)
                )
                temp_menu_structure = QMenu(hist_smi, self)
            else:
                build_function_image = factory_structure(
                    os.path.join(self.image_generated_dir, str(hist_img_file)),
                    hist_smi,
                    "image",
                )
                temp_menu_structure = QMenu(self)
                temp_menu_structure.setIcon(
                    QIcon(
                        str(Path(os.path.join(self.image_generated_dir, hist_img_file)))
                    )
                )
                temp_menu_structure.addAction(
                    QAction("Image", self, triggered=build_function_image)
                )

                self.history_item_smiles_menu.addAction(
                    QAction(
                        self,
                        icon=QIcon(
                            str(os.path.join(self.image_generated_dir, hist_img_file))
                        ),
                        triggered=build_function_smiles,
                    )
                )
                self.history_item_price_menu.addAction(
                    QAction(
                        self,
                        icon=QIcon(
                            str(os.path.join(self.image_generated_dir, hist_img_file))
                        ),
                        triggered=build_function_price,
                    )
                )

            temp_menu_structure.addAction(
                QAction("SMILES", self, triggered=build_function_smiles)
            )
            temp_menu_structure.addAction(
                QAction("CAS", self, triggered=build_function_cas)
            )
            temp_menu_structure.addAction(
                QAction("IUPAC", self, triggered=build_function_iupac)
            )
            temp_menu_structure.addAction(
                QAction("InChI", self, triggered=build_function_inchi)
            )
            temp_menu_structure.addAction(
                QAction("RDkit Mol", self, triggered=build_function_mol)
            )
            temp_menu_structure.addAction(
                QAction("Selfie", self, triggered=build_function_selfies)
            )

            self.history_structure_menu.addMenu(temp_menu_structure)

            if properties_dict:
                if hist_img_file == "False":
                    temp_menu_properties = QMenu(hist_smi, self)
                else:
                    temp_menu_properties = QMenu(self)
                    temp_menu_properties.setIcon(
                        QIcon(
                            str(os.path.join(self.image_generated_dir, hist_img_file))
                        )
                    )
                for key in properties_dict:
                    left_aligned = f"{key:<{18}}"
                    temp = f"{properties_dict[key]:.3f}"
                    right_aligned = f"{temp:>{5}}"
                    build_function_properties = factory_copy_to_clip(
                        properties_dict[key]
                    )
                    temp_menu_properties.addAction(
                        QAction(
                            left_aligned + right_aligned,
                            self,
                            triggered=build_function_properties,
                        )
                    )
            self.history_molecular_properties_menu.addMenu(temp_menu_properties)

    def remove_over_max_images(
        self,
    ):
        # delete old images if more than self.max_n_images:
        all_files = sorted(Path(self.image_output_dir).glob("*.png"))
        if len(all_files) > self.max_n_images:
            for old_file in all_files[: -self.max_n_images :]:
                os.remove(old_file)
        all_files = sorted(Path(self.image_generated_dir).glob("*.png"))
        if len(all_files) > self.max_n_images:
            for old_file in all_files[: -self.max_n_images :]:
                os.remove(old_file)

    def smiles_to_history(self, smiles):
        try:
            filename = self.converter.smiles_to_image(smiles)
            new_filename = os.path.basename(filename)
            timestamp = "".join(str(new_filename).split("_")[0])
        except:
            new_filename = "False"
            timestamp = "False"

        # Check for duplicates if the duplication option is enabled
        rows_to_write = []
        with open(self.csv_file_path, "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row and row[0] == smiles:
                    continue
                else:
                    rows_to_write.append(row)

        rows_to_write.append([smiles, new_filename, timestamp])

        with open(self.csv_file_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(rows_to_write)

        self.update_history_menu()
        self.remove_over_max_images()

    def open_url(self, url):
        if os.name == "nt":  # Windows
            os.system(f"start {url}")
        elif os.name == "posix":
            os.system(f'open "{url}"')

    def select_vendor(self, vendor):
        self.converter.vendors.select_vendor(vendor)
        for action in self.select_vendor_menu.actions():
            action.setChecked(action.text().replace(" ", "") == vendor)

    def start_clipboard_to_smiles_timer(self):
        self.clipboard_to_smiles_watcher_action.setText("Stop Watching Clipboard")
        self.clipboard_to_smiles_watcher_action.triggered.disconnect()
        self.clipboard_to_smiles_watcher_action.triggered.connect(
            self.stop_clipboard_to_smiles_timer
        )

        self.tray_icon.setIcon(self.recording_status_icon)

        self.clipboard_to_smiles_timer.timeout.connect(
            partial(self.clipboard_to, output_format="smiles")
        )
        self.clipboard_to_smiles_timer.start(500)

    def stop_clipboard_to_smiles_timer(self):
        self.clipboard_to_smiles_watcher_action.setText("Watch Clipboard for Molecules")

        self.clipboard_to_smiles_watcher_action.triggered.disconnect()
        self.clipboard_to_smiles_watcher_action.triggered.connect(
            self.start_clipboard_to_smiles_timer
        )

        self.tray_icon.setIcon(self.normal_status_icon)

        self.clipboard_to_smiles_timer.stop()

    def start_clipboard_image_queue_timer(self):
        self.clipboard_to_image_queue_action.setText("Stop Clipboard Monitoring")
        self.clipboard_to_image_queue_action.triggered.disconnect()
        self.clipboard_to_image_queue_action.triggered.connect(
            self.stop_clipboard_image_queue_timer
        )

        self.tray_icon.setIcon(self.recording_status_icon)

        # self.clipboard_to_image_timer = QTimer(self)
        self.clipboard_to_smiles_timer.timeout.connect(self.save_image_from_clipboard)
        self.clipboard_to_smiles_timer.start(500)

    def stop_clipboard_image_queue_timer(self):
        self.clipboard_to_image_queue_action.setText("Start Clipboard Monitoring")
        self.clipboard_to_image_queue_action.triggered.disconnect()
        self.clipboard_to_image_queue_action.triggered.connect(
            self.start_clipboard_image_queue_timer
        )

        self.tray_icon.setIcon(self.normal_status_icon)

        self.clipboard_to_smiles_timer.stop()
        self.batch_images_to_smiles()

    def copy_image_to_clipboard(self, image_path):
        try:
            temp_path = os.path.abspath(image_path)
            with open(temp_path, "rb") as image_file:
                image_data = image_file.read()

                if os.name == "nt":  # Windows
                    import ctypes

                    open_clip = ctypes.windll.user32.OpenClipboard
                    empty_clip = ctypes.windll.user32.EmptyClipboard
                    set_clip = ctypes.windll.user32.SetClipboardData
                    close_clip = ctypes.windll.user32.CloseClipboard

                    open_clip(None)
                    empty_clip()

                    set_clip(ctypes.c_uint(8), ctypes.c_char_p(image_data))
                    close_clip()

                elif os.name == "posix":
                    subprocess.run(
                        [
                            "osascript",
                            "-e",
                            'set the clipboard to (read (POSIX file "'
                            + str(temp_path)
                            + '") as «class PNGf»)',
                        ]
                    )

        except:
            self.tray_icon.showMessage(
                "{🚫 Failed to load Image into clipboard",
                str(image_path),
                self.tray_icon.Critical,
                5000,
            )

    def update_history_menu(self, sender=None):
        # Update the history menu with the last 10 SMILES strings
        self.read_history()
        self.show_history()

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
            return {"format": "image", "content": clip}

        clip = pyperclip.paste()

        if clip:
            try:
                mol = Chem.MolFromSmiles(clip)
                if mol:
                    return {"format": "smiles", "content": clip}
            except:
                print("Not SMILES")

            if clip.startswith("InChI="):
                try:
                    mol = Chem.MolFromInchi(clip)
                    if mol is not None:
                        return {"format": "inchi", "content": clip}
                except:
                    print("Invalid InChI")

            try:
                mol = Chem.MolFromMolBlock(clip)
                if mol is not None:
                    return {"format": "mol", "content": clip}
            except:
                print("Invalid MOL file")

            try:
                smi = sf.decoder(clip)
                if smi is not None and smi != "" and len(smi) > 0:
                    print("Selfies converted to SMILES:", smi)
                    return {"format": "selfies", "content": clip}
            except:
                print("Invalid Selfies")

            if ">" in clip:
                print('">" detected')
                try:
                    rxn = AllChem.ReactionFromSmarts(clip)
                    print("Reaction:", rxn)
                    if rxn:
                        return {"format": "smiles", "content": clip}
                except:
                    print("Not Reaction SMILES")

            if re.match(r"^\d+-\d+-\d+$", clip):
                return {"format": "cas", "content": clip}
            else:
                return {"format": "iupac", "content": clip}

    def clipboard_to(self, output_format):
        clipboard_input = self.clipboard_content_identification()

        # Only do the check if the last clipboard was the same as previous if the called by a Timer
        if (
            self.last_clipboard_cache == clipboard_input["content"]
            and isinstance(self.clipboard_to_smiles_timer, QTimer)
            and self.clipboard_to_smiles_timer.isActive()
        ):
            print("same cache", clipboard_input["content"], self.last_clipboard_cache)
            return
        else:
            self.last_clipboard_cache = clipboard_input["content"]
            output, smiles = self.converter.convert(
                clipboard_input["content"], clipboard_input["format"], output_format
            )

            if output:
                if output_format in [
                    "remove_atom_mapping",
                    "canonicalize",
                    "standardize",
                    "augment",
                ]:
                    pyperclip.copy(output)

                    self.tray_icon.showMessage(
                        f"✅ {output_format.capitalize()} copied to Clipboard for ",
                        f"{clipboard_input['content']}\n{output}",
                        self.tray_icon.Information,
                        5000,
                    )

                elif output_format == "price":
                    self.tray_icon.showMessage(
                        f"💰 Buy {output['item_name']}".format(output=output),
                        f"{output['amount']} for {output['price']} CHF\nor {output['price_per']:.2f} CHF/g".format(
                            output=output
                        ),
                        self.tray_icon.Information,
                        10000,
                    )

                    self.tray_icon.messageClicked.connect(self.open_url(output["link"]))

                else:
                    if output_format == "image":
                        self.copy_image_to_clipboard(output)
                    else:
                        pyperclip.copy(output)
                        self.last_clipboard_cache = output

                    self.smiles_to_history(smiles)
                    self.tray_icon.showMessage(
                        f"✅ {output_format.capitalize()} copied to Clipboard",
                        str(output),
                        self.tray_icon.Information,
                        5000,
                    )

            else:
                self.tray_icon.showMessage(
                    "Conversion Failed",
                    f"Couldn't convert Clipboard {output_format.capitalize()}. No {output_format.capitalize()} found for {pyperclip.paste()}",
                    self.tray_icon.Warning,
                    5000,
                )
                return None

        return output

    # Rest of the code stays the same
    def open_history_file(self):
        if os.name == "nt":  # Windows
            subprocess.run(["start", self.csv_file_path], shell=True)
        else:
            subprocess.run(["open", self.csv_file_path])

    def open_queue_folder(self):
        if os.name == "nt":  # Windows
            subprocess.run(["explorer", self.image_input_dir], shell=True)
        else:
            subprocess.run(["open", self.image_input_dir])

    def batch_images_to_smiles(self):
        for filename in os.listdir(self.image_input_dir):
            if filename.endswith(".png"):
                image_path = os.path.join(str(self.image_input_dir), str(filename))
                smiles = self.converter.image_to_smiles(
                    image_path, image_or_path="path"
                )
                if smiles:
                    self.smiles_to_history(smiles)

    def get_image_from_clipboard(self):
        image = ImageGrab.grabclipboard()
        if isinstance(image, Image.Image) and image != self.last_image_cache:
            # Save the image as a PNG file
            self.last_image_cache = image
            return image
        else:
            return False

    def save_image_from_clipboard(self):
        image = self.get_image_from_clipboard()
        if image:
            timestamp = datetime.now().strftime("%Y-%m-%d %H-%M-%S")
            file_name = f"{timestamp}_clipboard_image.png"
            file_path = os.path.join(self.image_input_dir, file_name)
            image.save(file_path, "PNG")


if __name__ == "__main__":
    myStyle = MyProxyStyle(
        "Fusion"
    )  # The proxy style should be based on an existing style,
    app = QApplication(sys.argv)
    app.setStyle(myStyle)

    app.setQuitOnLastWindowClosed(False)
    clipboard_app = ClipboardImageSaverApp()
    sys.exit(app.exec())
