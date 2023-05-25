# Clipboard-to-SMILES-Converter
Clipboard-to-SMILES-Converter is a lightweight macOS app designed to simplify the conversion process between images and SMILES notations for chemical structures. It offers a user-friendly interface, featuring a built-in history for easy reference, batch processing capabilities for efficient conversion of multiple images, continuous clipboard monitoring to automatically convert detected images, and all of this without requiring an internet connection.

![](paper/Clipboard2SMILES.png)
## Installation as App

1. Download the zip file:
 - Visit the provided Box link and download the zip file containing Clipboard-to-SMILES-Converter.
2. Unzip the file:
 - Locate the downloaded zip file and double-click on it to unzip its contents.
3. Execute the app:
 - After unzipping, you will find the “carlos_helper.app” file. Double-click on it to execute the app.
4. App not authorized popup:
 - When you try to run the app, a macOS popup might appear, indicating that the app is not authorized to run.
5. Authorize Clipboard-to-SMILES-Converter in Security settings:
 - Open the “System Preferences” by clicking on the Apple menu in the top-left corner of the screen and selecting “System Preferences”.
 - In the System Preferences window, click on “Security & Privacy”.
 - In the Security & Privacy settings, scroll down.
 - Allow to open carlos_helper.
6. (Optional) Change keyboard shortcut for screenshots:
 - Still, in the System Preferences window, click on “Keyboard”.
 - In the Keyboard settings, select the “Shortcuts” tab.
 - In the left sidebar, choose “Screen Shots” from the list.
 - Click on the existing keyboard shortcut for “Copy screenshot of a selected area to clipboard“ and assign a new shortcut, such as “Cmd + 4”, which will be more convenient since Clipboard-to-SMILES-Converter uses the default “Cmd + Shift +Option+ 4" shortcut for image conversion.
7. Relaunch Clipboard-to-SMILES-Converter:
 - Close the Security & Privacy and Keyboard settings.
 - Quit Clipboard-to-SMILES-Converter if it’s running.
 - Launch Clipboard-to-SMILES-Converter again to ensure that the changes and authorizations have taken effect, the startup can take up to a minute since it loads the model. (maybe a second restart is required) 
8. (Optional but recommended) Allow notifications:
 - When Clipboard-to-SMILES-Converter launches, you may be prompted to allow notifications. It is recommended to enable notifications since it is the main way the app communicates about the status of the conversion and other information from the app.


## Installation from Git

To use this project, please follow the steps below:

### Create and Activate Conda Environment
Create a new Conda environment using the following command:

```
conda create --name clipboard2smiles python=3.10
```
Activate the environment:
```
conda activate clipboard2smiles
```
git clone the repo in the desired path.
### Install Dependencies
Install the required dependencies by running the following command:
```
pip install -r requirements.txt 
```

To get more details about the underlaying model architecture refer to [MolScribe](https://github.com/thomas0809/MolScribe)

You can download and copy the required model with the following commands:

First you need to install huggingface_hub by running 
```
pip install huggingface_hub
```

then download the model with the following command

```
python - <<'EOF'
from huggingface_hub import hf_hub_download
import shutil

output_ckp_path = hf_hub_download('yujieq/MolScribe', 'swin_base_char_aux_1m.pth')
output_folder_path = '/'.join(output_ckp_path.split('/')[:-3])
shutil.copytree(output_folder_path, './models--yujieq--MolScribe')
EOF
```



### Lunching the App directly 
You can run the application directly using the following command:

```
python clipboard2smiles.py
```

### Build the App
To build the application, use the following command:

```
python setup.py py2app
```
This will create the necessary executable or bundle based on your operating system.
You may need to change the iconfile paramters in the setup.py file. If you are compiling the app on an Apple Silicon (M1, M2 etc.) see the FAQ of py2app to get the instruction to compile it for x86 architecture.


## Functionality Overview
A brief list of what each menu item is doing:
- Convert Clipboard Image to SMILES: Try to convert the image provided in the clipboard to SMILES. If their clipboard content is not a SMILES returns a notification with the current content of the clipboard, if it is an image it will save back the extracted molecular structure in the form of SMILES into the clipboard.
- Convert Clipboard SMILES to Image: If a valid SMILES is provided in the clipboard, the corresponding molecular structure is drawn and the image is directly saved back into the clipboard.
- Canonicalize Clipboard SMILES: If a valid SMILES is provided in the clipboard, the canonicalized form of the SMILES is saved back into the clipboard.
- Watch Clipboard for Molecules: The same functionality as "Convert Clipboard Image to SMILES" but done continuously without user input. Every image saved to the clipboard will be automatically scanned for a molecular structure, if one is detected with a high certainty the corresponding SMILES is saved back into the clipboard. If pressed the menu item changes into a stop button.
- Open History File: All converted SMILES are saved into History File with their corresponding timestamp.
- Copy SMILES from History: The ten most recent SMILES are displayed, when clicked on one the SMILES is saved into the clipboard.
- Copy Image from History: The ten most recent SMILES are displayed, when clicked on one the image of the molecule is saved into the clipboard.
- Start Clipboard monitoring: All images taken in the clipboard are saved into a queue folder. As soon as the user stops the clipboard monitoring, all images in the queue folder will be converted into SMILES. The main difference between "Watch Clipboard for Molecules" is that the conversion is done in batches and afterward. 
- Add images to queue: opens queue folder. In this folder the user can manually drag and drop, images to he/she/they would like to convert.
- Process image queue: starts converting all images in the queue folder.
