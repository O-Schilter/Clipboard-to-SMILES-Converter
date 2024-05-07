import PyInstaller.__main__


PyInstaller.__main__.run([
    'clipboard2smiles_qt.py',
    '--windowed',
    '--icon=/Users/oli/projects/carlos_helper/pictograms/carlos_helper_logo.icns',
    '--add-data=smiles.csv:.',
    '--add-data=pictograms:pictograms',
    '--add-data=image_generated:image_generated',
    '--add-data=image_output:image_output',
    '--add-data=image_queue:image_queue',
    '--add-data=molscribe:molscribe',
    '--add-data=models--yujieq--MolScribe:models--yujieq--MolScribe',
    '--collect-submodules=molscribe'
    '--hidden-import=molscribe',
    '--paths=/Users/oli/miniforge3/envs/clip2qt',
    '--noconfirm'
])