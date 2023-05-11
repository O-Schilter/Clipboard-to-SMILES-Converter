from setuptools import setup
import sys
sys.setrecursionlimit(5000)

APP = ['clipboard2smiles.py']
DATA_FILES = ['smiles.csv','pictograms','image_generated','image_output','image_queue','models--yujieq--MolScribe']
OPTIONS = {
    'argv_emulation': True,
    'plist': {
        'LSUIElement': True,
    },
    'packages': ['rumps','chardet'],
    'iconfile': '/Users/oli/projects/carlos_helper/pictograms/carlos_helper_logo.icns'}

setup(
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
) 