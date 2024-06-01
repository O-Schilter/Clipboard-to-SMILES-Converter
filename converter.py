import os
import torch
# import rumps
import chemdraw
from datetime import datetime


import selfies as sf
from rdkit import Chem
from smiles_operations import augemnt_smiles, canonicalize_smiles, remove_atom_mapping

from utils import create_filename_from_smiles
from api_conversions import smiles_to_cas_api, smiles_to_iupac_api, iupac_to_smiles_api, cas_to_smiles_api
from vendors import Vendors, Enamine, ChemieBrunschwieg


class Converter():
    def __init__(self, model, image_input_dir, image_generated_dir):
        self.model = model
        self.image_input_dir = image_input_dir
        self.image_generated_dir = image_generated_dir
        self.image_confidence_level = 0.66

        self.vendors = Vendors([ChemieBrunschwieg(), Enamine()])

        self.conversion_to_smiles = {
            'cas': self.cas_to_smiles,
            'iupac': self.iupac_to_smiles,
            'image': self.image_to_smiles,
            'smiles': self.smiles_to_smiles,
            'inchi': self.inchi_to_smiles,
            'mol': self.mol_to_smiles,
            'selfies': self.selfies_to_smiles,
        }
        self.conversion_from_smiles = {
            'cas': self.smiles_to_cas,
            'iupac': self.smiles_to_iupac,
            'image': self.smiles_to_image,
            'smiles': self.smiles_to_smiles,
            'price': self.smiles_to_price,
            'inchi': self.smiles_to_inchi,
            'mol': self.smiles_to_mol,
            'selfies': self.smiles_to_selfies,
            'remove_atom_mapping': self.smiles_remove_atom_mapping,
            'canonicalize': self.smiles_canonicalize,
            'standardize': self.smiles_standardize,
            'augment': self.smiles_augment,
        }

    def convert(self, content, input_format, output_format):
        smiles = self.conversion_to_smiles[input_format](content)
        if smiles:
            output = self.conversion_from_smiles[output_format](smiles)
        else:
            output = False
        print(output, smiles)
        return output, smiles

    # to SMILES conversions
    def cas_to_smiles(self, cas):
        try:
            smiles = cas_to_smiles_api(cas)
            return smiles
        except:
            return False

    def iupac_to_smiles(self, iupac):
        try:
            smiles = iupac_to_smiles_api(iupac)
            return smiles
        except:
            return False

    def image_to_smiles(self, image, image_or_path='image'):
        if image_or_path == 'image':
            timestamp = datetime.now().strftime('%Y-%m-%d %H-%M-%S')
            file_name = f'{timestamp}_clipboard_image.png'
            file_path = os.path.join(self.image_input_dir, file_name)
            image.save(file_path, "PNG")
        else:
            file_path = image
        with torch.no_grad():
            output = self.model.predict_image_file(
                str(file_path), return_atoms_bonds=False, return_confidence=True)
        smiles = output['smiles']
        confidence = output['confidence']

        os.remove(file_path)
        if confidence > self.image_confidence_level:
            return smiles
        else:
            # rumps.notification(
            #     f'🔁 Confidence lower than: {self.image_confidence_level} ({confidence:.2f})',
            #     'Retake screenshot', f"Proposed SMILES: {smiles}", sound=False, icon='pictograms/carlos_helper_bad.png')
            return False

    def smiles_to_smiles(self, smiles):
        return smiles

    def inchi_to_smiles(self, inchi):
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                return smiles
        except:
            return False

    def mol_to_smiles(self, mol):
        try:
            smiles = Chem.MolToSmiles(Chem.MolFromMolBlock(mol))
            print(smiles)
            return smiles

        except:
            return False

    def selfies_to_smiles(self, selfies):
        try:
            smiles = sf.decoder(selfies)
            return smiles
        except:
            return False

    # from SMILES conversions
    def smiles_to_cas(self, smiles):
        try:
            cas = smiles_to_cas_api(smiles)
            return cas
        except:
            return False

    def smiles_to_iupac(self, smiles):
        try:
            iupac = smiles_to_iupac_api(smiles)
            return iupac
        except:
            return False

    def smiles_to_image(self, smiles):
        file_name, timestamp = create_filename_from_smiles(smiles)
        drawer = chemdraw.Drawer(smiles)
        image_path = os.path.join(self.image_generated_dir,file_name)
        drawer.draw_img(image_path)
        return image_path

    def smiles_to_price(self, smiles):
        cas = self.smiles_to_cas(smiles)
        if cas:
            print('cas for price conversion found')
            price_dict = self.vendors.cas_to_price(cas)
            print('price_dict',price_dict)

            return price_dict

    def smiles_to_selfies(self, smiles):
        try:
            selfies = sf.encoder(smiles)
            return selfies
        except:
            False

    def smiles_to_mol(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol_string = Chem.MolToMolBlock(mol)
            return mol_string
        except:
            return False

    def smiles_to_inchi(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            inchi = Chem.MolToInchi(mol)
            return inchi
        except:
            return False

    # SMILES to SMILES conversions

    def smiles_remove_atom_mapping(self, smiles):
        try:
            return remove_atom_mapping(smiles)
        except:
            return False

    def smiles_canonicalize(self, smiles):
        # try:
        return canonicalize_smiles(smiles)
        # except:
        #     return False

    def smiles_standardize(self, smiles):
        try:
            return smiles
        except:
            return False

    def smiles_augment(self, smiles):
        try:
            return augemnt_smiles(smiles)
        except:
            return False
        