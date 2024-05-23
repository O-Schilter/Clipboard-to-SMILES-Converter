from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors

def create_filename_from_smiles(smiles):
    timestamp = datetime.now().strftime('%Y-%m-%d %H-%M-%S')
    smiles_for_filename = "".join(x for x in smiles if x.isalnum())
    new_filename = f'{timestamp}_{smiles_for_filename}.png'
    return new_filename, timestamp

def smiles_to_molecular_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol_weight = Descriptors.ExactMolWt(mol)
        return {
            "LogP": Chem.Descriptors.MolLogP(mol),
            "TPSA": Chem.Descriptors.TPSA(mol),
            "H-bond Donor": Chem.Descriptors.NumHDonors(mol),
            "H-bond Acceptor": Chem.Descriptors.NumHAcceptors(mol),
            "Ring Count": Chem.Descriptors.RingCount(mol),
            "Heavy Atom Count": Chem.Descriptors.HeavyAtomCount(mol),
            '[M]': mol_weight,
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
