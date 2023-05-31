import re
import requests
import pubchempy as pcp


def structure_converter(structure, output_format, input_format):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{structure}/{output_format}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        result = response.text.split('\n')
        if len(result) == 1:
            return result[0]
    except requests.exceptions.RequestException as e:
        print(f"An error occurred with cactus: {e}")
    print('Try conversion with PubChemy')
    if input_format == 'smiles':
        results = pcp.get_compounds(structure, 'smiles')
    else:
        results = pcp.get_compounds(structure, 'name')

    if results:
        if results[0].cid:
            if output_format == 'smiles':
                return results[0].canonical_smiles
            elif output_format == 'iupac':
                return results[0].iupac_name
            elif output_format == 'cas':
                synonyms = results[0].synonyms
                cas = next((synonym.split(
                    ' ')[-1] for synonym in synonyms if re.match(r'^\d+-\d+-\d+$', synonym)), None)
                if cas:
                    return cas
            return False
        else:
            return False
    else:
        return False


def smiles_to_cas_api(smiles):
    return structure_converter(smiles, 'cas', 'smiles')


def smiles_to_iupac_api(smiles):
    return structure_converter(smiles, 'iupac', 'smiles')


def cas_to_iupac_api(cas):
    return structure_converter(cas, 'iupac', 'cas')


def iupac_to_cas_api(iupac):
    return structure_converter(iupac, 'cas', 'iupac')


def iupac_to_smiles_api(iupac):
    return structure_converter(iupac, 'smiles', 'iupac')


def cas_to_smiles_api(cas):
    return structure_converter(cas, 'smiles', 'cas')
