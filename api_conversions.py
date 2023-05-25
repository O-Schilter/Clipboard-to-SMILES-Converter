import requests
import pubchempy as pcp


import json


import re




def structure_converter(structure,output_format,input_format):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{structure}/{output_format}"
    try:
        response = requests.get(url,verify=False)
        response.raise_for_status()
        result = response.text.split('\n')
        return result[0]
    except requests.exceptions.RequestException as e:
        print(f"An error occurred with cactus: {e}")
        print('conversion with PubChemy')
        if input_format == 'smiles':
            results = pcp.get_compounds(structure, 'smiles')
        else: 
            results = pcp.get_compounds(structure, 'name')

        if results:
            if results[0].cid:
                if output_format=='smiles':
                    return results[0].canonical_smiles
                elif output_format=='iupac':
                    return results[0].iupac_name
                elif output_format == 'cas':
                    synonyms = results[0].synonyms
                    cas = next((synonym.split(' ')[-1] for synonym in synonyms if re.match(r'^\d+-\d+-\d+$', synonym)), None)
                    if cas:
                        return cas
                return False
            else:
                return False
        else:
            return False


def smiles_to_cas_api(smiles):
    return structure_converter(smiles,'cas','smiles')

def smiles_to_iupac_api(smiles):
    return structure_converter(smiles,'iupac','smiles')

def cas_to_iupac_api(cas):
    return structure_converter(cas,'iupac','cas')

def iupac_to_cas_api(iupac):
    return structure_converter(iupac,'cas','iupac')

def iupac_to_smiles_api(iupac):
    return structure_converter(iupac,'smiles','iupac')

def cas_to_smiles_api(cas):
    return structure_converter(cas,'smiles','cas')



def cas_to_price_api(cas_number, vendor):
    print('finding price for:',cas_number)
    if vendor == 'ChemieBrunschwieg':
        url = 'https://www.chemie-brunschwig.ch/actions/sapShop/products/get-products'
        headers = {
            'Content-Type': 'application/json',
            'Accept': 'application/json, text/plain, */*',
            'X-Requested-With': 'XMLHttpRequest',
        }
        data = {
            "fulltext": str(cas_number),
            'brandid': '',
            'skip': 0,
            'take': 500,
            'startWithSearch': True,
            'useOciSession': False
        }
        
        response = requests.post(url, headers=headers, data=json.dumps(data),verify=False)
        response.raise_for_status()
        products = response.json()['Items']

        cheapest_product = None
        fastest_product = None
        lowest_price_per_amount = float('inf')
        lowest_lead_time = float('inf')
        try:
            for product in products:
                product_cas_number = product['CAS']
                pack_size = product['Packsize']
                price = product['Preis']
                lead_time = product['LeadTime']

                if product_cas_number == cas_number and pack_size and price:
                    amount = extract_amount(pack_size)
                    price_per_amount = price / amount

                    if price_per_amount < lowest_price_per_amount:
                        lowest_price_per_amount = price_per_amount
                        cheapest_product = product
                        cheapest_product['price_per'] = price_per_amount
                        cheapest_product['link'] = 'https://www.chemie-brunschwig.ch/shop/?term='+cheapest_product['ItemCode']
                    if lead_time < lowest_lead_time or (lead_time == lowest_lead_time and price_per_amount < lowest_price_per_amount):
                        lowest_lead_time = lead_time
                        fastest_product = product
            return {k: v for k, v in cheapest_product.items() if k in ['ItemName','Preis','link','price_per','Packsize']}

        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
            return None

import re

def extract_amount(pack_size):
    units = {
        'mg': 0.001,
        'g': 1,
        'gr':1,
        'kg': 1000,
        'ml': 1,
        'l': 1000
    }

    regex = r'(\d+(?:\.\d+)?)\s*(?:x\s*)?(\d+(?:\.\d+)?)?\s*([a-zA-Z]+)\s*'
    matches = re.findall(regex, pack_size)
    grams = 0

    if matches:
        for match in matches:
            quantity = float(match[0])
            factor = float(match[1]) if match[1] else 1
            unit = match[2].lower()

            if unit in units:
                conversion_factor = units[unit]
                grams += quantity * factor * conversion_factor

    return grams




if __name__ == '__main__':

    # SMILES to IUPAC name translation

    # SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    # IUPAC_name = translate_forward(SMILES)
    # print("IUPAC name of "+SMILES+" is: "+IUPAC_name)

    # Create an instance of the ClipboardImageSaverApp class
    cas_list = smiles_to_cas_api('CC(=O)Oc1ccccc1C(=O)O')
    # print(cas_list)
    print(cas_to_price_api(cas_list,vendor='ChemieBrunschwieg'))





    # print(smiles_to_iupac('O=C(C)Oc1ccccc1C(=O)O'))
