
import requests
import json
import re


class Vendor():
    def __init__(self, name, url):
        self.name = name
        self.url = url


class ChemieBrunschwieg(Vendor):
    def __init__(self, name='ChemieBrunschwieg', url='https://www.chemie-brunschwig.ch/actions/sapShop/products/get-products'):
        super().__init__(name, url)

    def cas_to_price(self, cas_number):
        print('finding price for:', cas_number)
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

        response = requests.post(
            self.url, headers=headers, data=json.dumps(data), verify=False)
        response.raise_for_status()
        products = response.json()['Items']

        cheapest_product = None
        lowest_price_per_amount = float('inf')
        try:
            for product in products:
                product_cas_number = product['CAS']
                pack_size = product['Packsize']
                price = product['Preis']

                if product_cas_number == cas_number and pack_size and price:
                    amount = self.extract_amount(pack_size)
                    price_per_amount = price / amount

                    if price_per_amount < lowest_price_per_amount:
                        lowest_price_per_amount = price_per_amount
                    # if lead_time < lowest_lead_time or (lead_time == lowest_lead_time and price_per_amount < lowest_price_per_amount):
                    #     lowest_lead_time = lead_time
                    #     fastest_product = product

                        cheapest_product = {
                            "item_name": product['ItemName'],
                            "price": product['Preis'],
                            "price_per": price_per_amount,
                            "amount": product['Packsize'],
                            "link": 'https://www.chemie-brunschwig.ch/shop/?term='+product['ItemCode']
                        }
            return cheapest_product

        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
            return None

    def extract_amount(self, pack_size):
        units = {
            'mg': 0.001,
            'g': 1,
            'gr': 1,
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


class Enamine(Vendor):
    def __init__(self, name='Enamine', url='https://new.enaminestore.com/api/v1/'):
        super().__init__(name, url)

    def cas_to_price(self, cas_number):
        # Construct the first API endpoint to get catalog information
        endpoint_1 = f"{self.url}catalog/?q={cas_number}&pageSize=500"
        response_1 = requests.get(endpoint_1)

        if response_1.status_code == 200:
            data_1 = response_1.json()
            search_results = data_1["searchResults"]

            # Check if there are search results for the CAS number
            if search_results:
                # Get the code from the first search result
                code = search_results[0]["code"]

                # Construct the second API endpoint to get price information
                endpoint_2 = f"{self.url}catalog/price?id={code}&cat=BB&cur=USD"
                response_2 = requests.get(endpoint_2)

                if response_2.status_code == 200:
                    data_2 = response_2.json()
                    samples = data_2["samples"]

                    # Find the cheapest price per amount based on grams (g) or kilograms (kg)
                    cheapest_price_per_gram = float("inf")

                    for sample in samples:
                        amount = sample["amount"]
                        measure = sample["measure"]
                        price = sample["price"]

                        if measure == "mg":
                            amount /= 1000  # Convert milligrams to grams
                            unit = "g"
                        elif measure == "kg":
                            amount *= 1000  # Convert kilograms to grams
                            unit = "g"
                        else:
                            unit = measure

                        # Calculate the price per gram
                        price_per_gram = price / amount

                        if price_per_gram < cheapest_price_per_gram:
                            cheapest_price_per_gram = price_per_gram

                            cheapest_product = {
                                "item_name": data_1["searchResults"][0]["name"],
                                "price": price,
                                "price_per": price_per_gram,
                                "amount": f"{amount} {unit}",
                                "link": f"https://new.enaminestore.com/catalog/{code}"
                            }

                    return cheapest_product

                else:
                    # Handle API error for the second endpoint
                    print("Error retrieving price information")
            else:
                # No search results for the CAS number
                print("No catalog information found for the CAS number")
        else:
            # Handle API error for the first endpoint
            print("Error retrieving catalog information")

        return None


class Vendors():
    def __init__(self, vendors):
        self.vendors = vendors
        self.cas_to_price = self.vendors[0].cas_to_price

    def select_vendor(self, vendor):
        for ven in self.vendors:
            if ven.name == vendor:
                self.cas_to_price = ven.cas_to_price


if __name__ == '__main__':

    vend = Vendors([Enamine(), ChemieBrunschwieg()])
    vend.select_vendor('ChemieBrunschwieg')
    print(vend.cas_to_price('64-17-5'))

