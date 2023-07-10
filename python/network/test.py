import json
from python.classes.Species import Species, species_dict, species_list
from json import JSONEncoder


class Employee:
    def __init__(self, name, salary, address):
        self.name = name
        self.salary = salary
        self.address = address


class Address:
    def __init__(self, city, street, pin):
        self.city = city
        self.street = street
        self.pin = pin

# subclass JSONEncoder


class EmployeeEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class Protein:

    def __init__(self, gene_id: str):
        self.gene_id = gene_id
        self.t_ids = set()
        self.p_ids = set()
        self.ncbis = set()
        self.swissprots = set()
        self.trembls = set()
        self.names = set()
        self.refseqs = set()


class ProteinEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, Protein):
            return ProteinEncoder.default(self, obj)
        return json.JSONEncoder.default(self, obj)


address = Address("Alpharetta", "7258 Spring Street", "30004")
employee = Employee("John", 9000, address)

print("Printing to check how it will look like")
print(EmployeeEncoder().encode(employee))

print("Encode Employee Object into JSON formatted Data using custom JSONEncoder")
employeeJSONData = json.dumps(employee, indent=4, cls=EmployeeEncoder)
print(employeeJSONData)

# Let's load it using the load method to check if we can decode it or not.
print("Decode JSON formatted Data")
employeeJSON = json.loads(employeeJSONData)
print(employeeJSON)

# Protein('Q0085',{'Q0085_mRNA'},{'Q0085'},{'854601'},{'P00854'},set(),{'ATP6'},{'NP_009313.1'})
p = Protein('Q0085')
print('test')
print(SetEncoder().encode(p))

# json_string = json.dumps(p.__dict__, cls=SetEncoder)


for species1 in species_list:
    for species2 in species_list:
        species_short1 = species1.name
        species_short2 = species2.name

        if species1.name > species2.name:
            species_short1, species_short2 = species_short2, species_short1

        print(
            f'./sana -fg1 ../PPI-Network-Alignment/networks/SANA/{species1.short_name.lower().replace(" ", "_")}.network-109-4.4.222.el -fg2 ../PPI-Network-Alignment/networks/SANA/{species2.short_name.lower().replace(" ", "_")}.network-109-4.4.222.el -ec 1.0 -o {species_short1}_{species_short2}')
