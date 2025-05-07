#----------SOME UNIVERSAL DEFINITIONS----------#
OxygenNum = {'Al2O3': 3,'CaO': 1, 'Cl': 0, 'CO2': 2, 'CoO': 1, 'Cr2O3': 3,'CuO': 1, 'F': 0,
             'Fe2O3': 3,'FeO': 1, 'Ga2O3': 3,'H2O': 1, 'K2O': 1, 'MgO': 1, 'MnO': 1, 'MoO2': 2,
             'Na2O': 1, 'Nb2O5': 5,'NiO': 1, 'O': 1, 'P2O5': 5, 'PbO': 1, 'ReO3': 3, 'S': 0,
             'SiO2': 2, 'Ta2O5': 5,'TeO2': 2, 'ThO': 1, 'TiO2': 2, 'UO2': 2, 'V2O3': 3, 'WO3': 3,
             'ZnO': 1}

CationNum = {'Al2O3': 2,'CaO': 1, 'Cl': 1, 'CO2': 1, 'CoO': 1, 'Cr2O3': 2,'CuO': 1, 'F': 1,
             'Fe2O3': 2,'FeO': 1, 'Ga2O3': 2,'H2O': 2, 'K2O': 2, 'MgO': 1, 'MnO': 1, 'MoO2': 1,
             'Na2O': 2, 'Nb2O5': 2,'NiO': 1, 'O': 0, 'P2O5': 2, 'PbO': 1, 'ReO3': 1, 'S': 1,
             'SiO2': 1, 'Ta2O5': 2,'TeO2': 1, 'ThO': 1, 'TiO2': 1, 'UO2': 1, 'V2O3': 2, 'WO3': 1,
             'ZnO': 1}

"""Names"""
oxides = ['Al2O3','CaO', 'Cl', 	'CO2', 'CoO', 'Cr2O3','CuO', 'F', 'Fe2O3','FeO', 'Ga2O3','H2O', 
		  'K2O', 'MgO', 'MnO', 'MoO2','Na2O', 'Nb2O5','NiO', 'O', 'P2O5', 'PbO', 'ReO3', 'S',
          'SiO2', 'Ta2O5','TeO2', 'ThO', 'TiO2', 'UO2', 'V2O3', 'WO3', 'ZnO']

elements = ['Al','Ca', 'Cl', 'C', 'Co', 'Cr','Cu', 'F', 'Fe3','Fe', 'Ga','H', 'K', 'Mg', 'Mn', 'Mo',
            'Na', 'Nb','Ni', 'O', 'P', 'Pb', 'Re', 'S', 'Si', 'Ta','Te', 'Th', 'Ti', 'U', 'V', 'W',
            'Zn']

major_oxides = ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O', 'MnO', 'P2O5']

volatiles = ['CO2', 'H2O', 'F', 'Cl', 'S']

"""Transformations"""
oxides_to_cations = dict(zip(oxides, elements))
cations_to_oxides = dict(zip(elements, oxides))

"""Masses"""
trunc_elements = elements.copy()
trunc_elements.remove('Fe3')
elementMass = {i : element(i).atomic_weight for i in trunc_elements}
elementMass.update({'Fe3': element('Fe').atomic_weight})
oxideMass = {k: CationNum[k]*elementMass[v]+OxygenNum[k]*element('O').atomic_weight for k, v in oxides_to_cations.items()}

"""Other standard values"""
corgne_species = ['Ni', 'Cu', 'Si', 'Mn', 'Cr', 'Ga', 'Nb', 'Ta']
standard_interactions = ['S', 'C', 'O', 'Ni', 'Cu', 'Si', 'Mn', 'Cr', 'Ga', 'Nb', 'Ta']
norris_interactions = ['Fe', 'C', 'O', 'Si', 'P', 'S', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Ga',
                       'Ge', 'Zr', 'Nb', 'Mo', 'Hf', 'Ta', 'W']

#---------ERROR HANDLING--------#
class Error(Exception):
	"""Base class for exceptions in this module."""
	pass


class InputError(Error):
	"""Exception raised for errors in the input.

	Attributes:
		expression -- input expression in which the error occurred
		message -- explanation of the error
	"""

	def __init__(self, message):
		self.message = message

class GeneralError(Error):
	"""Exception raised for errors in the input.

	Attributes:
		expression -- input expression in which the error occurred
		message -- explanation of the error
	"""

	def __init__(self, message):
		self.message = message