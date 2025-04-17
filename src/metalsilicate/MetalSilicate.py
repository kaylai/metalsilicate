import pandas as pd
import numpy as np
import interactionparameters as ip
import activitycoefficients as ac
import opticalbasicity as ob
# import statsmodels.api as sm
# import statsmodels.formula.api as smf
# from patsy import dmatrices
from mendeleev import element
from mendeleev.fetch import fetch_table
import warnings as w
import sys
from time import sleep

#---------SAND BOX----------#
"""
SAND BOX
========

data = pandas.read_excel('data.xlsx', sheet_name="ForPython", header=[0,1]) #import excel file 'data.xlsx'
data = data.dropna(axis=1, how='all') #drop all columns that contain no data
data = data.loc[:, (data != 0).any(axis=0)] #drop all columns that contain all zeroes

formula = '{} ~ {} -1'.format(data.columns[0], ' + '.join(data.columns[1:])) #defines the formula needed to put into the dmatrices call without knowing how many X variables there will be but assuming the y variable is always the first column


# print(formula)
# y, X = dmatrices(formula, data=data, return_type='dataframe') #needed to align matrices

# model = sm.OLS(y, X)
# res = model.fit()

# print(res.summary())


ORDER OF OPERATIONS
-------------------
Correct for interactions of solutes in Fe base alloys (after Wade and Wood, 2005). 
		NOTE: do we want to calculate activities for a bunch of elements in the metal? Then we could use this instead of
		the straight up mole fractions or wt% conentrations in the regression. This might be best if we want to compare
		data with a large range of metal compositions (e.g. in terms of S, C, Si). So we can calculate the activities of
		at least those components, and use X for the others?

	1. Lookup EPSILON_i_j at a reference T
	2. Calculate EPSIOLON_i_j at sample T
	3. Lookup GAMMA_i_0 at reference T
	4. Calculate GAMMA_i_0 at sample T
	5. Calculate activity coefficient GAMMA_i from ref value GAMMA_i_0 and interaction parameter EPSILON_i_j at sample T
	6. Calculate the activity of the solute as mole fraction times gamma.

Calculate the fO2:
	1. Calcualte Fe activity in the metal
	2. Calculate FeO activity in the silicate
	3. Calculate deltaIW fO2 from activity from mole fraction and GAMMA of FeO in silicate and Fe in metal (see Corgne et al., 2007)
"""

#----------SOME UNIVERSAL DEFINITIONS----------#
OxygenNum = {'Al2O3': 3,'CaO': 1, 'Cl': 0, 'CO2': 2, 'CoO': 1, 'Cr2O3': 3,'CuO': 1, 'F': 0, 'Fe2O3': 3,'FeO': 1, 'Ga2O3': 3,'H2O': 1, 
		  'K2O': 1, 'MgO': 1, 'MnO': 1, 'MoO2': 2,'Na2O': 1, 'Nb2O5': 5,'NiO': 1, 'O': 1, 'P2O5': 5, 'PbO': 1, 'ReO3': 3, 'S': 0, 'SiO2': 2, 
		  'Ta2O5': 5,'TeO2': 2, 'ThO': 1, 'TiO2': 2, 'UO2': 2, 'V2O3': 3, 'WO3': 3, 'ZnO': 1}

CationNum = {'Al2O3': 2,'CaO': 1, 'Cl': 1, 'CO2': 1, 'CoO': 1, 'Cr2O3': 2,'CuO': 1, 'F': 1, 'Fe2O3': 2,'FeO': 1, 'Ga2O3': 2,'H2O': 2, 
		  'K2O': 2, 'MgO': 1, 'MnO': 1, 'MoO2': 1,'Na2O': 2, 'Nb2O5': 2,'NiO': 1, 'O': 0, 'P2O5': 2, 'PbO': 1, 'ReO3': 1, 'S': 1, 'SiO2': 1, 
		  'Ta2O5': 2,'TeO2': 1, 'ThO': 1, 'TiO2': 1, 'UO2': 1, 'V2O3': 2, 'WO3': 1, 'ZnO': 1}

"""Names"""
oxides = ['Al2O3','CaO', 'Cl', 	'CO2', 'CoO', 'Cr2O3','CuO', 'F', 'Fe2O3','FeO', 'Ga2O3','H2O', 
		  'K2O', 'MgO', 'MnO', 'MoO2','Na2O', 'Nb2O5','NiO', 'O', 'P2O5', 'PbO', 'ReO3', 'S', 'SiO2', 
		  'Ta2O5','TeO2', 'ThO', 'TiO2', 'UO2', 'V2O3', 'WO3', 'ZnO']
elements = ['Al','Ca', 'Cl', 'C', 'Co', 'Cr','Cu', 'F', 'Fe3','Fe', 'Ga','H', 
		  'K', 'Mg', 'Mn', 'Mo','Na', 'Nb','Ni', 'O', 'P', 'Pb', 'Re', 'S', 'Si', 
		  'Ta','Te', 'Th', 'Ti', 'U', 'V', 'W', 'Zn']
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
norris_interactions = ['Fe', 'C', 'O', 'Si', 'P', 'S', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Ga', 'Ge', 'Zr', 'Nb', 'Mo', 'Hf', 'Ta', 'W']

def status_bar(percent, sample_name, barLen=20):
	"""
	Prints a status bar to the terminal.

	percent: float
		Percent value of progress from 0 to 1

	barLen: int
		Length of bar to print
	""" 
	sys.stdout.write("\r")
	sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent), barLen, percent * 100))
	sys.stdout.write("  Working on sample " + str(sample_name))
	if percent == 1.0:
		sys.stdout.write("\n")
	sys.stdout.flush()

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

#---------SOME DATA TRANSFORMATION METHODS--------#
def get_oxides(sample):
	"""
	Returns a sample composition with only compositional oxide data, removing any extranneous data.
	Useful when passing a self-defined sample (e.g. dict or pandas Series) to a some other function.

	Parameters
	----------
	sample: pandas Series, dictionary
		A sample composition plus other sample information

	Returns
	-------
	Sample passed as > Returned as
			pandas Series > pandas Series
			dictionary > dictionary

				Sample composition with extranneous information removed.
	"""
	_sample = sample.copy()
	for oxide in oxides:
		if oxide in _sample.keys():
			pass
		else:
			_sample[oxide] = 0.0

	clean = {oxide:  _sample[oxide] for oxide in oxides}

	if isinstance(_sample, dict):
		return clean
	if isinstance(_sample, pd.core.series.Series):
		return pd.Series(clean)

def get_elements(sample):
	"""
	Returns a sample composition with only compositional element data, removing any extranneous data.
	Useful when passing a self-defined sample (e.g. dict or pandas Series) to a some other function.

	Parameters
	----------
	sample: pandas Series, dictionary
		A sample composition plus other sample information

	Returns
	-------
	Sample passed as > Returned as
			pandas Series > pandas Series
			dictionary > dictionary

				Sample composition with extranneous information removed.
	"""
	_sample = sample.copy()
	for element in elements:
		if element in _sample.keys():
			pass
		else:
			_sample[element] = 0.0

	clean = {element:  _sample[element] for element in elements}

	if isinstance(_sample, dict):
		return clean
	if isinstance(_sample, pd.core.series.Series):
		return pd.Series(clean)

def normalize(sample):
	"""Normalizes an input composition to 100%. This is the 'standard' normalization routine.

	Parameters
	----------
	sample:    pandas Series, dictionary, pandas DataFrame, or ExcelFile object
		A single composition can be passed as a dictionary. Multiple compositions can be passed either as
		a pandas DataFrame or an ExcelFile object. Compositional information as oxides must be present.

	Returns
	-------
	Sample passed as > Returned as
		pandas Series > pandas Series
		dictionary > dictionary
		pandas DataFrame > pandas DataFrame
		ExcelFile object > pandas DataFrame

			Normalized major element oxides.
	"""
	def single_normalize(sample):
		single_sample = sample
		return {k: 100.0 * v / sum(single_sample.values()) for k, v in single_sample.items()}

	def multi_normalize(sample):
		multi_sample = sample.copy()
		multi_sample["Sum"] = sum([multi_sample[oxide] for oxide in oxides])
		for column in multi_sample:
			if column in oxides:
				multi_sample[column] = 100.0*multi_sample[column]/multi_sample["Sum"]

		del multi_sample["Sum"]
		return multi_sample

	if isinstance(sample, dict):
		_sample = sample.copy()
		return single_normalize(_sample)
	elif isinstance(sample, pd.core.series.Series):
		_sample = pd.Series(sample.copy())
		sample_dict = sample.to_dict()
		return pd.Series(single_normalize(sample_dict))
	elif isinstance(sample, ExcelFile):
		_sample = sample
		data = _sample.data
		return multi_normalize(data)
	elif isinstance(sample, pd.DataFrame):
		return multi_normalize(sample)

def wtpercentElements_to_molElements(sample_elements):
	""" Takes in a pandas Series or dict containing major elements in wt%, and converts it
	to molar proportions (normalised to 1).

	Parameters
	----------
	sample_elements         dict or pandas Series
		Major elements in wt% (as elements, not oxides)

	Returns
	-------
	dict or pandas Series
		Molar proportions of major elements, normalised to 1.
	"""
	molElements = {}
	if type(sample_elements) == dict or type(sample_elements) == pd.core.series.Series:
		if type(sample_elements) == dict:
			elementslist = list(sample_elements.keys())
		elif type(sample_elements) == pd.core.series.Series:
			elementslist = list(sample_elements.index)

		for element in elementslist:
			molElements[element] = sample_elements[element]/elementMass[element]

		if type(sample_elements) == pd.core.series.Series:
			molElements = pd.Series(molElements)
			molElements = molElements/molElements.sum()
		else:
			total = np.sum(list(molElements.values()))
			for element in elementslist:
				molElements[element] = molElements[element]/total

		return molElements

	elif isinstance(sample, pd.DataFrame):
		data = sample
		for key, value in elementMass.items():
			data.loc[:, key] /= value

		data["MPSum"] = sum([data[element] for element in sample_elements])

		for element in sample_elements:
			data.loc[:, element] /= data['MPSum']
		del data['MPSum']

		return data

	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")

def wtpercentOxides_to_molOxides(sample_oxides):
	""" Takes in a pandas Series or dict containing major element oxides in wt%, and converts it
	to molar proportions (normalised to 1).

	Parameters
	----------
	sample_oxides         dict or pandas Series
		Major element oxides in wt%

	Returns
	-------
	dict or pandas Series
		Molar proportions of major element oxides, normalised to 1.
	"""
	molOxides = {}
	if type(sample_oxides) == dict or type(sample_oxides) == pd.core.series.Series:
		if type(sample_oxides) == dict:
			oxideslist = list(sample_oxides.keys())
		elif type(sample_oxides) == pd.core.series.Series:
			oxideslist = list(sample_oxides.index)

		for ox in oxideslist:
			molOxides[ox] = sample_oxides[ox]/oxideMass[ox]

		if type(sample_oxides) == pd.core.series.Series:
			molOxides = pd.Series(molOxides)
			molOxides = molOxides/molOxides.sum()
		else:
			total = np.sum(list(molOxides.values()))
			for ox in oxideslist:
				molOxides[ox] = molOxides[ox]/total

		return molOxides

	elif isinstance(sample, pd.DataFrame):
		data = sample
		for key, value in oxideMass.items():
			data.loc[:, key] /= value

		data["MPOSum"] = sum([data[oxide] for oxide in sample_oxides])

		for oxide in sample_oxides:
			data.loc[:, oxide] /= data['MPOSum']
		del data['MPOSum']

		return data

	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")

def wtpercentOxides_to_molSingleO(sample_oxides,exclude_volatiles=False):
	""" Takes in a pandas Series containing major element oxides in wt%, and constructs
	the chemical formula, on a single oxygen basis.

	Parameters
	----------
	sample_oxides         dict or pandas Series
		Major element oxides in wt%

	Returns
	-------
	dict or pandas Series
		The chemical formula of the composition, on a single oxygen basis. Each element is
		a separate entry in the Series.
	"""
	molCations = {}
	if type(sample_oxides) == dict:
		oxideslist = list(sample_oxides.keys())
	elif type(sample_oxides) == pd.core.series.Series:
		oxideslist = list(sample_oxides.index)
	else:
		raise InputError("The composition input must be a pandas Series or dictionary.")

	total_O = 0.0
	for ox in oxideslist:
		if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
			cation = oxides_to_cations[ox]
			molCations[cation] = CationNum[ox]*sample_oxides[ox]/oxideMass[ox]
			total_O += OxygenNum[ox]*sample_oxides[ox]/oxideMass[ox]
	if type(sample_oxides) == pd.core.series.Series:
		molCations = pd.Series(molCations)
		molCations = molCations/total_O
	else:
		# total = np.sum(list(molCations.values()))
		for ox in oxideslist:
			if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
				cation = oxides_to_cations[ox]
				molCations[cation] = molCations[cation]/total_O

	return molCations

def oxides_to_elements(sample_oxides):
	"""
	Converts composition in oxides to elements

	Parameters
	----------
	sample_oxides: dict
		Major element oxides in wt%
	"""
	converted_sample = {}
	original_sum = sum(sample_oxides.values())

	for oxide in sample_oxides:
		element = oxides_to_cations[oxide]
		conversion_factor = CationNum[oxide]*elementMass[element]/oxideMass[oxide]
		converted_sample[element] = sample_oxides[oxide] * conversion_factor

	new_sum = sum(converted_sample.values())
	converted_sample["O"] = original_sum - new_sum

	return converted_sample

def elements_to_oxides(sample_elements):
	"""
	Converts composition in elements to oxides

	Parameters
	----------
	sample_elements: dict
		Major element composition in wt%
	"""

	converted_sample = {}

	for element in sample_elements:
		oxide = cations_to_oxides[element]
		conversion_factor = CationNum[oxide]*elementMass[element]/oxideMass[oxide]
		converted_sample[oxide] = sample_elements[element] / conversion_factor

	return converted_sample

def preprocess_dataset(data, how):
		"""
		Adds 0.0 values to any oxide or element data not passed.

		Parameters
		----------
		data: pandas DataFrame
			Composition of samples in wt%

		how: str
			Describes if data are in terms of oxides or elements. Can be one of 'oxides', 'elements'.

		Returns
		-------
		pandas DataFrame
		"""
		_data = data.copy()

		if how == "oxides":
			datatype = oxides
		if how == "elements":
			datatype = elements

		for item in datatype:
			if item in _data.columns:
				pass
			else:
				_data[item] = 0.0

		return _data

def clean(data):
	"""
	Takes a pandas dataframe (e.g. myfile.data, myfile.silicate_data) and removes any columns with all 0's, any non-numeric data.

	Parameters
	----------
	data: pandas DataFrame
		A pandas DataFrame object.

	Returns
	-------
	pandas DataFrame
	"""
	_data = data.copy()

	_data = _data.apply(pd.to_numeric, errors='coerce')
	_data = _data.fillna(0) #fill in any missing data with 0's
	_data = _data.dropna(axis=1, how='all') #drop all columns that contain no data
	_data = _data.loc[:, (data != 0).any(axis=0)] #drop all columns that contain all zeroes
	_data = _data.loc[(_data!=0).any(axis=1)] #drop all rows that contain all zeroes

	return _data

def rename_duplicates(df, suffix='-duplicate-'):
	appendents = (suffix + df.groupby(level=0).cumcount().astype(str).replace('0','')).replace(suffix, '')
	return df.set_index(df.index.astype(str) + appendents)

#----------IMPORT FILE AND PROCESS---------#

class ExcelFile(object):
	"""An excel file with sample names and variables. File must have two header rows. The first header row breaks data into the following: data about
	the experiments themselves (e.g., pressure, temperature, duration, experimental apparatus); data describing the silicate melt composition in wt% oxides;
	data describing the first metal composition in wt% elements; data describing a second metal composition in wt% elements; any notes about the sample. The
	first column will become the index of all of these datasets, so by convention this should be the sample name.

	Attributes
	----------
		filename: str
			Path to the excel file, e.g., "my_file.xlsx". This always needs to be passed, even if the user is passing a pandas DataFrame
			rather than an Excel file.	

		sheet_name: str
			OPTIONAL. Default value is 0 which gets the first sheet in the excel spreadsheet file. This implements the pandas.
			read_excel() sheet_name parameter. But functionality to read in more than one sheet at a time (e.g., pandas.read_excel(sheet_name=None))
			is not yet imlpemented in VESIcal. From the pandas 1.0.4 documentation:
				Available cases:
					- Defaults to 0: 1st sheet as a DataFrame
					- 1: 2nd sheet as a DataFrame
					- "Sheet1": Load sheet with name “Sheet1”
		
		experimental_data: str
			Name of the header label corresponding to experimental metadata

		silicate_data: str
			Name of the header label corresponding to silicate melt composition in wt% oxides

		metal_data: str
			Name of the header label corresponding to the metal composition in wt% elements

		metal_2_data: str
			Name of the header label corresponding to a second metal composition in wt% elements

		notes: str
			Name of the header label corresponding to any notes about the sample
	"""

	def __init__(self, filename, sheet_name=0, 
					experimental_data="Experimental", silicate_data="Silicate_Melt_wt_oxides", 
					metal_data="Metal_1_wt_elements", metal_2_data="Metal_2_wt_elements", notes="Notes"):
		"""Return an ExcelFile object whoes parameters are defined here."""

		if isinstance(sheet_name, str) or isinstance(sheet_name, int):
			pass
		else:
			raise InputError("If sheet_name is passed, it must be of type str or int. Cannot import more than one sheet at a time.")

		data = pd.read_excel(filename, sheet_name=sheet_name, header=[0,1], index_col=0)
		data = rename_duplicates(data) #handle any duplicated sample names
		data = data.fillna(0) #fill in any missing data with 0's
		data = data.loc[(data!=0).any(axis=1)] #drop all rows that contain all zeroes

		experimental_data = data[experimental_data]
		silicate_data = data[silicate_data]
		silicate_data = clean(silicate_data)
		metal_data = data[metal_data]
		metal_data = clean(metal_data)
		try:
			metal_2_data = data[metal_2_data]
			metal_2_data = clean(metal_2_data)
		except:
			metal_2_data = pd.DataFrame(index=metal_data.index)
		notes = data[notes]

		self.data = data
		self.experimental_data = experimental_data
		self.silicate_data = silicate_data
		self.metal_data = metal_data
		self.metal_2_data = metal_2_data
		self.notes = notes

	def save_excelfile(self, filename, calculations, sheet_name=None):
		"""
		Saves data calculated by the user in batch processing mode (using the ExcelFile class methods) to an organized
		excel file, with the original user data plus any calculated data.

		Parameters
		----------
		filename: string
			Name of the file. Extension (.xlsx) should be passed along with the name itself, all in quotes (e.g., 'myfile.xlsx').

		calculations: pandas DataFrame or list of pandas DataFrames
			A single variable or list of variables containing calculated outputs from any of the core ExcelFile functions

		sheet_name: None, string, or list
			OPTIONAL. Default value is None. Allows user to set the name of the sheet or sheets written to the Excel file.

		Returns
		-------
		Excel File
			Creates and saves an Excel file with data from each calculation saved to its own sheet.
		"""
		if isinstance(calculations, list):
			if isinstance(sheet_name, list) or sheet_name is None:
				pass
		else:
			calculations = [calculations]
		with pd.ExcelWriter(filename) as writer:
			self.data.to_excel(writer, 'Original_User_Data')
			if sheet_name is None:
				for n, df in enumerate(calculations):
					df.to_excel(writer, 'Calc%s' % n)
			elif isinstance(sheet_name, list):
				pass
			else:
				sheet_name = [sheet_name]
			if isinstance(sheet_name, list):
				if len(sheet_name) == len(calculations):
					pass
				else:
					raise InputError("calculations and sheet_name must have the same length")

				for i in range(len(calculations)):
					if isinstance(sheet_name[i], str):
						calculations[i].to_excel(writer, sheet_name[i])
					else:
						raise InputError("if sheet_name is passed, it must be list of strings")

		return print("Saved " + str(filename))

	def get_silicate_comp(self, sample, norm='none'):
		"""
		Returns wt% oxide composition of a single silicate sample from a user-imported excel file as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied to the sample.

			'standard' normalizes the entire input composition (including any volatiles) to 100%.

			'fixedvolatiles' normalizes oxides to 100%, including volatiles. The volatile
			wt% will remain fixed, whilst the other major element oxides are reduced proportionally
			so that the total is 100 wt%.

			'additionalvolatiles' normalizes oxides to 100%, assuming it is volatile-free. If
			H2O or CO2 are passed to the function, their un-normalized values will be retained
			in addition to the normalized non-volatile oxides, summing to >100%.

			'none' returns the value-for-value un-normalized composition.

		Returns
		-------
		dictionary
			Composition of the silicate phase in wt% oxides
		"""
		if norm == 'none' or norm == 'standard':
			pass
		else:
			raise InputError('norm must be either none or standard.')

		data = preprocess_dataset(self.silicate_data, how='oxides')
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])

		sample_comp = get_oxides(sample_dict)

		if norm == 'standard':
			return normalize(sample_comp)
		if norm == 'none':
			return sample_comp

	def get_metal_comp(self, sample, metal_1=True, metal_2=False, norm='none'):
		"""
		Returns wt% element composition of a single metal sample from a user-imported excel file as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		metal_1: bool
			OPTIONAL. Detault is True. If True, returns composition of first metal "metal 1" (self.metal_data).
		
		metal_2: bool
			OPTIONAL. Detault is False. If True, returns composition of second metal "metal 2" (self.metal_2_data).
			Note if both metal_1 and metal_2 are set to True, this function returns to objects.

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied to the sample.

			'standard' normalizes the entire input composition (including any volatiles) to 100%.

			'fixedvolatiles' normalizes oxides to 100%, including volatiles. The volatile
			wt% will remain fixed, whilst the other major element oxides are reduced proportionally
			so that the total is 100 wt%.

			'additionalvolatiles' normalizes oxides to 100%, assuming it is volatile-free. If
			H2O or CO2 are passed to the function, their un-normalized values will be retained
			in addition to the normalized non-volatile oxides, summing to >100%.

			'none' returns the value-for-value un-normalized composition.

		Returns
		-------
		dictionary or dictionaries
			Composition of the metal phase in wt%. If both metal_1 and metal_2 are set to True, two objects
			are returned.
		"""
		if norm == 'none' or norm == 'standard':
			pass
		else:
			raise InputError('norm must be either none or standard.')

		if metal_1 == True:
			data = preprocess_dataset(self.metal_data, how='elements')
			my_sample = pd.DataFrame(data.loc[sample])
			sample_dict = (my_sample.to_dict()[sample])

			sample_comp = get_elements(sample_dict)

			if norm == 'standard':
				metal_1_return = normalize(sample_comp)
			if norm == 'none':
				metal_1_return = sample_comp

			if metal_2 == False:
				return metal_1_return
		
		if metal_2 == True:
			data2 = preprocess_dataset(self.metal_2_data, how='elements')
			my_sample = pd.DataFrame(data2.loc[sample])
			sample_dict = (my_sample.to_dict()[sample])

			sample_comp = get_elements(sample_dict)

			if norm == 'standard':
				metal_2_return = normalize(sample_comp)
			if norm == 'none':
				metal_2_return = sample_comp

			if metal_1 == False:
				return metal_2_return

		if metal_1 == True and metal_2 == True:
			return metal_1_return, metal_2_return

	def get_experimental_metadata(self, sample):
		"""
		Returns experimental metadata of a single metal sample from a user-imported excel file as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		Returns
		-------
		dictionary
			Metadata of the sample
		"""
		if norm == 'none' or norm == 'standard':
			pass
		else:
			raise InputError('norm must be either none or standard.')

		data = self.experimental_data
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])

		return sample_dict

	def process(self, pressure, temperature, species=corgne_species, interactions=standard_interactions, filename=None, print_status=True):
		"""
		Process the excel file and save a standardized set of data to an excel file that can be used for
		regression modelling via python or other programs.

		Parameters
		----------
		pressure: float or str
			Pressure, in GPa. Can be passed as float, in which case the passed value is used as the pressure for all samples. 
			Alternatively, pressure information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object.

		temperature: float or str
			Temperature, in degrees C. Can be passed as float, in which case the passed value is used as the temperature for all samples. 
			Alternatively, temperature information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object.

		species: str or list
			List of elements for which to calculate gamma values, if that list is different than the list of interactions.
			Default value is corgne_species, which is the set of species chosen that returns values most closely matching
			the reported values from the Corgne paper.

		interactions: list
			OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
			the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
			temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.

		filename: str or None
			OPTIONAL. If none, no Excel file will be created. To save the calculated data to an excel spreadsheet file, enter the 
			filename here as a string, including the extension (.xlsx).

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame + saved Excel file (optional)
			Dataframe containing original user input data plus:
				- gammaFe in metal
				- gamma solutes in metal
				- gammaFeO in silicate
				- NBO/T for silicate
				- fO2 as dIW
				- fO2 as log(fO2) absolute
				- X_silicate
				- X_metal
				- activity FeO silicate as XFeO*gammaFeO
				- activity Si in metal
			#TODO:
				- activity SiO2 silicate

		#TODO - could add option to output warnings from the gamma calcs. Could write them to be informative as to what elements
		were used for the gamma calculation (e.g., which interaction parameters).

		"""
		return_data = pd.DataFrame(index=self.data.index)

		#Convert silicate from wt% to mol fraction oxides
		if print_status == True:
			print("Converting silicate comp from wt% to X...")
		molOxides_Silicate = self.wtpercentOxides_to_molOxides_Silicate(print_status=print_status)

		#Convert metal 1 and metal 2 data from wt% to mole fraction elements
		if print_status == True:
			print("Converting metal comps from wt% to X...")
		molElements_metal_1, molElements_metal_2 = self.wtpercentElements_to_molElements_Metal(metal_1=True, metal_2=True, print_status=print_status)

		#Compute gammaFe in metal
		if print_status == True:
			print("Calculating gammaFe values...")
		Fe_gammas = self.calc_gamma_Fe_metal(temperature=temperature, interactions=interactions, print_status=print_status)

		#Compute gamma of solutes in metal
		if print_status == True:
			print("Calculating gamma values for solutes...")
		calcd_gammaFes = {index: row["gamma_Fe_metal"] for index, row in Fe_gammas.iterrows()}
		gamma_sol = self.calc_gamma_solute_metal(temperature=temperature, species=species, interactions=interactions, gammaFe=calcd_gammaFes, print_status=print_status)

		#Compute gammaFeO in silicate
		if print_status == True:
			print("Calculating gammaFeO in the silicate...")
		gammaFeO_calc = self.calc_gamma_FeO_silicate()

		#Compute aFeO in silicate as XFeO*gammaFeO
		if print_status == True:
			print("Calculating FeO activity in the melt as XFeO * gammaFeO...")
		aFeO_calc = pd.DataFrame(index=self.data.index)
		aFeO_calc["aFeO_Silicate"] = molOxides_Silicate["FeO"] * gammaFeO_calc["gamma_FeO_silicate"] 

		#Compute aSi in metal
		if print_status == True:
			print("Calculating Si activity in the metal...")
		aSi_calc = self.metal_activity_from_composition(species="Si", temperature=temperature, interactions=interactions)

		#Compute NBO/T for silicate
		if print_status == True:
			print("Calculting NBO/T values for silicates...")
		NBO_T = self.calc_NBO_T(metadata=True)

		#Compute fO2 as dIW
		if print_status == True:
			print("Calculating dIW values...")
		calcd_gammaFeOs = {index: row["gamma_FeO_silicate"] for index, row in gammaFeO_calc.iterrows()}
		dIW = self.calc_dIW(temperature=temperature, interactions=interactions, gammaFe=calcd_gammaFes, gammaFeO=calcd_gammaFeOs, print_status=print_status)

		#Compute fO2 as log(fO2) absolute
		if print_status == True:
			print("Calculating log(fO2) values...")
		calcd_dIWs = dIW["fO2_dIW"]
		log_fO2 = self.calc_logfO2(pressure=pressure, temperature=temperature, dIW=calcd_dIWs, print_status=print_status)

		return_data = pd.concat([return_data, Fe_gammas, gamma_sol, gammaFeO_calc, aFeO_calc, aSi_calc, NBO_T, dIW, 
								 log_fO2, molOxides_Silicate, molElements_metal_1, molElements_metal_2], axis=1)
		thermo_params = pd.concat([Fe_gammas, gamma_sol, gammaFeO_calc, aFeO_calc, aSi_calc, NBO_T, dIW, log_fO2], axis=1)

		if filename is not None:
			if isinstance(filename, str):
				self.save_excelfile(filename=filename, 
									calculations=[molOxides_Silicate, molElements_metal_1, molElements_metal_2, thermo_params], 
									sheet_name=["molOxides_Silicate", "molElements_Metal_1", "molElements_Metal_2", "Thermo_Params"])
			else:
				raise InputError("filename must be type str")

		return return_data

	def wtpercentOxides_to_molOxides_Silicate(self, print_status=False):
		"""
		Converts all silicate data in an ExcelFile object from wt% oxides to mol fraction oxides.

		Parameters
		----------
		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame
			DataFrame with converted silicate data.
		"""
		data = self.silicate_data
		user_data = data.copy()

		molOxides_vals = []
		iterno = 0
		for index, row in user_data.iterrows():
			iterno += 1
			if print_status == True:
				percent = iterno/len(data.index)
				status_bar(percent, index)
			wtpercentOxides = self.get_silicate_comp(index)
			calc = wtpercentOxides_to_molOxides(wtpercentOxides)
			molOxides_vals.append(calc)

		return_data = pd.DataFrame(molOxides_vals, index=user_data.index)

		return return_data

	def wtpercentElements_to_molElements_Metal(self, metal_1=True, metal_2=False, print_status=False):
		"""
		Converts all metal data in an ExcelFile object from wt% elements to mol fraction elements.

		Parameters
		----------
		metal_1: bool
			OPTIONAL. Detault is True. If True, returns composition of first metal "metal 1" (self.metal_data).
		
		metal_2: bool
			OPTIONAL. Detault is False. If True, returns composition of second metal "metal 2" (self.metal_2_data).
			Note if both metal_1 and metal_2 are set to True, this function returns to objects.

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		two pandas DataFrames
			DataFrames with converted metal data for metals 1 and 2.
		"""
		if metal_1 == True:
			data = self.metal_data
			user_data = data.copy()

			molElements_vals = []
			iterno = 0
			for index, row in user_data.iterrows():
				iterno += 1
				if print_status == True:
					percent = iterno/len(data.index)
					status_bar(percent, index)
				wtpercentElements = self.get_metal_comp(index, metal_1=True, metal_2=False)
				calc = wtpercentElements_to_molElements(wtpercentElements)
				molElements_vals.append(calc)

			metal_1_return = pd.DataFrame(molElements_vals, index=user_data.index)

		if metal_2 == True:
			data2 = self.metal_2_data
			user_data2 = data2.copy()

			molElements_vals2 = []
			iterno = 0
			for index, row in user_data2.iterrows():
				iterno += 1
				if print_status == True:
					percent = iterno/len(data2.index)
					status_bar(percent, index)
				wtpercentElements2 = self.get_metal_comp(index, metal_1=False, metal_2=True)
				calc2 = wtpercentElements_to_molElements(wtpercentElements2)
				molElements_vals2.append(calc2)

			metal_2_return = pd.DataFrame(molElements_vals2, index=user_data2.index)

		if metal_1 == True and metal_2 == False:
			return metal_1_return

		if metal_1 == False and metal_2 == True:
			return metal_2_return

		if metal_1 == True and metal_2 == True:
			return metal_1_return, metal_2_return

	def calc_NBO_T(self, method=4, metadata=False):
		"""
		Calcultes the NBO/T (non-bridging oxygens over tetrahedrally coordinated ions) for all samples in a dataset. Written after
		excel spreadsheet by FM McCubbin.

		Parameters
		----------
		method: int
			Default is 4. Integer from 1 to 4 inclusive. Use to choose which NBO/T parameterization to use as defined here:
				1: NF = Si, Al-mixed, Ti
				2: NF = Si, Al-mixed, Ti, P
				3: NF = Si, Al-mixed, Fe3+, Ti, P
				4: NF = Si, Al-mixed, Fe3+, Cr3+, Ti, P

		metadata: bool
			Default is False. If set to true, a dict is returned with more information calculated about the sample.

		Returns
		-------
		pandas DataFrame
		"""
		data = self.silicate_data
		user_data = data.copy()

		NBO_vals = []
		for index, row in user_data.iterrows():
			bulk_comp = self.get_silicate_comp(index)
			try:
				calc = calc_NBO_T(sample=bulk_comp, method=method, metadata=metadata)
			except:
				if metadata == True:
					calc = {"NBO/T": "Could not calculate NBO/T",
							"Mg#": "Could not calculate NBO/T",
							"Alkalinity Index": "Could not calculate NBO/T",
							"Metaluminous": "Could not calculate NBO/T",
							"Peraluminous": "Could not calculate NBO/T",
							"Peralkaline": "Could not calculate NBO/T"}
				if metadata == False:
					calc = "Could not calculate NBO/T"
			NBO_vals.append(calc)

		return_data = pd.DataFrame(index=user_data.index)

		if metadata == True:
			NBO_T_vals = []
			Mg_num_vals = []
			Alk_vals = []
			Meta_vals = []
			Pera_vals =[]
			Peralk_vals  = []
			for result_dict in NBO_vals:
				NBO_T_vals.append(result_dict["NBO/T"])
				Mg_num_vals.append(result_dict["Mg#"])
				Alk_vals.append(result_dict["Alkalinity Index"])
				Meta_vals.append(result_dict["Metaluminous"])
				Pera_vals.append(result_dict["Peraluminous"])
				Peralk_vals.append(result_dict["Peralkaline"])
			return_data["NBO_T"] = NBO_T_vals
			return_data["Mg#"] = Mg_num_vals
			return_data["Alkalinity Index"] = Alk_vals
			return_data["Metaluminous"] = Meta_vals
			return_data["Peraluminous"] = Pera_vals
			return_data["Peralkaline"] = Peralk_vals
		if metadata == False:
			return_data["NBO_T"] = NBO_vals

		return return_data

	def calc_gamma_Fe_metal(self, temperature, interactions=standard_interactions, print_warnings=False, print_status=True):
		"""
		Batch process calculation of gamma Fe. Calculates the activity coefficient, gamma, for iron. Interaction parameters epsilon are computed for all
		elements passed to interactions, so long as interaction parameter values are known.

		Parameters
		----------
		metal_comp: dict
			Dictionary of compositional information only for a metal, in terms of wt% elements

		temperature: float, or str
			Temperature, in degrees C. Can be passed as float, in which case the passed value is used as the temperature for all samples. 
			Alternatively, temperature information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object.

		interactions: list
			OPTIONAL.List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
			the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
			temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library. It is assumed
			that the temperature column is in the experimental_data portion of the file.
		
		print_warnings: bool
			OPTIONAL. Default is False. If set to True, any warnings related to the lack of compositional data or interaction
			parameters will be printed for each sample.

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame
		"""

		data = self.metal_data
		user_data = data.copy()

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		gamma_Fe_vals = []
		temp_vals = []
		iterno = 0
		for index, row in user_data.iterrows():
			iterno +=1
			if print_status == True:
				percent = iterno/len(data.index)
				status_bar(percent, index)
			if file_has_temp == True:
					temperature = self.experimental_data.loc[index][temp_name]
			bulk_comp = self.get_metal_comp(index)
			calc = calc_gamma_Fe_metal(metal_comp=bulk_comp, temperature=temperature, interactions=interactions, print_warnings=print_warnings)
			gamma_Fe_vals.append(calc)
			temp_vals.append(temperature)
		
		user_data["Temperature_C_Modeled"] = temp_vals
		user_data["gamma_Fe_metal"] = gamma_Fe_vals

		return_data = user_data.filter(['Temperature_C_Modeled', 'gamma_Fe_metal'], axis=1)

		return return_data

	def calc_gamma_solute_metal(self, temperature, species=None, interactions=standard_interactions, gammaFe="calculate", print_warnings=False, print_status=False):
		"""
		Calculates the activity coefficient, gamma, for all given solutes in an Fe-rich metal alloy. Interaction parameters epsilon are 
		computed for all elements passed to interactions, so long as interaction parameter values are known.

		Parameters
		----------
		temperature: float
			Temperature at which to perform the calculation, in degrees C.

		species: str or list
			List of elements for which to calculate gamma values, if that list is different than the list of interactions.
			Default value is None, in which case the elements list = interactions.

		interactions: list
			OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
			the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
			temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.
		
		gammaFe: dict
			OPTIONAL. Default is "calculate" in which case the gammaFe value will be calculated here.
			If the gammaFe value is already known for all samples, a dict or Pandas Series mapping the sample name
			to the gammaFe value can be passed here to avoid duplicating this calculation. Index of pandas Series or
			keys of the dict must match those of the user data (e.g. they must be sample names).

		print_warnings: bool
			OPTIONAL. Default is False. If set to True, any warnings related to the lack of compositional data or interaction
			parameters will be printed for each sample.

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame
			Dataframe with index as sample names plus the temperature used for modelling and all calculated gamma values
		"""
		data = self.metal_data
		user_data = data.copy()

		if isinstance(gammaFe, dict):
			has_gammaFe = True
		elif gammaFe == "calculate":
			gammaFe_single = gammaFe
			has_gammaFe = False
		else:
			raise InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		if species == None:
			species = interactions

		if isinstance(species, str):
			species = [species]

		return_data = pd.DataFrame(index=user_data.index)
		specno = 0
		iterno = 0
		for spec in species:
			specno += 1
			gamma_solute_vals = []
			for index, row in user_data.iterrows():
				iterno += 1
				percno = iterno*specno
				if print_status == True:
					length = len(data.index) * len(species)
					percent = iterno/length
					status_bar(percent, index)
				if file_has_temp == True:
					temperature = self.experimental_data.loc[index][temp_name]
				bulk_comp = self.get_metal_comp(index)
				if has_gammaFe == True:
					gammaFe_single = gammaFe[index]
				calc = calc_gamma_solute_metal(metal_comp=bulk_comp, species=spec, temperature=temperature, interactions=interactions, gammaFe=gammaFe_single, print_warnings=print_warnings)
				gamma_solute_vals.append(calc)
			column_head_title = "gamma_" + spec + "_metal"
			return_data[column_head_title] = gamma_solute_vals

		if file_has_temp == True:
			return_data["Temperature_C_Modeled"] = self.experimental_data[temp_name]
		else:
			return_data["Temperature_C_Modeled"] = [temperature]*len(return_data.index)

		return return_data

	def calc_gamma_FeO_silicate(self):
		"""
		Calculates the gammaFeO in the silicate melt for all samples.

		Returns
		-------
		pandas DataFrame
		"""

		data = self.silicate_data
		user_data = data.copy()

		gamma_FeO_vals = []
		for index, row in user_data.iterrows():
			bulk_comp = self.get_silicate_comp(index)
			calc = calc_gamma_FeO_silicate(silicate_comp=bulk_comp)
			gamma_FeO_vals.append(calc)

		return_data = pd.DataFrame(index=user_data.index)
		return_data["gamma_FeO_silicate"] = gamma_FeO_vals

		return return_data

	def metal_activity_from_composition(self, species, temperature, interactions=standard_interactions):
		"""
		Returns the activity of the given species in an Fe-rich metal.

		Parameters
		----------
		species: string or list
			Name of desired species for which to calculate the activity. Must match form of elements used in MetalSilicate (e.g., 'Fe',
			'W', 'Ti'). Can input multiple species as list of strings.

		temperature: float, or str
			Temperature, in degrees C. Can be passed as float, in which case the passed value is used as the temperature for all samples. 
			Alternatively, temperature information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object. Required to calculate gamma_Fe_metal.

		interactions: list
			OPTIONAL.List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
			the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
			temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.

		Returns
		-------
		pandas DataFrame
			Dataframe with index as sample names plus the temperature used for modelling and all calculated activity values
		"""
		data = self.metal_data
		user_data = data.copy()

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		if isinstance(species, str):
			species = [species] #if string, turn into list of length 1
		if isinstance(species, list):
			pass
		else:
			raise InputError("Species must be input as string or list of strings. You have input " + str(type(species)))

		return_data = pd.DataFrame(index=user_data.index)
		for s in species:
			activity_vals = []
			for index, row in user_data.iterrows():
				if file_has_temp == True:
					temperature = self.experimental_data.loc[index][temp_name]
				bulk_comp = self.get_metal_comp(index)
				try:
					calc = metal_activity_from_composition(metal_comp=bulk_comp, species=s, temperature=temperature, interactions=interactions)
					activity_vals.append(calc)
				except:
					activity_vals.append("Could not compute.")
					w.warn("Could not compute activity for " + str(s) + " for some samples.",RuntimeWarning,stacklevel=2)

			column_head_title = "a" + s
			return_data[column_head_title] = activity_vals

		if file_has_temp == True:
			return_data["Temperature_C_Modeled"] = self.experimental_data[temp_name]
		else:
			return_data["Temperature_C_Modeled"] = [temperature]*len(return_data.index)

		return return_data

	def calc_dIW(self, temperature, interactions=standard_interactions, gammaFe="calculate", gammaFeO="calculate", print_warnings=False, print_status=True):
		"""
		Calculates fO2 in terms of delta Iron-Wustite for all samples. Calculation is performed using mole fractions 
		and activity coefficients of Fe in the metal and FeO in the silicate

		Parameters
		----------
		temperature: float, or str
			Temperature, in degrees C. Can be passed as float, in which case the passed value is used as the temperature for all samples. 
			Alternatively, temperature information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object. Required to calculate gamma_Fe_metal.

		interactions: list
			OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
			the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
			temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.

		gammaFe and gammaFeO: dict
			OPTIONAL. Default is "calculate" in which case the gammaFe or gammaFeO value will be calculated here.
			If the gammaFe or gammaFeO value is already known for all samples, a dict or Pandas Series mapping the sample name
			to the gammaFe or gammaFeO value can be passed here to avoid duplicating this calculation. Index of pandas Series or
			keys of the dict must match those of the user data (e.g. they must be sample names).

		print_warnings: bool
			OPTIONAL. Default is False. If set to True, any warnings related to the lack of compositional data or interaction
			parameters will be printed for each sample.

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame
		"""
		silicate_data = self.silicate_data.copy()
		metal_data = self.metal_data.copy()

		if isinstance(gammaFe, dict):
			has_gammaFe = True
		elif gammaFe == "calculate":
			gammaFe_single = gammaFe
			has_gammaFe = False
		else:
			raise InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")

		if isinstance(gammaFeO, dict):
			has_gammaFeO = True
		elif gammaFeO == "calculate":
			gammaFeO_single = gammaFeO
			has_gammaFeO = False
		else:
			raise InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")


		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temp must be type str or float or int")

		temp_vals = []
		fO2_vals = []
		iterno = 0
		for index, row in silicate_data.iterrows():
			iterno += 1
			if print_status == True:
				percent = iterno/len(silicate_data.index)
				status_bar(percent, index)
			if file_has_temp == True:
				temperature = self.experimental_data.loc[index][temp_name]
			silicate_comp = self.get_silicate_comp(index)
			metal_comp = self.get_metal_comp(index)
			if has_gammaFe == True:
				gammaFe_single = gammaFe[index]
			if has_gammaFeO == True:
				gammaFeO_single = gammaFeO[index]
			calc = calc_dIW(silicate_comp=silicate_comp, metal_comp=metal_comp, temperature=temperature, interactions=interactions, gammaFe=gammaFe_single, gammaFeO=gammaFeO_single, print_warnings=print_warnings)
			fO2_vals.append(calc)

		return_data = pd.DataFrame(index=silicate_data.index)
		return_data["fO2_dIW"] = fO2_vals

		return return_data

	def calc_logfO2(self, pressure, temperature, dIW="calculate", print_status=False, **kwargs):
		"""
		Calculates the absolute fO2 value (as log(fO2)) based on deltaIW value, pressure, and temperature.

		Parameters
		----------
		pressure: float or str
			Pressure in GPa. Can be passed as float, in which case the passed value is used as the pressure for all samples. 
			Alternatively, pressure information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object.

		temperature: float or str
			Temperature in degrees C. Can be passed as float, in which case the passed value is used as the temperature for all samples. 
			Alternatively, temperature information for each individual sample may already be present in the ExcelFile object. If so, pass 
			the str value corresponding to the columntitle in the ExcelFile object.

		dIW: float or str or pandas Series
			fO2 in terms of deltaIW. If "calculate" is passed, the dIW value will be calculated using the silicate and metal information 
			given. If "calculate" is passed, user may also pass any arguments relevant to calc_dIW() method to **kwargs. dIW can be passed 
			as float, in which case the passed value is used as the dIW for all samples. Alternatively, dIW information for each individual 
			sample may already be present in the ExcelFile object. If so, pass the str value corresponding to the columntitle in the 
			ExcelFile object.

		print_status: bool
			OPTIONAL. Default is True. If set to True, the progress of the calculation will be printed to the terminal.

		Returns
		-------
		pandas DataFrame
			Column name for computed values is "log(fO2)"
		"""
		experimental_data = self.experimental_data.copy()
		dIW_to_be_calcd = False
		dIW_is_series = False
		file_has_dIW = False
		file_has_press = False
		file_has_temp = False

		if isinstance(pressure, str):
			file_has_press = True
			press_name = pressure
		elif isinstance(pressure, float) or isinstance(pressure, int):
			file_has_press = False
		else:
			raise InputError("pressure must be type str or float or int")

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise InputError("temperature must be type str or float or int")

		if isinstance(dIW, str):
			if dIW == "calculate":
				file_has_dIW = False
				dIW_to_be_calcd = True
				if 'interactions' not in kwargs:
					interactions=standard_interactions
				if 'gammaFe' not in kwargs:
					gammaFe = "calculate"
				if 'gammaFeO' not in kwargs:
					gammaFeO = "calculate"
				if 'print_warnings' not in kwargs:
					print_warnings=False
				if print_status == True:
					print("Calculating dIW values...")
				dIW_calcd_vals = self.calc_dIW(temperature=temperature, interactions=interactions, gammaFe=gammaFe, gammaFeO=gammaFeO, print_warnings=print_warnings, print_status=print_status)
			else:
				file_has_dIW = True
				dIW_name = dIW
		elif isinstance(dIW, float) or isinstance(dIW, int):
			file_has_dIW = False
		elif isinstance(dIW, pd.core.series.Series):
			dIW_is_series = True
			dIW_series = dIW
		else:
			raise InputError("dIW must be type str or float or int or pandas Series")

		logfO2_vals = []
		dIW_calcd_list = []
		iterno = 0
		if print_status == True and dIW_to_be_calcd == True:
			print("Calculating log(fO2) values...")
		for index, row in experimental_data.iterrows():
			iterno += 1
			if print_status == True:
				percent = iterno/len(experimental_data.index)
				status_bar(percent, index)
			if file_has_press == True:
				pressure = experimental_data.loc[index][press_name]
			if file_has_temp == True:
				temperature = experimental_data.loc[index][temp_name]
			if file_has_dIW == True:
				dIW = experimental_data.loc[index][dIW_name]
			if dIW_to_be_calcd == True:
				dIW = dIW_calcd_vals.loc[index]["fO2_dIW"]
				dIW_calcd_list.append(dIW)
			if dIW_is_series == True:
				dIW = dIW_series.loc[index]

			calc = calc_logfO2_from_IW(pressure=pressure, temperature=temperature, dIW=dIW)
			logfO2_vals.append(calc)

		return_data = pd.DataFrame(index=experimental_data.index)
		return_data["log(fO2)"] = logfO2_vals
		if dIW_to_be_calcd == True:
			return_data["fO2_dIW"] = dIW_calcd_list

		return return_data

	def calc_optical_basicity(self):
		"""
		Calculates the optical basicity foe every sample based on composition. Uses optical
		basicity values from Table 4 of Leboutellier A, Courtine P (1998) Improvement of a
		bulk optical basicity table for oxidic systems. J Solid State Chem 137:94–103.
		https://doi.org/10.1006/jssc.1997.7722

		Uses method after Crisp, L.J. and Berry, A.J. (2022) A new model for zircon saturaiton in
		silicate melts. Contributions to Mineralogy and Petrology, 177:71.
		https://doi.org/10.1007/s00410-022-01925-6

		Lambda = Sigma(m_i*n_i*Lambda_i)/Sigma(m_i*n_i)
		where m_i is the number of moles of each oxide, n_i is the number of oxygens in the oxide
		(for example 3 for Al2O3), and Lambda_i is the optical basicity coefficient of each oxide
		from Table 4 of Leboutellier and Courtine (1998). The optical basicity for all Fe-bearing
		melt compositions here is calculated assuming Fe 3+ /ΣFe = 0, (where ΣFe = Fe 2+ + Fe 3+ ).

		Parameters
		----------
		self: ExcelFile object

		Returns
		-------
		pandas DataFrame
			Column name for computed values is "Lambda"	
		"""
		silicate_data = self.silicate_data.copy()

		optical_bas_list = []
		for index, row in silicate_data.iterrows():
			silicate_comp = self.get_silicate_comp(index)
			calc = calc_optical_basicity(silicate_comp)
			optical_bas_list.append(calc)

		return_data = pd.DataFrame(index=silicate_data.index)
		return_data["Lambda"] = optical_bas_list

		return return_data

#----------PREP DATA---------#
#Write these to perform calculation on one sample at a time, use ExcelFile class methods to do iteration over the entire dataset.

def calc_ln_gamma_naught_at_temperature(species, temperature):
	"""
	Calculates the reference value for the activity coefficient gamma at the given temperature based on a known value for gamma at a reference 
	temperature. NOTA BENE: if there is no tabulated value for the reference gamma, ideality will be assumed (the reference gamma will be
	set equal to 1).

	Parameters
	----------
	species: str
		String of the name of the element for which to calculate gamma_naught

	temperature: float
		Temperature at which to calculate the reference gamma, in degrees C.
	"""
	user_T_K = temperature + 273.15

	try:
		ref_gamma = ac.activity_coeffs.loc[species]['gamma_naught_i']
		ref_T_K = ac.activity_coeffs.loc[species]['Temp_K']

		if isinstance(ref_gamma, float) == False and isinstance(ref_gamma, int) == False:
			ln_gamma_naught = 1
		else:
			ln_gamma_naught = (ref_T_K/user_T_K)*np.log(ref_gamma)
	except:
		ln_gamma_naught = 1

	return ln_gamma_naught

def calc_epsilon_at_temperature(i, j, temperature):
	"""
	Calculates the value for the interaction parameter epsilon between elements i and j at the given temperature based on a known value for e at a 
	reference temperature.

	Parameters
	----------
	i: str
		String of the name of the first of two elements for which to calculate epsilon

	j: str
		String of the name of the second of two elements for which to calculate epsilon

	temperature: float
		Temperature at which to calculate the reference gamma, in degrees C.
	"""
	i_j = i + ' ' + j
	ref_e = ip.interaction_params.loc[i_j]['ei(j,k)']

	#convert e_i_j to epsilon_i_j using Equation 20 in Steelmaking Data Sourcebook
	epsilon_i_j_ref = (ref_e/0.00434 - 1)*(element(j).atomic_weight/55.85)+1

	#Convert epsilon to value at sample temp using equation 14 from Corgne et al. (2008)
	ref_T_K = ip.interaction_params.loc[i_j]['Temp_K']
	user_T_K = temperature + 273.15
	epsilon_i_j = (ref_T_K/user_T_K)*epsilon_i_j_ref

	return epsilon_i_j

def calc_gamma_Fe_metal(metal_comp, temperature, interactions=standard_interactions, print_warnings=True):
	"""
	Calculates the activity coefficient, gamma, for iron. Interaction parameters epsilon are computed for all
	elements passed to interactions, so long as interaction parameter values are known.

	Parameters
	----------
	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
		the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
		temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.
	
	print_warnings: bool
			OPTIONAL. Default is True. If set to True, any warnings related to the lack of compositional data or interaction
			parameters will be printed.
	"""
	# Note: Exclude Fe from this calculation by default

	def calc_indiv_term1(i, molfrac, epsilon_i_i):
		return epsilon_i_i*(molfrac[i] + np.log(1-molfrac[i]))

	def calc_indiv_term2(j, k, molfrac, epsilon_j_k):
		return epsilon_j_k*molfrac[j]*molfrac[k]*(1+((np.log(1-molfrac[j]))/molfrac[j]) + (np.log(1-molfrac[k])/molfrac[k]))

	def calc_indiv_term3(i, k, molfrac, epsilon_i_k):
		return epsilon_i_k*molfrac[i]*molfrac[k]*(1+((np.log(1-molfrac[k]))/molfrac[k]) - (1/(1-molfrac[i])))

	def calc_indiv_term4(j, k, molfrac, epsilon_j_k):
		return epsilon_j_k*molfrac[j]**2*molfrac[k]**2*((1/(1-molfrac[j])) + (1/(1-molfrac[k])) - 1)

	def calc_indiv_term5(i, k, molfrac, epsilon_i_k):
		return epsilon_i_k*molfrac[i]**2*molfrac[k]**2*((1/(1-molfrac[i])) + (1/(1-molfrac[k])) + (molfrac[i]/(2*(1-molfrac[i])**2)) - 1)

	user_interactions = interactions.copy()
	interactions = []
	skipped_elements = []
	for element in user_interactions:
		if element in metal_comp and metal_comp[element]>0:
			interactions.append(element)
		else:
			skipped_elements.append(element)

	
	molfrac = wtpercentElements_to_molElements(metal_comp)

	interactions_NminusOne = interactions[:-1]
	interactions_JplusOne = interactions[1:]

	term1 = 0
	term2 = 0
	term3 = 0
	term4 = 0
	term5 = 0

	no_interaction_param = []
	for i in interactions:
		skip_this = False
		try:
			epsilon_i_i = calc_epsilon_at_temperature(i, i, temperature)
		except:
			skip_this = True
			param_string = str(i) + " and " + str(i)
			no_interaction_param.append(param_string)
		if skip_this == False:
			term1 += calc_indiv_term1(i, molfrac, epsilon_i_i)

	for j in interactions_NminusOne:
		for k in interactions_JplusOne:
			skip_this = False
			try:
				epsilon_j_k = calc_epsilon_at_temperature(j, k, temperature)
			except:
				skip_this = True
				param_string = str(j) + " and " + str(k)
				no_interaction_param.append(param_string)
			if skip_this == False:
				term2 += calc_indiv_term2(j, k, molfrac, epsilon_j_k)
				term4 += calc_indiv_term4(j, k, molfrac, epsilon_j_k)
			

	for i in interactions:
		for k in interactions:
			skip_this = False
			if i == k:
				pass
			else:
				try:
					epsilon_i_k = calc_epsilon_at_temperature(i, k, temperature)
				except:
					skip_this = True
					param_string = str(i) + " and " + str(k)
					no_interaction_param.append(param_string)
				if skip_this == False:
					term3 += calc_indiv_term3(i, k, molfrac, epsilon_i_k)
					term5 += calc_indiv_term5(i, k, molfrac, epsilon_i_k)

	ln_gamma_Fe = term1 - term2 + term3 + 0.5*term4 - term5
	gamma_Fe = np.exp(ln_gamma_Fe)

	if print_warnings == True:
		if len(skipped_elements) > 0:
			w.warn("No compositional information given for: " + ', '.join(skipped_elements) + ". These elements were skipped.",RuntimeWarning,stacklevel=2)
		if len(no_interaction_param) > 0:
			w.warn("No interaction parameters known for: " + ', '.join(no_interaction_param),RuntimeWarning,stacklevel=2)

	return gamma_Fe

def calc_gamma_solute_metal(metal_comp, species, temperature, gammaFe="calculate", interactions=standard_interactions, print_warnings=True, **kwargs):
	"""
	Calculates the activity coefficient, gamma, for any solutes in an Fe-rich metal alloy. Interaction parameters epsilon are computed for all
	elements passed to interactions, so long as interaction parameter values are known.

	Parameters
	----------
	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	species: str
		String of the name of the element for which to calculate gamma

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	gammaFe: float
		OPTIONAL. Default is "calculate" in which case the gammaFe value will be calculated here.
		If the gammaFe value is already known, it can be passed here to avoid having to calculate it again.

	elements: list
		OPTIONAL. List of elements for which to calculate gamma values, if that list is different than the list of interactions.
		Default value is None, in which case the elements list = interactions.

	interactions: list
		OPTIONAL.List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
		the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
		temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.
	
	print_warnings: bool
			OPTIONAL. Default is True. If set to True, any warnings related to the lack of compositional data or interaction
			parameters will be printed.
	"""
	# Note: Should Fe be included in the interactions for this calculation?

	user_interactions = interactions.copy()
	interactions = []
	no_comp_interact = []
	no_interaction_param = []
	for element in user_interactions:
		if element in metal_comp and metal_comp[element]>0:
			interactions.append(element)
		else:
			no_comp_interact.append(element)

	if gammaFe == "calculate":
		gammaFe = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature, interactions=interactions, print_warnings=print_warnings)
	ln_gamma_Fe_metal = np.log(gammaFe)
	
	i = species
	if metal_comp[i]>0:	
		ln_gamma_species_naught_ref = calc_ln_gamma_naught_at_temperature(species=species, temperature=temperature)
		molfrac = wtpercentElements_to_molElements(metal_comp)
		try:
			epsilon_i_i = calc_epsilon_at_temperature(i=species, j=species, temperature=temperature)
		except:
			epsilon_i_i = 0 #TODO should this be 1?

		term1 = 0
		term2 = 0

		for j in interactions:
			skip_this = False
			if j == i:
				pass
			else:
				try:
					epsilon_i_j = calc_epsilon_at_temperature(i=i, j=j, temperature=temperature)
				except:
					skip_this = True
					param_string = str(i) + " and " + str(j)
					no_interaction_param.append(param_string)

				if skip_this == False:
					term1 += epsilon_i_j*molfrac[j]*(1 + ((np.log(1-molfrac[j]))/molfrac[j]) - (1/(1-molfrac[i])))
					term2 += epsilon_i_j*molfrac[j]**2*molfrac[i]*((1/(1-molfrac[i])) + (1/(1-molfrac[j])) + (molfrac[i]/(2*(1-molfrac[i])**2))-1)

		ln_gamma_i = ln_gamma_Fe_metal + ln_gamma_species_naught_ref - epsilon_i_i*np.log(1-molfrac[species]) - term1 + term2

		gamma_i = np.exp(ln_gamma_i)

		if print_warnings == True:
			if len(no_comp_interact) > 0:
				w.warn("No compositional information given for: " + ', '.join(no_comp_interact) + ". Interaction with this element not calculated.",RuntimeWarning,stacklevel=2)
			if len(no_interaction_param) > 0:
				w.warn("No interaction parameter known for: " + ', '.join(no_interaction_param) + ".")

		return gamma_i
	else:
		if print_warnings == True:
			w.warn("No compositional information given for species " + species + ". Cannot calculate gamma without concentration.",RuntimeWarning,stacklevel=2)
		return 0

def calc_gamma_FeO_silicate(silicate_comp):
	"""
	Returns a value for gammaFeO in the silicate. Parameterization is based on Holdzheid, where gammaFeO is taken as a constant value from
	1.7-3, dependent only upon MgO content. 

	Parameters
	----------
	silicate_comp: dict
		Dictionary of the composition of the silicate in wt% oxides.

	Returns
	-------
	float
		gammaFeO in the silicate melt
	"""

	if silicate_comp['MgO'] <= 20:
		return 1.7
	else:
		gamma_value =  1.7 + 0.1*(silicate_comp['MgO']-20)
		if gamma_value <= 3:
			return gamma_value
		else:
			return 3

# def calc_gamma_FeO_silicate_test(silicate_comp):
# 	"""
# 	Returns a value for gammaFeO in the silicate. Parameterization is based on Holdzheid, where gammaFeO is taken as a constant value from
# 	1.7-3, dependent only upon MgO content. 

# 	Parameters
# 	----------
# 	silicate_comp: dict
# 		Dictionary of the composition of the silicate in wt% oxides.

# 	Returns
# 	-------
# 	float
# 		gammaFeO in the silicate melt
# 	"""
# 	gamma_value =  1.7 + 0.1*(silicate_comp['MgO']-20)
# 	return gamma_value


def calc_activity(X, gamma):
	"""
	Returns the value of the activity of any species given the concentration of that species in mol fraction, X, and the
	activity coefficient gamma.

	Parameters
	----------
	X: float
		Concentration of the given species in mol fraction

	gamma: float
		Activity coefficient of the given species

	Returns
	-------
	float
		Activity of the given species
	"""

	return X * gamma

def metal_activity_from_composition(metal_comp, species, temperature, interactions=standard_interactions):
	"""
	Returns the activity of the given species in an Fe-rich metal, calculated as X times gamma.

	Parameters
	----------
	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	species: string
		Name of desired species for which to calculate the activity. Must match form of elements used in MetalSilicate (e.g., 'Fe',
		'W', 'Ti')

	temperature: float
		Temperature at which to perform the calculation, in degrees C.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
		the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
		temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.

	Returns
	-------
	float
		Activity of the given species in an Fe-rich metal
	"""
	if isinstance(species, str):
		pass
	else:
		raise InputError("Species must be passed as string. You have passed " + str(type(species)))

	molfrac_metal = wtpercentElements_to_molElements(metal_comp)
	try:
		X = molfrac_metal[species]
	except:
		raise InputError("No compositional information given for " + species)

	if X == 0:
		raise InputError("No compositional information given for " + species)

	if species == "Fe":
		gamma = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature, interactions=interactions, print_warnings=False)
	elif species in elements:
		gamma = calc_gamma_solute_metal(metal_comp=metal_comp, species=species, temperature=temperature, 
										gammaFe="calculate", interactions=standard_interactions, print_warnings=False)
	else:
		raise InputError("Species must be one of: " + str(elements))

	if gamma == 0:
		raise GeneralError("Something went wrong. Unable to calculate an activity coefficient gamma for " + species)

	return calc_activity(X=X, gamma=gamma)

def calc_dIW(silicate_comp, metal_comp, temperature=None, gammaFe="calculate", gammaFeO="calculate", interactions=standard_interactions, print_warnings=False):
	"""
	Returns fO2 in terms of delta Iron-Wustite. Calculation is performed using mole fractions and activity coefficients
	of Fe in the metal and FeO in the silicate

	Parameters
	----------
	silicate_comp: dict
		Dictionary of the composition of the silicate in wt% oxides.

	metal_comp: dict
		Dictionary of compositional information only for a metal, in terms of wt% elements

	temperature: float
		Temperature in degrees C.

	gammaFe and gammaFeO: float
		OPTIONAL. Default is "calculate" in which case the gammaFe or gammaFeO value will be calculated here.
		If the gammaFe or gammaFeO value is already known, it can be passed here to avoid having to calculate it again.

	interactions: list
		OPTIONAL. List of strings of element names. Elements are solutes in a metal Fe liquid alloy for which interaction parameters are known and for which
		the user wishes to calculate the effects of interaction within the alloy. Elements need not be infinitely dilute. Compositional and
		temperature ranges for which interaction parameters are known are given in the interaction_parameters script within this library.

	Returns
	-------
	float
		logfO2 in terms of delta Iron-Wustite
	"""
	molfrac_silicate = wtpercentOxides_to_molOxides(silicate_comp)
	molfrac_metal = wtpercentElements_to_molElements(metal_comp)

	X_FeO_silicate = molfrac_silicate['FeO']
	X_Fe_metal = molfrac_metal['Fe']

	if gammaFeO == "calculate":
		gamma_FeO_silicate = calc_gamma_FeO_silicate(silicate_comp)
	else:
		gamma_FeO_silicate = gammaFeO

	if gammaFe == "calculate":
		if temperature == None:
			raise InputError("Error: no temperature given. Cannot calculate gammaFe without temperature.")
		gamma_Fe_metal = calc_gamma_Fe_metal(metal_comp=metal_comp, temperature=temperature, interactions=interactions, print_warnings=print_warnings)
	else:
		gamma_Fe_metal = gammaFe

	return 2*np.log10(X_FeO_silicate/X_Fe_metal) + 2*np.log10(gamma_FeO_silicate/gamma_Fe_metal)

def calc_IW(pressure, temperature):
	"""
	Fe-FeO (Iron-Wustite)
	=====================
	Define IW buffer value at P
	
	References
	----------
	Campbell et al. (2009) High-pressure effects on the iron-iron oxide and nickel-nickel oxide oxygen fugacity buffers

	Parameters
	----------
	pressure: float
		Pressure in GPa

	temperature: float or numpy array
		Temperature in degrees C

	Returns
	-------
	float or numpy array
		log(fO2)

	Polynomial coefficients
	-----------------------
	log fO2  =  (a0+a1*P) + (b0+b1*P+b2*P^2+b3*P^3)/T
	a0: 6.54106 
	a1: 0.0012324
	b0: -28163.6
	b1: 546.32
	b2: -1.13412
	b3: 0.0019274               
	"""
	#convert T from C to K
	T_K = temperature+273.15

	log_fO2 = (6.54106+0.0012324*pressure) + (-28163.6+546.32*pressure-1.13412*pressure**2+0.0019274*pressure**3)/T_K

	return log_fO2

def calc_logfO2_from_IW(pressure, temperature, dIW="calculate"):
	"""
	Calculates the absolute fO2 value (as log(fO2)) based on deltaIW value, pressure, and temperature.

	Parameters
	----------
	pressure: float
		Pressure in GPa

	temperature: float
		Temperature in degrees C

	dIW: float
		fO2 in terms of deltaIW

	Returns
	-------
	float
		log(fO2) absolute value
	"""

	IW_at_PT = calc_IW(pressure,temperature)
	log_fO2 = IW_at_PT + dIW

	return log_fO2

def calc_NBO_T(sample, method=4, metadata=False):
	"""
	Calcultes the NBO/T (non-bridging oxygens over tetrahedrally coordinated ions) for a given composition. Written after
	excel spreadsheet by FM McCubbin.

	Parameters
	----------
	sample: dict
		Dictionary of compositional information only, in terms of wt%

	method: int
		Default is 4. Integer from 1 to 4 inclusive. Use to choose which NBO/T parameterization to use as defined here:
			1: NF = Si, Al-mixed, Ti
			2: NF = Si, Al-mixed, Ti, P
			3: NF = Si, Al-mixed, Fe3+, Ti, P
			4: NF = Si, Al-mixed, Fe3+, Cr3+, Ti, P

	metadata: bool
		Default is False. If set to true, a dict is returned with more information calculated about the sample.

	Returns
	-------
	float or dict
		If metadata=False (default), returns float of NBO/T value
		If metadata=True, returns dict with keys #TODO add keys here
	"""
	oxides_for_calc = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "Fe2O3", "FeO", "MgO", "CaO",
						"Na2O", "K2O", "MnO", "P2O5"]

	_sample = get_oxides(sample) #only get oxide data, nothing extranneous
	_sample = {oxide:  sample[oxide] for oxide in oxides_for_calc} #only get oxides we care about for this calc

	#STEP 1: Convert to mole fraction
	sample_X = wtpercentOxides_to_molOxides(_sample)

	#STEP 2: Convert wt% to chemical forumla on a single oxygen basis
	atoms = wtpercentOxides_to_molSingleO(_sample)

	#CALCULATE METADATA
	Mg_number = (100*atoms["Mg"]/(atoms["Mg"]+atoms["Fe"]))
	alkalinity_index = (atoms["Na"]+atoms["K"])/atoms["Al"]

	if atoms["Al"] < (atoms["Ca"]+atoms["Na"]+atoms["K"]) and atoms["Al"] > (atoms["Na"]+atoms["K"]):
		metaluminous = True
	else:
		metaluminous = False

	if atoms["Al"] > atoms["Ca"] + atoms["Na"] + atoms["K"]:
		peraluminous = True
	else:
		peraluminous = False

	if atoms["Al"] < atoms["Na"] + atoms["K"]:
		peralkaline = True
	else:
		peralkaline = False

	#STEP 3: Calculate NBO/T
	sum_atoms = sum(atoms.values())

	if method == 3:
		atoms["Cr"] = 0
	if method == 2:
		atoms["Cr"] = 0
		atoms["Fe3"] = 0
	if method == 1:
		atoms["P"] = 0

	a_alpha = (atoms["Fe3"]+atoms["Al"]+atoms["Cr"])
	a_beta = (atoms["K"]+atoms["Na"]+0.5*(atoms["Fe"]+atoms["Mg"]+atoms["Mn"] + (atoms["Ca"]-1.5*atoms["P"])))
	if a_alpha > a_beta:
		a = a_beta
		b = a_alpha - a_beta
	else:
		a = a_alpha
		b = 0

	if sum_atoms < a_alpha:
		NBO_over_T = (1/(atoms["Si"]+atoms["Ti"]+atoms["P"]+(a*b*3+atoms["Fe"]*2+atoms["Mg"]*2 + 
						(atoms["Ca"]-(atoms["P"]*1.5)-(0.5*((a_alpha)-atoms["Na"]-atoms["K"])))*2+atoms["Mn"]*2)))
	else:
		NBO_over_T = (1/(atoms["Si"]+atoms["Ti"]+atoms["Al"]+atoms["Fe3"]+atoms["P"]+atoms["Cr"])*
						(atoms["Fe"]*2+atoms["Mg"]*2+(atoms["Ca"]-(atoms["P"]*1.5))*2+
						(atoms["Na"]+atoms["K"]-atoms["Al"]-atoms["Fe3"]-atoms["Cr"])+atoms["Mn"]*2))

	if metadata == True:
		return {"NBO/T": NBO_over_T,
				"Mg#": Mg_number,
				"Alkalinity Index": alkalinity_index,
				"Metaluminous": metaluminous,
				"Peraluminous": peraluminous,
				"Peralkaline": peralkaline}
	else:
		return NBO_over_T


def calc_optical_basicity(sample):
	"""
	Calculates the optical basicity of a sample based on composition. Uses optical basicity values
	from Table 4 of Leboutellier A, Courtine P (1998) Improvement of a bulk optical basicity table
	for oxidic systems. J Solid State Chem 137:94–103. https://doi.org/10.1006/jssc.1997.7722

	Uses method after Crisp, L.J. and Berry, A.J. (2022) A new model for zircon saturaiton in
	silicate melts. Contributions to Mineralogy and Petrology, 177:71.
	https://doi.org/10.1007/s00410-022-01925-6

	Lambda = Sigma(m_i*n_i*Lambda_i)/Sigma(m_i*n_i)
	where m_i is the number of moles of each oxide, n_i is the number of oxygens in the oxide (for
	example 3 for Al2O3), and Lambda_i is the optical basicity coefficient of each oxide from
	Table 4 of Leboutellier and Courtine (1998). The optical basicity for all Fe-bearing melt
	compositions here is calculated assuming Fe 3+ /ΣFe = 0, (where ΣFe = Fe 2+ + Fe 3+ ).

	Parameters
	----------
	sample: dict
		Dictionary of compositional information only, in terms of wt%

	Returns
	-------
	float
		Optical basicity value, Lambda.	
	"""
	oxides_for_calc = ob.optical_basicity.index.tolist()

	# only get oxide data, nothing extranneous
	_sample = get_oxides(sample)

	# STEP 1: Convert to mole fraction
	sample_X = wtpercentOxides_to_molOxides(_sample)

	# STEP 2: Calculate numerator
	numer_sum = 0
	for ox, val in sample_X.items():
		if ox in oxides_for_calc:
			numer = val * OxygenNum[ox] * ob.optical_basicity.loc[ox]['Lambda']
			numer_sum += numer
		else:
			pass

	# STEP 3: Calculate denominator
	denom_sum = 0
	for ox, val in sample_X.items():
		if ox in oxides_for_calc:
			denom = val * OxygenNum[ox]
			denom_sum += denom
		else:
			pass

	# STEP 4: divide
	optical_bas = numer_sum / denom_sum

	return optical_bas
