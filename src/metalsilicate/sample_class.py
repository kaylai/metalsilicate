from metalsilicate import core, sample_class, batchfile
import pandas as pd
import numpy as np

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
	for oxide in core.oxides:
		if oxide in _sample.keys():
			pass
		else:
			_sample[oxide] = 0.0

	clean = {oxide:  _sample[oxide] for oxide in core.oxides}

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
	for element in core.elements:
		if element in _sample.keys():
			pass
		else:
			_sample[element] = 0.0

	clean = {element:  _sample[element] for element in core.elements}

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
		multi_sample["Sum"] = sum([multi_sample[oxide] for oxide in core.oxides])
		for column in multi_sample:
			if column in core.oxides:
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
	elif isinstance(sample, batchfile.ExcelFile):
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
			molElements[element] = sample_elements[element]/core.elementMass[element]

		if type(sample_elements) == pd.core.series.Series:
			molElements = pd.Series(molElements)
			molElements = molElements/molElements.sum()
		else:
			total = np.sum(list(molElements.values()))
			for element in elementslist:
				molElements[element] = molElements[element]/total

		return molElements

	elif isinstance(sample_class.sample, pd.DataFrame):
		data = sample_class.sample
		for key, value in core.elementMass.items():
			data.loc[:, key] /= value

		data["MPSum"] = sum([data[element] for element in sample_elements])

		for element in sample_elements:
			data.loc[:, element] /= data['MPSum']
		del data['MPSum']

		return data

	else:
		raise core.InputError("The composition input must be a pandas Series or dictionary.")

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
			molOxides[ox] = sample_oxides[ox]/core.oxideMass[ox]

		if type(sample_oxides) == pd.core.series.Series:
			molOxides = pd.Series(molOxides)
			molOxides = molOxides/molOxides.sum()
		else:
			total = np.sum(list(molOxides.values()))
			for ox in oxideslist:
				molOxides[ox] = molOxides[ox]/total

		return molOxides

	elif isinstance(sample_class.sample, pd.DataFrame):
		data = sample_class.sample
		for key, value in core.oxideMass.items():
			data.loc[:, key] /= value

		data["MPOSum"] = sum([data[oxide] for oxide in sample_oxides])

		for oxide in sample_oxides:
			data.loc[:, oxide] /= data['MPOSum']
		del data['MPOSum']

		return data

	else:
		raise core.InputError("The composition input must be a pandas Series or dictionary.")

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
		raise core.InputError("The composition input must be a pandas Series or dictionary.")

	total_O = 0.0
	for ox in oxideslist:
		if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
			cation = sample_class.oxides_to_cations[ox]
			molCations[cation] = core.CationNum[ox]*sample_oxides[ox]/core.oxideMass[ox]
			total_O += core.OxygenNum[ox]*sample_oxides[ox]/core.oxideMass[ox]
	if type(sample_oxides) == pd.core.series.Series:
		molCations = pd.Series(molCations)
		molCations = molCations/total_O
	else:
		# total = np.sum(list(molCations.values()))
		for ox in oxideslist:
			if exclude_volatiles == False or (ox != 'H2O' and ox != 'CO2'):
				cation = sample_class.oxides_to_cations[ox]
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
		element = sample_class.oxides_to_cations[oxide]
		conversion_factor = core.CationNum[oxide]*core.elementMass[element]/core.oxideMass[oxide]
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
		oxide = sample_class.cations_to_oxides[element]
		conversion_factor = core.CationNum[oxide]*core.elementMass[element]/core.oxideMass[oxide]
		converted_sample[oxide] = sample_elements[element] / conversion_factor

	return converted_sample
