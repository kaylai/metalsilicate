import pandas as pd
import warnings as w
import sys

from metalsilicate import core
from metalsilicate import sample_class
from metalsilicate import calculate

#----------IMPORT FILE AND PROCESS---------#
def preprocess_dataset(data, how):
		"""
		Adds 0.0 values to any oxide or element data not passed.

		Parameters
		----------
		data: pandas DataFrame
			Composition of samples in wt%

		how: str
			Describes if data are in terms of oxides or elements. Can be one of 'oxides',
			'elements'.

		Returns
		-------
		pandas DataFrame
		"""
		_data = data.copy()

		if how == "oxides":
			datatype = core.oxides
		if how == "elements":
			datatype = core.elements

		for item in datatype:
			if item in _data.columns:
				pass
			else:
				_data[item] = 0.0

		return _data

def clean(data):
	"""
	Takes a pandas dataframe (e.g. myfile.data, myfile.silicate_data) and removes any columns with
	all 0's, any non-numeric data.

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
	"""
	Looks for duplicate names in the index of a dataframe. If found, each duplicate name is
	appended with '-duplicate-n' where n is the nth duplicate.
	
	Parameters
	----------
	df: pandas DataFrame
        A pandas DataFrame object with an index of strings
		
	suffix: a string to append to the name of each duplicate string in the index. A number 'n' will
	always be appended after suffix indicating the nth duplicate
	
	Returns
	-------
	pandas DataFrame
	"""
	appendents = (suffix +
			      df.groupby(level=0).cumcount().astype(str).replace('0','')).replace(suffix, '')
	return df.set_index(df.index.astype(str) + appendents)

class ExcelFile(object):
	"""An excel file with sample names and variables. File must have two header rows. The first
	header row breaks data into the following: data about the experiments themselves (e.g.,
	pressure, temperature, duration, experimental apparatus); data describing the silicate melt
	composition in wt% oxides; data describing the first metal composition in wt% elements; data
	describing a second metal composition in wt% elements; any notes about the sample. The first
	column will become the index of all of these datasets, so by convention this should be the
	sample name.

	Attributes
	----------
		filename: str
			Path to the excel file, e.g., "my_file.xlsx". This always needs to be passed, even if
			the user is passing a pandas DataFrame rather than an Excel file.	

		sheet_name: str
			OPTIONAL. Default value is 0 which gets the first sheet in the excel spreadsheet file.
			This implements the pandas. read_excel() sheet_name parameter. But functionality to read
			in more than one sheet at a time (e.g., pandas.read_excel(sheet_name=None)) is not yet
			imlpemented in VESIcal. From the pandas 1.0.4 documentation:
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
					metal_data="Metal_1_wt_elements", metal_2_data="Metal_2_wt_elements",
					notes="Notes"):
		"""Return an ExcelFile object whoes parameters are defined here."""

		if isinstance(sheet_name, str) or isinstance(sheet_name, int):
			pass
		else:
			raise core.InputError("If sheet_name is passed, it must be of type str or int. Cannot" \
			"import more than one sheet at a time.")

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
		Saves data calculated by the user in batch processing mode (using the ExcelFile class
		methods) to an organized excel file, with the original user data plus any calculated data.

		Parameters
		----------
		filename: string
			Name of the file. Extension (.xlsx) should be passed along with the name itself, all in
			quotes (e.g., 'myfile.xlsx').

		calculations: pandas DataFrame or list of pandas DataFrames
			A single variable or list of variables containing calculated outputs from any of the
			core ExcelFile functions

		sheet_name: None, string, or list
			OPTIONAL. Default value is None. Allows user to set the name of the sheet or sheets
			written to the Excel file.

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
					raise core.InputError("calculations and sheet_name must have the same length")

				for i in range(len(calculations)):
					if isinstance(sheet_name[i], str):
						calculations[i].to_excel(writer, sheet_name[i])
					else:
						raise core.InputError("if sheet_name is passed, it must be list of strings")

		return print("Saved " + str(filename))

	def get_silicate_comp(self, sample, norm='none'):
		"""
		Returns wt% oxide composition of a single silicate sample from a user-imported excel file
		as a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied
			to the sample.

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
			raise core.InputError('norm must be either none or standard.')

		data = preprocess_dataset(self.silicate_data, how='oxides')
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])

		sample_comp = sample_class.get_oxides(sample_dict)

		if norm == 'standard':
			return sample_class.normalize(sample_comp)
		if norm == 'none':
			return sample_comp

	def get_metal_comp(self, sample, metal_1=True, metal_2=False, norm='none'):
		"""
		Returns wt% element composition of a single metal sample from a user-imported excel file as
		a dictionary

		Parameters
		----------
		sample: string
			Name of the desired sample

		metal_1: bool
			OPTIONAL. Detault is True. If True, returns composition of first metal "metal 1"
			(self.metal_data).
		
		metal_2: bool
			OPTIONAL. Detault is False. If True, returns composition of second metal "metal 2"
			(self.metal_2_data). Note if both metal_1 and metal_2 are set to True, this function
			returns to objects.

		norm_style: string
			OPTIONAL. Default value is 'standard'. This specifies the style of normalization applied
			to the sample.

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
			Composition of the metal phase in wt%. If both metal_1 and metal_2 are set to True, two
			objects are returned.
		"""
		if norm == 'none' or norm == 'standard':
			pass
		else:
			raise core.InputError('norm must be either none or standard.')

		if metal_1 == True:
			data = preprocess_dataset(self.metal_data, how='elements')
			my_sample = pd.DataFrame(data.loc[sample])
			sample_dict = (my_sample.to_dict()[sample])

			sample_comp = sample_class.get_elements(sample_dict)

			if norm == 'standard':
				metal_1_return = sample_class.normalize(sample_comp)
			if norm == 'none':
				metal_1_return = sample_comp

			if metal_2 == False:
				return metal_1_return
		
		if metal_2 == True:
			data2 = preprocess_dataset(self.metal_2_data, how='elements')
			my_sample = pd.DataFrame(data2.loc[sample])
			sample_dict = (my_sample.to_dict()[sample])

			sample_comp = sample_class.get_elements(sample_dict)

			if norm == 'standard':
				metal_2_return = sample_class.normalize(sample_comp)
			if norm == 'none':
				metal_2_return = sample_comp

			if metal_1 == False:
				return metal_2_return

		if metal_1 == True and metal_2 == True:
			return metal_1_return, metal_2_return

	def get_experimental_metadata(self, sample):
		"""
		Returns experimental metadata of a single metal sample from a user-imported excel file as a
		dictionary

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
			raise core.InputError('norm must be either none or standard.')

		data = self.experimental_data
		my_sample = pd.DataFrame(data.loc[sample])
		sample_dict = (my_sample.to_dict()[sample])

		return sample_dict

	def process(self, pressure, temperature, species=core.corgne_species,
			 interactions=core.standard_interactions, filename=None, print_status=True):
		"""
		Process the excel file and save a standardized set of data to an excel file that can be used
		for
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
				raise core.InputError("filename must be type str")

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
			calc = sample_class.wtpercentOxides_to_molOxides(wtpercentOxides)
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
				calc = sample_class.wtpercentElements_to_molElements(wtpercentElements)
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
				calc2 = sample_class.wtpercentElements_to_molElements(wtpercentElements2)
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
				calc = calculate.calc_NBO_T(sample=bulk_comp, method=method, metadata=metadata)
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

	def calc_gamma_Fe_metal(self, temperature, interactions=core.standard_interactions, print_warnings=False, print_status=True):
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
			raise core.InputError("temp must be type str or float or int")

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
			calc = calculate.calc_gamma_Fe_metal(metal_comp=bulk_comp, temperature=temperature, interactions=interactions, print_warnings=print_warnings)
			gamma_Fe_vals.append(calc)
			temp_vals.append(temperature)
		
		user_data["Temperature_C_Modeled"] = temp_vals
		user_data["gamma_Fe_metal"] = gamma_Fe_vals

		return_data = user_data.filter(['Temperature_C_Modeled', 'gamma_Fe_metal'], axis=1)

		return return_data

	def calc_gamma_solute_metal(self, temperature, species=None, interactions=core.standard_interactions, gammaFe="calculate", print_warnings=False, print_status=False):
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
			raise core.InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise core.InputError("temp must be type str or float or int")

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
				calc = calculate.calc_gamma_solute_metal(metal_comp=bulk_comp, species=spec, temperature=temperature, interactions=interactions, gammaFe=gammaFe_single, print_warnings=print_warnings)
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
			calc = calculate.calc_gamma_FeO_silicate(silicate_comp=bulk_comp)
			gamma_FeO_vals.append(calc)

		return_data = pd.DataFrame(index=user_data.index)
		return_data["gamma_FeO_silicate"] = gamma_FeO_vals

		return return_data

	def metal_activity_from_composition(self, species, temperature, interactions=core.standard_interactions):
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
			raise core.InputError("temp must be type str or float or int")

		if isinstance(species, str):
			species = [species] #if string, turn into list of length 1
		if isinstance(species, list):
			pass
		else:
			raise core.InputError("Species must be input as string or list of strings. You have input " + str(type(species)))

		return_data = pd.DataFrame(index=user_data.index)
		for s in species:
			activity_vals = []
			for index, row in user_data.iterrows():
				if file_has_temp == True:
					temperature = self.experimental_data.loc[index][temp_name]
				bulk_comp = self.get_metal_comp(index)
				try:
					calc = calculate.metal_activity_from_composition(metal_comp=bulk_comp, species=s, temperature=temperature, interactions=interactions)
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

	def calc_dIW(self, temperature, interactions=core.standard_interactions, gammaFe="calculate", gammaFeO="calculate", print_warnings=False, print_status=True):
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
			raise core.InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")

		if isinstance(gammaFeO, dict):
			has_gammaFeO = True
		elif gammaFeO == "calculate":
			gammaFeO_single = gammaFeO
			has_gammaFeO = False
		else:
			raise core.InputError("If passed, gammaFe argument must be 'Calculate' or type dict. Not sure what you want here? Don't pass this argument.")


		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise core.InputError("temp must be type str or float or int")

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
			calc = calculate.calc_dIW(silicate_comp=silicate_comp, metal_comp=metal_comp, temperature=temperature, interactions=interactions, gammaFe=gammaFe_single, gammaFeO=gammaFeO_single, print_warnings=print_warnings)
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
			raise core.InputError("pressure must be type str or float or int")

		if isinstance(temperature, str):
			file_has_temp = True
			temp_name = temperature
		elif isinstance(temperature, float) or isinstance(temperature, int):
			file_has_temp = False
		else:
			raise core.InputError("temperature must be type str or float or int")

		if isinstance(dIW, str):
			if dIW == "calculate":
				file_has_dIW = False
				dIW_to_be_calcd = True
				if 'interactions' not in kwargs:
					interactions=core.standard_interactions
				if 'gammaFe' not in kwargs:
					gammaFe = "calculate"
				if 'gammaFeO' not in kwargs:
					gammaFeO = "calculate"
				if 'print_warnings' not in kwargs:
					print_warnings=False
				if print_status == True:
					print("Calculating dIW values...")
				dIW_calcd_vals = self.calc_dIW(temperature=temperature, interactions=interactions,
								               gammaFe=gammaFe, gammaFeO=gammaFeO,
											   print_warnings=print_warnings,
											   print_status=print_status)
			else:
				file_has_dIW = True
				dIW_name = dIW
		elif isinstance(dIW, float) or isinstance(dIW, int):
			file_has_dIW = False
		elif isinstance(dIW, pd.core.series.Series):
			dIW_is_series = True
			dIW_series = dIW
		else:
			raise core.InputError("dIW must be type str or float or int or pandas Series")

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

			calc = calculate.calc_logfO2_from_IW(pressure=pressure, temperature=temperature, dIW=dIW)
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
			calc = calculate.calc_optical_basicity(silicate_comp)
			optical_bas_list.append(calc)

		return_data = pd.DataFrame(index=silicate_data.index)
		return_data["Lambda"] = optical_bas_list

		return return_data

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