"""
metalsilicate

Python library for calculating thermodynamic parameters of geochemical systems with both silicate
melt and Fe-rich metal in equilibrium.

Originally developed by Iacovino et al. to process a partitioning database and data from completed
gas levitation laser melting experiments to develop a thermodynamic model for Si partitioning
between silicate melt and Fe-rich metal in support of assessing the cataclysmic reduction hypothesis
for Mercury’s origin
"""

__version__ = "0.1.0"
__author__ = "Kayla Iacovino"

# ----------------- IMPORTS ----------------- #
import warnings as w
import pandas as pd

import metalsilicate.core
import metalsilicate.batch_calculate
import metalsilicate.batchfile
import metalsilicate.calculate
import metalsilicate.sample_class

# -------------- BATCH PROCESSING ------------ #
class ExcelFile(metalsilicate.batchfile.ExcelFile):
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