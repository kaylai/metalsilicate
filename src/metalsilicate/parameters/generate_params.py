import pandas as pd
import numpy as np
interaction_parameter_df = pd.read_excel('ActivityCoeffs_and_InteractionParams.xlsx', sheet_name='InteractionParameters_Liquid')
interaction_parameter_df.fillna("No Data")

activity_coeffs_df = pd.read_excel('ActivityCoeffs_and_InteractionParams.xlsx', sheet_name='ActivityCoefficients')
activity_coeffs_df.fillna("No Data")

import math
f= open("hard_coded_calibrations.py","w+")
cols = ['i,j,k', 'ei(j,k)', 'Temp_K', 'ConcRange_MassPercent', 'TempDependency', 'TempRange_K', 'Reference', 'Note']

f.write("\n")
f.write("interaction_params = pd.DataFrame({ \n")
for col in cols:
	iterno = 1
	f.write("'" + (str(col)+"': ["))
	for index, row in interaction_parameter_df.iterrows():
		try:
			value = float(row[col])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[col]) + "'")
		if iterno < len(interaction_parameter_df.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }).set_index('i,j,k') \n")

f = open("hard_coded_activity_coefficients.py","w+")
cols = ['i', 'Element_State', 'Fe_State', 'gamma_naught_i', 'Temp_K', 'DeltaG_J_per_g-atom', 'TempRange_K', 'Ref', 'Year', 'Note']

f.write ("\n")
f.write("activity_coeffs = pd.DataFrame({\n")
for col in cols:
	iterno = 1
	f.write("'" + (str(col)+"': ["))
	for index, row in activity_coeffs_df.iterrows():
		try:
			value = float(row[col])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[col]) + "'")
		if iterno < len(activity_coeffs_df.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }).set_index('i') \n")