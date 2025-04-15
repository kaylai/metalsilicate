# Python 3.5
# Script for MetalSilicate python library
# Copyright Kayla Iacovino
# This script contains pandas dataframes with activity coefficient data digitized from the Steelmaking Data Source Book

import pandas as pd

# These hard-coded pandas DataFrames are generated from an Excel Spreadsheet and hardcoded for easier access.
"""
The code to generate these:
import pandas as pd
import numpy as np
activity_coeffs_df = pd.read_excel('ActivityCoeffs_and_InteractionParams.xlsx', sheet_name='ActivityCoefficients')
activity_coeffs_df.fillna("No Data")

import math
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


##NOTE! Delete the final comma at the end of the final column.

##Reminder: to look up a value, use syntax: activity_coeffs.loc['Ag']['gamma_naught_i']
"""


activity_coeffs = pd.DataFrame({
'i': ['Ag','Al','Al_delta','B','C','C_austenite','C_ferrite','Ca','Ce','Co','Cr','Cr_solid','Cu','H','H_delta','H_austenite','H_ferrite','La','Mn','Mo','Mo_solid','N','N_delta','N_austenite','N_ferrite','Nb','Nb','Ni','O','P','Pb','S','Si','Sn','Ta','Ti','Ti_solid','U','V','V_solid','W','W_solid','Zr','Zr_solid'], 
'Element_State': ['l','l','l','s','gr','gr','gr','l','l','l','l','s','l','g','g','g','g','l','l','l','s','g','g','g','g','l','s','l','g','g','l','g','l','l','l','l','s','l','l','s','l','s','l','s'], 
'Fe_State': ['liquid iron','liquid iron','delta iron','liquid iron','liquid iron','austenite','ferrite','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','delta iron','austenite','ferrite','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','delta iron','austenite','ferrite','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron','liquid iron'], 
'gamma_naught_i': [200.0,0.049,0.027,0.022,0.538,7.49,313.0,2270.0,0.322,0.55,1.0,1.14,8.58,'–','–','–','–',9.3,1.44,1.0,2.2,'–','–','–','–',0.2,1.4,0.66,'–','–',837.0,'–',0.0013,2.58,0.04,0.004,0.009,0.027,0.08,0.1,1.0,7.6,0.037,0.043], 
'Temp_K': [1873.0,1873.0,1673.0,1873.0,1873.0,1273.0,1073.0,1880.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1773.0,1273.0,1073.0,1873.0,1843.0,1823.0,1823.0,1873.0,1773.0,1273.0,1073.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0,1873.0], 
'DeltaG_J_per_g-atom': ['82400-43.76T','-17100-19.4T','-17900-19.3T','-65300-21.55T','17230-39.87T','72180+52.05T*logT-227.03T','83720-55.76T','No Data','-16700-46.44T','No Data','-37.70T','19200-46.86T','33500-39.37T','36460+30.46T','28650+45.35T','22630+45.35T','28650+45.35T','125,900-94.6T','No Data','No Data','No Data','9916+20.17T','29090+19.91T','-8620+37.42T','29090+19.1T','No Data','23000-52.3T','-20900-31.1T','-117100-3.39T','-157700+5.4T','No Data','-125100+18.5T','-131500-17.24T','No Data','No Data','No Data','No Data','-56100-50.3T','-42300-36.0T','-20700-45.6T','-48.1T','No Data','-51000-42.38T','-34700-50.00T'], 
'TempRange_K': ['No Data','m.p.–1873','1673–m.p.','No Data','1773–2073','1773–1673','923–1073','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data','No Data'], 
'Ref': [2.0,6.0,6.0,2.0,90.0,76.0,100.0,17.0,102.0,203.0,2.0,2.0,2.0,22.0,3434.0,434.0,434.0,282.0,288.0,436.0,436.0,437.0,438.0,439.0,438.0,440.0,2.0,2.0,360.0,441.0,41.0,442.0,2.0,513.0,440.0,440.0,385.0,2.0,2.0,2.0,2.0,436.0,2.0,2.0], 
'Year': [1974.0,1977.0,1977.0,1974.0,1962.0,1970.0,1967.0,1964.0,1977.0,1978.0,1974.0,1974.0,1974.0,1963.0,1950.0,1950.0,1950.0,1977.0,1982.0,1984.0,1984.0,1982.0,1973.0,1955.0,1973.0,1973.0,1974.0,1974.0,1959.0,1980.0,1971.0,1981.0,1974.0,1981.0,1973.0,1973.0,1984.0,1974.0,1974.0,1974.0,1974.0,1984.0,1974.0,1974.0], 
'Note': ['From solubility data, regular solution assumed','From solubility data, regular solution assumed','No Data','No Data','No Data','No Data','No Data','No Data','No Data','log(gammaCO) = -0.257N^2_Fe','Ideal solution assumed','No Data','No Data','1/2H2','1/2H2','1/2H2','1/2H2','No Data','No Data','No Data','±0.17','1/2N2','1/2N2','1/2N2','1/2N2','No Data','No Data','No Data','1/2O2; Data were combined with deltaG_0 = -251877+58.325T for reaction H2(g) + 1/2O2(g) = H2O(g) [14]','1/2P2','No Data','1/2S2','No Data','Mass analysis','Calculated value','No Data','Calculated value','Regular solution assumed','No Data','No Data','No Data','±1.1','Gamma_0 Zr assumed equal to Gamma_0 Ti and regular solution assumed','No Data'] 
 }).set_index('i') 


