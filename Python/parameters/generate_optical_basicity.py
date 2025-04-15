import pandas as pd
import numpy as np
optical_basicity = pd.read_excel('optical_basicity.xlsx',)
optical_basicity.fillna("No Data")

import math
f= open("opticalbasicity.py","w+")
cols = ['Phase', 'Coordination', 'Lambda']

f.write("\n")
f.write("optical_basicity = pd.DataFrame({ \n")
for col in cols:
	iterno = 1
	f.write("'" + (str(col)+"': ["))
	for index, row in optical_basicity.iterrows():
		try:
			value = float(row[col])
			if str(value) == 'nan':
				f.write("'No Data'")
			else:
				f.write(str(value))
		except:
			f.write("'" + str(row[col]) + "'")
		if iterno < len(optical_basicity.index):
			f.write(",")
		iterno += 1
	f.write("], \n")
f.write(" }).set_index('Phase') \n")