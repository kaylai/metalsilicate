import MetalSilicate as ms
import pandas as pd
from datetime import datetime

pd.set_option('display.max_rows', None, 'display.max_columns', None)

myfile = ms.ExcelFile('master.xlsx', sheet_name='Data')

#Settings for Saved File Name
todays_date = datetime.today().strftime('%Y-%m-%d')
filename = 'MS_out_' + todays_date + '.xlsx'

#Process data
processed = myfile.process(pressure="P (GPa)", temperature="T_C", filename=filename)

