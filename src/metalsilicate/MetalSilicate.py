
from mendeleev.fetch import fetch_table


# import statsmodels.api as sm
# import statsmodels.formula.api as smf
# from patsy import dmatrices


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
		NOTE: do we want to calculate activities for a bunch of elements in the metal? Then we could
		use this instead of the straight up mole fractions or wt% conentrations in the regression.
		This might be best if we want to compare data with a large range of metal compositions (e.g.
		in terms of S, C, Si). So we can calculate the activities of at least those components, and
		use X for the others?

	1. Lookup EPSILON_i_j at a reference T
	2. Calculate EPSIOLON_i_j at sample T
	3. Lookup GAMMA_i_0 at reference T
	4. Calculate GAMMA_i_0 at sample T
	5. Calculate activity coefficient GAMMA_i from ref value GAMMA_i_0 and interaction parameter
	EPSILON_i_j at sample T
	6. Calculate the activity of the solute as mole fraction times gamma.

Calculate the fO2:
	1. Calcualte Fe activity in the metal
	2. Calculate FeO activity in the silicate
	3. Calculate deltaIW fO2 from activity from mole fraction and GAMMA of FeO in silicate and Fe in
	metal (see Corgne et al., 2007)
"""



# def calc_gamma_FeO_silicate_test(silicate_comp):
# 	"""
# 	Returns a value for gammaFeO in the silicate. Parameterization is based on Holdzheid, where
#	gammaFeO is taken as a constant value from 1.7-3, dependent only upon MgO content. 

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



