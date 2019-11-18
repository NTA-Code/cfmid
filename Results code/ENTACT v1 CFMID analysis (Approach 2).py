# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:44:37 2019

@author: AChao
"""


############################################################################
# This script is to process: 
# 1. ENTACT v1 MS2 PCDL information
# 2. ENTACT v1 MS2 CFMID search results (using summed score data)
############################################################################

############################################################################
# Debugging key
# df1 = The spiked list (manual review file) for each mixture's spiked compounds. Contains manual review information (Pass/fail, found in PCDL matching etc.)
# df2 = CFMID search results (from CFMID batch search)
# df3 = Agilent's PCDL-exported compound list
############################################################################

import os
import time
import pandas as pd

start_time = time.time()

os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/Results code")

########################################################################
# Read in Spiked list of chemicals (which also has PCDL information)
########################################################################

xls = pd.ExcelFile('Python Combined ENTACT v1 manual review results(re-PCDL searched with less stringent parameters).xlsx')

sheet_list = ['499', '500', '501', '502', '503', '504', '505', '506', '507', '508']
#sheet_list = ['499']

counter = 1

for i in sheet_list:
    print('Reading spiked list sheet', i)
    if counter == 1:
        df1 = pd.read_excel(xls, sheetname=i)
        df1['mixture'] = i
    else:
        df_temp = pd.read_excel(xls, sheetname=i) ## Also tried sheet_name = 1
        df_temp['mixture'] = i
        df1 = pd.concat([df1, df_temp])
    counter += 1

df1.rename(columns={'MS_Ready_DTXCID':'OLDCID'}, inplace=True) # Rename DTXCID column to "DTXCID"
df1.drop('in PCDL?', axis=1, inplace=True)
df1['cfmid_match'] = 'Y'

df1['MSMS present pos?'].fillna('N', inplace=True)
df1['MSMS present neg?'].fillna('N', inplace=True)

# This below section combines whether MSMS was seen in either pos/neg mode as a general was it seen in any mode
def find_MSMS (row):
    if row['MSMS present pos?'][0] == 'Y' :
        return 'Y'
    if row['MSMS present neg?'][0] == 'Y' :
        return 'Y'
    return 'N'
    
df1['MSMS present?'] = df1.apply (lambda row: find_MSMS (row), axis=1)


########################################################################
# Read in the PCDL file (containing list of compounds within PCDLs)
# This file is already filtered for only those compounds with MSMS spectra > 1, and de-duplicated

# Read in SID -> CID file (DSSToxMS-Ready 09102018)
# Convert SID's to CID's
########################################################################

df3 = pd.read_csv('PCDL_DTXSID_with_MS.csv')
print('PCDL file read')
df3.rename(columns={'DTXSID_combined':'DTXSID'}, inplace=True) # Rename DTXSID column
df3['in PCDL'] = 'Y'

# This needs to be uncommented and run once. Once the DSSTox is loaded into memory it can be commented out

df_conv = pd.read_csv('DSSToxMS-Ready 09102018_SIDCID_only.csv')


# Add DTXCID's from 09102018 database to PCDL file
df3 = pd.merge(df3, df_conv[['DTXSID', 'DTXCID']], how='left', on = 'DTXSID')

# Add DTXCID's from 09102018 database to spiked list file
df1 = pd.merge(df1, df_conv[['DTXSID', 'DTXCID']], how='left', on = 'DTXSID')


########################################################################
# Merge PCDL list with spiked list to see which compounds are in PCDL's
########################################################################

df1 = pd.merge(df1, df3[['DTXCID', 'in PCDL']], how='left', on = 'DTXCID')

# Filter out df1 spiked list by whether it was a Pass
df1 = df1.loc[df1['Decision'] == 'Pass']


########################################################################
# Read in CFMID search results - Score summed across CE CFMID results files
########################################################################

############2/6/2019 Adjust read-in to not have CE column... keep old code, because it might be useful for CE comparisons

os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/Results code/CFMID search results")

counter = 1

# Iterate through all the mixtures selected above
# Also, iterate through each energy level for each mixture

for i in sheet_list:
    print('Reading cfmid file', i)
    name_1 = i + '_pos_C_01_CFMID_OneScore_AllHits.xlsx'
    name_2 = i + '_pos_C_02_CFMID_OneScore_AllHits.xlsx'
    name_3 = i + '_pos_C_03_CFMID_OneScore_AllHits.xlsx'
    name_4 = i + '_neg_C_01_CFMID_OneScore_AllHits.xlsx'
    name_5 = i + '_neg_C_02_CFMID_OneScore_AllHits.xlsx'
    name_6 = i + '_neg_C_03_CFMID_OneScore_AllHits.xlsx'
    if counter == 1:
        # Read in Positive
        df2 = pd.read_excel(name_1) # Read in first file, assign it CE 10
        df2['CE'] = 10
        
        df_temp = pd.read_excel(name_2) # Read in second file, assign it CE 20
        df_temp['CE'] = 20     
        df2 = pd.concat([df2, df_temp])
    
        df_temp = pd.read_excel(name_3) # Read in third file, assign it CE 40
        df_temp['CE'] = 40
  
        df2 = pd.concat([df2, df_temp])
        df2['mixture'] = i
        df2['cfmid_mode'] = 'Esi+'

        # Read in Negative
        df_temp = pd.read_excel(name_4) # Read in first neg file, assign it CE 10
        df_temp['CE'] = 10
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_5) # Read in second neg file, assign it CE 20
        df_temp['CE'] = 20
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'         
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_6) # Read in third neg file, assign it CE 40
        df_temp['CE'] = 40
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'
        df2 = pd.concat([df2, df_temp])        
        
    else:
        df_temp = pd.read_excel(name_1) # Read in first file, assign it CE 10
        df_temp['CE'] = 10
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi+'         
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_2) # Read in second file, assign it CE 20
        df_temp['CE'] = 20
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi+'         
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_3) # Read in third file, assign it CE 40
        df_temp['CE'] = 40
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi+'        
        df2 = pd.concat([df2, df_temp])
        
        # Read in Negative
        df_temp = pd.read_excel(name_4) # Read in first neg file, assign it CE 10
        df_temp['CE'] = 10
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'        
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_5) # Read in second neg file, assign it CE 20
        df_temp['CE'] = 20
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'       
        df2 = pd.concat([df2, df_temp])
        
        df_temp = pd.read_excel(name_6) # Read in third neg file, assign it CE 40
        df_temp['CE'] = 40
        df_temp['mixture'] = i
        df_temp['cfmid_mode'] = 'Esi-'         
        df2 = pd.concat([df2, df_temp])              
        
    counter += 1



# Match the spiked list compounds to the CFMID results by DTXCID. Anything in the spiked list will get a "y" in cfmid_match column
df2_new = pd.merge(df2, df1[['DTXCID', 'cfmid_match', 'mixture']], how='left', on = ['DTXCID', 'mixture'])
df2_new = df2_new.drop_duplicates() # Duplicates can occur if the compound appears in multiple rows in the spiked list

# The following lines rank the compounds (grouped by a single MGF mass) for each energy level
df2_new['rank_bymass'] = df2_new.groupby(['MASS_in_MGF', 'CE'])['SCORE'].rank(ascending=False)

# The following lines generate percentile values (grouped by a single MGF mass) for each energy level
df2_new['percentile_by_mass'] = df2_new.groupby(['MASS_in_MGF', 'CE'])['SCORE'].rank(pct=True)
df2_new['percentile_by_mass'] = df2_new['percentile_by_mass'] * 100
# by formula
df2_new['percentile_by_formula'] = df2_new.groupby(['MASS_in_MGF', 'CE', 'FORMULA'])['SCORE'].rank(pct=True)
df2_new['percentile_by_formula'] = df2_new['percentile_by_formula'] * 100


# Generate quotient values for each energy level
df2_new['quot_by_mass_max'] = df2_new.groupby(['MASS_in_MGF', 'mixture', 'CE'])['SCORE'].transform('max')
df2_new['quot_by_mass'] = df2_new['SCORE']/df2_new['quot_by_mass_max']
df2_new['quot_by_form_max'] = df2_new.groupby(['MASS_in_MGF', 'FORMULA', 'mixture', 'CE'])['SCORE'].transform('max')
df2_new['quot_by_form'] = df2_new['SCORE']/df2_new['quot_by_form_max']



# In the case that there is no energy score (i.e. no spectrum match), make the rank -1 as a flag for this condition
df2_new['rank_bymass'].loc[df2_new['SCORE'] == 0] = -1 # These are mass-matched rankings
df2_new['RANK'].loc[df2_new['SCORE'] == 0] = -1 # These are formula-matched rankings

df2_new.rename(columns={'MATCHES':'formula_matches'}, inplace=True) # Rename matches to formula matches
df2_new['total_matches'] = df2_new.groupby(['MASS_in_MGF', 'CE'])['DTXCID'].transform('count') # Update the matches column for all matches to the mass (Hussein's code does matches to the formula)

df2_match = df2_new.loc[df2_new['cfmid_match'] == 'Y']
df2_match_pos = df2_match.loc[df2_match['cfmid_mode'] == 'Esi+']
df2_match_neg = df2_match.loc[df2_match['cfmid_mode'] == 'Esi-']


########################################################################
# Merge cfmid results information into the spiked list dataframe
########################################################################
df1.drop('cfmid_match', axis=1, inplace=True)
df1 = df1.rename(columns = {'Re-search match?':'pcdl_match'})
df1_cfmid = pd.merge(df1[['DTXSID', 'Preferred_Name', 'Python MSMS present', 'DTXCID', 'mixture', 'MS_Ready_Mol_Formula', 
                          'MS_Ready_Monoisotopic_Mass', 'Ionization_Mode', 'Mass', 'Est_mz', 'Retention_Time', 'MSMS present pos?', 
                          'MSMS present neg?', 'MSMS present?', 'pcdl_match', 'in PCDL', 'Comment', 'Decision', 'Star_Rating']], 
    df2_match[['DTXCID', 'mixture', 'cfmid_mode', 'CE', 'cfmid_match', 'rank_bymass', 'total_matches', 'RANK', 'formula_matches', 'SCORE',
               'MASS_in_MGF', 'percentile_by_mass', 'percentile_by_formula', 'quot_by_mass', 'quot_by_form']], 
               how='left', on = ['DTXCID', 'mixture'])

########################################################################
# Break up CFMID results into specific energy levels to generate results from
########################################################################

df1_cfmid_CE10 = df1_cfmid.loc[df1_cfmid['CE'] == 10]
df1_cfmid_CE20 = df1_cfmid.loc[df1_cfmid['CE'] == 20]
df1_cfmid_CE40 = df1_cfmid.loc[df1_cfmid['CE'] == 40]
df1_cfmid_NA = df1_cfmid.loc[df1_cfmid['CE'].isnull()]


# De-duplicate CFMID results via CID and select for only passes
df1_cfmid_CE10 = df1_cfmid_CE10[df1_cfmid_CE10['Decision'] == 'Pass']
df1_cfmid_CE20 = df1_cfmid_CE20[df1_cfmid_CE20['Decision'] == 'Pass']
df1_cfmid_CE40 = df1_cfmid_CE40[df1_cfmid_CE40['Decision'] == 'Pass']

df1_cfmid_CE10 = df1_cfmid_CE10.drop_duplicates(subset=['DTXCID'])
df1_cfmid_CE20 = df1_cfmid_CE20.drop_duplicates(subset=['DTXCID'])
df1_cfmid_CE40 = df1_cfmid_CE40.drop_duplicates(subset=['DTXCID'])


'''
os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/Results code")
df1_cfmid_CE10.to_csv('Approach 2 exp CE10 ENTACT compound cfmid results.csv', sep=',', encoding='utf-8', index=False)
df1_cfmid_CE20.to_csv('Approach 2 exp CE20 ENTACT compound cfmid results.csv', sep=',', encoding='utf-8', index=False)
df1_cfmid_CE40.to_csv('Approach 2 exp CE40 ENTACT compound cfmid results.csv', sep=',', encoding='utf-8', index=False)
'''
    
print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))

