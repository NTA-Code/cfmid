# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:07:47 2019

@author: AChao
"""

############################################################################
# This script is to process: 
# 1. ENTACT v1 MS2 PCDL information
# 2. ENTACT v1 MS2 CFMID search results (using one score data)
############################################################################

import os
import time
import pandas as pd
import numpy as np

time_list = []
start_time = time.time()
time_list.append(time.time())


os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/Results code")


########################################################################
# Read in Spiked list of chemicals (which also has PCDL information)
########################################################################

xls = pd.ExcelFile('Python Combined ENTACT v1 manual review results(re-PCDL searched with less stringent parameters).xlsx')

sheet_list = ['499', '500', '501', '502', '503', '504', '505', '506', '507', '508']
#sheet_list = ['499']


counter = 1

for i in sheet_list:
    #print(i)
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

# Rename spiked list formula for matching later
df1.rename(columns={'Mol_Formula':'spiked_formula'}, inplace=True) 

########################################################################
# Read in the PCDL file (containing list of compounds within PCDLs)
# This file is already filtered for only those compounds with MSMS spectra > 1, and de-duplicated

# Read in SID -> CID file (DSSToxMS-Ready 09102018)
# Convert SID's to CID's
########################################################################

time_list.append(time.time())
print ("\n\nCurrent total runtime: --- %s seconds ---" % round((time_list[len(time_list)-1] - time_list[0]), 1))
print('Reading in PCDLs...')


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

########################################################################
# Read in CFMID search results - This is the file from the other script:
# ENTACT v1 CFMID blinded analyais of data.py
# OneScores across three exp. CE's summed
# So essentially 9 scores summed into one at the very end
########################################################################
time_list.append(time.time())
print ("\n\nTime since last point: --- %s seconds ---" % round((time_list[len(time_list)-1] - time_list[len(time_list)-2]), 1))
print ("Current total runtime: --- %s seconds ---" % round((time_list[len(time_list)-1] - time_list[0]), 1))

print('reading in CFMID results...\n')
############2/6/2019 Adjust read-in to not have CE column... keep old code, because it might be useful for CE comparisons

# Below csv generated from "ENTACT v1 CFMID blinded analysis of data.py"
df2 = pd.read_csv('ENTACT_CFMID_all_mixtures_merged_results_rounded.csv')

# Need to convert df1 and df2 mixture data types to match (otherwise they can't be merged on mixture)
df1['mixture'] = df1['mixture'].astype(np.int64)

# Pull out first character of pcdl match column into new column
df1 = df1.rename(columns = {'Re-search match?':'pcdl_match'})
df1.loc[df1['pcdl_match'].notnull(), 'pcdl_match_1'] = df1['pcdl_match'].str[0]

# Get just PCDL matches
df1_pcdl = df1.loc[df1['pcdl_match_1'] == 'Y']
# de-duplicate PCDL matches by CID
df1_pcdl = df1_pcdl.drop_duplicates(subset=['DTXCID'])


# Filter out df1 spiked list by whether it was a Pass
df1 = df1.loc[df1['Decision'] == 'Pass']
df1_pass = df1.loc[df1['Decision'] == 'Pass']
df1_pass = df1_pass.drop_duplicates(subset=['DTXSID'])

# Pull out pcdl matches only from the passes
df1_pcdl_pass = df1_pass.loc[df1_pass['pcdl_match_1'] == 'Y']

# Pull out those with python MS2 verified presence from the passes
df1_pass_python_ms =  df1_pass.loc[df1_pass['Python MSMS present'] == 'Y']
df1_pass_python_ms = df1_pass_python_ms.drop_duplicates(subset=['DTXSID'])

# Pull out those with python MS2 verified presence AND in PCDL from the passes
df1_pass_python_ms_in_pcdl = df1_pass_python_ms.loc[df1_pass_python_ms['in PCDL'] == 'Y']



# Store the spiked compound DTXCID in a new column
df1['spiked_DTXCID'] = df1['DTXCID']


# Match the spiked list compounds to the CFMID results by DTXCID. Anything in the spiked list will get a "y" in cfmid_match column
df2_new = pd.merge(df2, df1[['DTXCID', 'spiked_formula', 'spiked_DTXCID', 'cfmid_match', 'mixture']], how='left', on = ['DTXCID', 'mixture'])
df2_new = df2_new.drop_duplicates() # Duplicates can occur if the compound appears in multiple rows in the spiked list


# The following lines rank the compounds (grouped by a single MGF mass) for each energy level
df2_new['rank_bymass'] = df2_new.groupby(['MASS_in_MGF', 'mixture'])['SCORE'].rank(ascending=False)
df2_new['RANK'] = df2_new.groupby(['FORMULA', 'mixture', 'MASS_in_MGF'])['SCORE'].rank(ascending=False)

# Get percentiles for the compounds (grouped by a single MGF mass)
df2_new['percentile_by_mass'] = df2_new.groupby(['MASS_in_MGF', 'mixture'])['SCORE'].rank(pct=True)
df2_new['percentile_by_mass'] = df2_new['percentile_by_mass'] * 100
#df2_new['percentile_by_formula'] = df2_new.groupby(['FORMULA', 'mixture'])['SCORE'].rank(pct=True)
df2_new['percentile_by_formula'] = df2_new.groupby(['MASS_in_MGF', 'FORMULA', 'mixture'])['SCORE'].rank(pct=True)
df2_new['percentile_by_formula'] = df2_new['percentile_by_formula'] * 100

# Generate quotient values for the compounds (for both mass and formula matching)
df2_new['cfmid_score_mass_max'] = df2_new.groupby(['MASS_in_MGF', 'mixture'])['SCORE'].transform('max')
df2_new['cfmid_score_mass_quot'] = df2_new['SCORE']/df2_new['cfmid_score_mass_max']
df2_new['cfmid_score_form_max'] = df2_new.groupby(['MASS_in_MGF', 'FORMULA', 'mixture'])['SCORE'].transform('max')
df2_new['cfmid_score_form_quot'] = df2_new['SCORE']/df2_new['cfmid_score_form_max']


# In the case that there is no energy score (i.e. no spectrum match), make the rank 10000000 as a flag for this condition
df2_new['rank_bymass'].loc[df2_new['SCORE'] == 0] = -1 # These are mass-matched rankings
df2_new['RANK'].loc[df2_new['SCORE'] == 0] = -1 # These are formula-matched rankings

df2_new.rename(columns={'MATCHES':'formula_matches'}, inplace=True) # Rename matches to formula matches
df2_new['total_matches'] = df2_new.groupby(['MASS_in_MGF', 'mixture'])['DTXCID'].transform('count') # Update the matches column for all matches to the mass (Hussein's code does matches to the formula)



df2_match = df2_new.loc[df2_new['cfmid_match'] == 'Y']
df2_match_pos = df2_match.loc[df2_match['cfmid_mode'] == 'Esi+']
df2_match_neg = df2_match.loc[df2_match['cfmid_mode'] == 'Esi-']

########################################################################
# Merge cfmid results information into the spiked list dataframe
########################################################################
df1.drop('cfmid_match', axis=1, inplace=True)
df1_cfmid = pd.merge(df1[['DTXSID', 'Preferred_Name', 'Python MSMS present', 'DTXCID', 'mixture', 'MS_Ready_Mol_Formula', 'MS_Ready_Monoisotopic_Mass', 'Ionization_Mode', 'Mass', 'Est_mz', 'Retention_Time', 'MSMS present pos?', 'MSMS present neg?', 'pcdl_match', 'in PCDL', 'Comment', 'Decision', 'Star_Rating']], df2_match[['DTXCID', 'mixture', 'cfmid_mode', 'cfmid_match', 'percentile_by_mass', 'percentile_by_formula', 'cfmid_score_mass_quot', 'cfmid_score_form_quot', 'rank_bymass', 'total_matches', 'RANK', 'formula_matches', 'SCORE', 'MASS_in_MGF']], how='left', on = ['DTXCID', 'mixture'])

########################################################################
# Break up CFMID results into specific energy levels to generate results from
########################################################################
df1_cfmid_pos = df1_cfmid.loc[df1_cfmid['cfmid_mode'] == 'Esi+']
df1_cfmid_neg = df1_cfmid.loc[df1_cfmid['cfmid_mode'] == 'Esi-']

# De-duplicate CFMID results via both CID and MASS_in_MGF (one row per mgf precursor/compound combination)
df1_cfmid_pos = df1_cfmid_pos.drop_duplicates(subset=['DTXCID', 'MASS_in_MGF', 'mixture'])
df1_cfmid_neg = df1_cfmid_neg.drop_duplicates(subset=['DTXCID', 'MASS_in_MGF', 'mixture'])

df1_cfmid_both = pd.concat([df1_cfmid_pos, df1_cfmid_neg]) # Combine positive and negative results into one dataframe
df1_cfmid_both = df1_cfmid_both.drop_duplicates(subset=['DTXCID']) #Drop duplicates of DTXCID

# Export dataframe to csv
#df1_cfmid_both.to_csv('Approach 3 ENTACT compound cfmid results.csv', sep=',', encoding='utf-8', index=False)


print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))