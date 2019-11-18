# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 13:26:54 2019

@author: AChao
"""

############################################################################
# This script reads in three separate experimental CE CFMID searched files, and then sums 
# the scores for DTXCID's found in multiple energies into one score
############################################################################


import os
import time
import pandas as pd
start_time = time.time()

sheet_list = ['499', '500', '501', '502', '503', '504', '505', '506', '507', '508']
#sheet_list = ['499']



########################################################################
# Read in CFMID search results - One Score
########################################################################

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
    
    
df2_pos = df2.loc[df2['cfmid_mode'] == 'Esi+'] # Separate out CFMID results by ionization mode
df2_neg = df2.loc[df2['cfmid_mode'] == 'Esi-']

df2_pos_10 = df2_pos.loc[df2_pos['CE'] == 10] # Separate the results by the individual energy experiments
df2_pos_20 = df2_pos.loc[df2_pos['CE'] == 20]
df2_pos_40 = df2_pos.loc[df2_pos['CE'] == 40]
df2_neg_10 = df2_neg.loc[df2_neg['CE'] == 10]
df2_neg_20 = df2_neg.loc[df2_neg['CE'] == 20]
df2_neg_40 = df2_neg.loc[df2_neg['CE'] == 40]

df2_pos_10.rename(columns={'SCORE':'score_10'}, inplace=True) #Give each score column a name indicating the energy so these can be merged
df2_pos_20.rename(columns={'SCORE':'score_20'}, inplace=True) 
df2_pos_40.rename(columns={'SCORE':'score_40'}, inplace=True) 
df2_neg_10.rename(columns={'SCORE':'score_10'}, inplace=True)
df2_neg_20.rename(columns={'SCORE':'score_20'}, inplace=True)
df2_neg_40.rename(columns={'SCORE':'score_40'}, inplace=True)


df2_pos_merge = pd.merge(df2_pos_10, df2_pos_20[['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode', 'score_20', 'MASS_in_MGF']], 
                         how='outer', on = ['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode'])
df2_pos_merge = pd.merge(df2_pos_merge, df2_pos_40[['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode', 'score_40', 'MASS_in_MGF']], 
                         how='outer', on = ['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode'])

# The below is because sometimes a compound is found at 40 and not in 10, but we need to capture its mgf mass information.
# Without these lines, the mgf mass info gets lost because it merges on 10, for which in some cases there is no mgf mass information.
df2_pos_merge['test'] = df2_pos_merge[['MASS_in_MGF', 'MASS_in_MGF_x', 'MASS_in_MGF_y']].mean(axis=1)

df2_pos_merge = df2_pos_merge.drop(['MASS_in_MGF', 'MASS_in_MGF_x', 'MASS_in_MGF_y', 'CE', 'RANK'], 1)
df2_pos_merge.rename(columns={'test':'MASS_in_MGF'}, inplace=True) 



df2_neg_merge = pd.merge(df2_neg_10, df2_neg_20[['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode', 'score_20', 'MASS_in_MGF']], 
                         how='outer', on = ['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode'])
df2_neg_merge = pd.merge(df2_neg_merge, df2_neg_40[['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode', 'score_40', 'MASS_in_MGF']], 
                         how='outer', on = ['DTXCID', 'MASS', 'FORMULA', 'MATCHES', 'mixture', 'cfmid_mode'])

df2_neg_merge['test'] = df2_neg_merge[['MASS_in_MGF', 'MASS_in_MGF_x', 'MASS_in_MGF_y']].mean(axis=1)

df2_neg_merge = df2_neg_merge.drop(['MASS_in_MGF', 'MASS_in_MGF_x', 'MASS_in_MGF_y', 'CE', 'RANK'], 1)
df2_neg_merge.rename(columns={'test':'MASS_in_MGF'}, inplace=True)

df2_pos_merge['SCORE'] = df2_pos_merge[['score_10', 'score_20', 'score_40']].sum(axis='columns') # Sum the three energy scores into one score
df2_pos_merge = df2_pos_merge.drop(['score_10', 'score_20', 'score_40'], 1) # drop the individual scores

df2_neg_merge['SCORE'] = df2_neg_merge[['score_10', 'score_20', 'score_40']].sum(axis='columns') # Sum the three energy scores into one score
df2_neg_merge = df2_neg_merge.drop(['score_10', 'score_20', 'score_40'], 1) # drop the individual scores

df2_final = pd.concat([df2_pos_merge, df2_neg_merge])

df2_final['MASS_in_MGF'] = df2_final['MASS_in_MGF'].apply(lambda x:round(x,3))

'''
os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/Results code")
df2_final.to_csv('ENTACT_CFMID_all_mixtures_merged_results_rounded.csv', sep=',', encoding='utf-8', index=False)
'''

print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))