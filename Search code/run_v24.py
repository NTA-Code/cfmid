# -*- coding: utf-8 -*-
"""
Created on Tue May 01 12:34:34 2018

@author: HALGhoul
"""

import mgf_parser_v24 as mg
import time
import os

# Directory in which exported mgf files are stored
os.chdir("C:/users/achao/OneDrive - Environmental Protection Agency (EPA)/profile/Desktop/pythontemp/CFMID code for manuscript submission/")


# Filename of the exported mgf file to read, process, match against CFM-ID database, and score
filename = '499_neg_C_01'


# This snippet of code processes the mgf file into a csv file
# De-duplicates spectra by precursor mass, keeping the spectrum with the highest sum intensity of fragments
mgfile=filename+'.mgf'
mg.parseMGF(mgfile)


# Generate the filename of the csv to read
file = filename+'.csv'
fpcdl = os.getcwd()+'/499 master spiked list.csv' # Placeholder file- functionality removed


# The below code searches the CFMID database for the spectra contained within the csv file
# Parameters:
# file- name of the csv file
# filename- name for outputting purposes
# fpcdl- Ignore (removed feature)
# 10 - The accuracy window (in ppm) for matching precursor mass (from csv) to potential candidate compounds in CFMID database
# 0.02 - Accuracy window (in Da) for matching fragment ions in CFMID spectrum to experimental spectrum (csv)
# POSMODE - Select whether data is in positive mode (= True) or negative mode (= False)
# filtering - Ignore (removed feature)
# energy - Ignore (removed feature)
# bymass - Ignore (removed feature)
# pcdl_mode - Ignore (removed feature)
t0=time.clock()
mg.compare_mgf_df(file,filename,fpcdl,10,0.02,POSMODE=False,filtering=False,energy=False,bymass=True,pcdl_mode=False)
t1=time.clock()
print("time to Process is " + str(t1-t0))


# Final output is in two files:
# One with the suffix "_CFMID_Multiscores_ALLHits.xlsx" - This contains individual scores against CFMID at CE 10, 20, and 40
# One with the suffix "_CFMID_Onescore_ALLHits.xlsx" - This contains the CFMID scores summed across CFMID CE's into one score (i.e. Approach 2 and 3)