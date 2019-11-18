# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:11:13 2017

@author: Hussein Al Ghoul
"""

import os
import pandas as pd
import numpy as np

import CosineDotProduct_v24 as cpd
import pymysql as mysql
import time
#from sqlalchemy import create_engine
#cnx = create_engine('mysql://root:password@localhost/db')


dfc = None
dfn = None

Adduct_Mass = 1.007825


def read_NTA_data(file=''):  # read a csv file into a DataFrame
    global dfn
    dfn = pd.read_csv(file)
    #dfn.columns = dfn.columns.str.replace(' ','_')
    energy_column = 'Collision Energy A (Fragment Mass, Predicted Formula, Abundance, Abundance %)'
    regexI = "[\^,]*([0-9]+\.*[0-9]*)$" # a regex to find the score looking for a pattern of "=something_to_find ]" 
    regexE = "(.*?)\,{2}"
    dfn['INTENSITY0N'] = dfn[energy_column].str.extract(regexI,expand=True).astype('float64')
    dfn['PMASS'] = dfn[energy_column].str.extract(regexE,expand=True).astype('float64')
    dfn = dfn[['CASRN','PMASS','INTENSITY0N']]
    dfn.sort_values(['CASRN','INTENSITY0N'],ascending=[True,False],inplace=True)    
    dfn.to_csv("nta.csv",index=False)
    #print dfn
    return dfn





def parseMGF(file=''):
    NFile = file.rsplit('/',1)[-1]
    NewFile = NFile.rsplit('.',1)[0] + ".csv"
    with open(file) as f:
        RESULT = list()
        for line in f:
            if line.startswith('TITLE'):
                result = list()
                fields = line.split(' ')
                title, MS, of, pmass, charge, at, RT, mins, delimeter = fields
                #print(fields)   # AC
                #print (pmass, charge, RT)
                #result.append([pmass,RT])
            if line.startswith('RTINSECONDS'):
                RTS = line.split('=')[1]
                for line in f:
                    if line.split(' ')[0] == 'END':
                        break
                    a, b  = line.split('\t')
                    result.append([float(pmass), float(RT), charge, float(a),float(b)])
                RESULT.append(result)
        #print RESULT[0]
    categories = [ "RUN %s" %i for i in range(0,len(RESULT))]
    dfg = pd.concat([pd.DataFrame(d) for d in RESULT], keys=categories)
    dfg.columns = ["MASS", "RETENTION TIME", "CHARGE", "PMASS_y","INTENSITY"]
    dfg.sort_values(['MASS','RETENTION TIME'],ascending=[True,True],inplace=True) 
    df1 = dfg.reset_index()
    df1['TOTAL INTENSITY'] = df1.groupby(['MASS','RETENTION TIME'])['INTENSITY'].transform(sum)
    df1.sort_values(['MASS','TOTAL INTENSITY'],ascending=[True,True],inplace=True)
    df1 = df1.groupby('MASS').apply(lambda x: x[x['TOTAL INTENSITY'] == x['TOTAL INTENSITY'].max()])
    df1.loc[df1['PMASS_y']>df1['MASS'],'INTENSITY']=None
    df1.dropna(inplace=True)
    df1.sort_values(['MASS','INTENSITY'],ascending=[True,False],inplace=True)
    #df1['INTENSITY0M'] = df1.groupby('MASS')['INTENSITY'].apply(lambda x: (x/x.nlargest(2).min())*100.0)
    df1['INTENSITY0M'] = (df1['INTENSITY']/df1.groupby('MASS')['INTENSITY'].transform(np.max))*100.0
    df1.loc[df1['INTENSITY0M']>100,'INTENSITY0M']=None
    #df1.loc[df1['INTENSITY0M']<0.1,'INTENSITY0M']=None
    df1.reset_index(drop=True, inplace=True) # reset index
    #df1['ENERGY'] = 'energy0'
    df1.to_csv(NewFile,float_format='%.5f',index=False)
    dfg = df1
    #dfg.to_csv("CE10d_mgf.csv",index=False)
    return dfg


def spectrum_reader(file=''):
    dfg = pd.read_csv(file)
    #dfg = dfg.groupby(['MASS','RETENTION TIME']).head(30)
    dfg = dfg[(dfg['INTENSITY0M']<=100) & (dfg['INTENSITY0M']>0.0)]
    return dfg


''' A SQL query to get all the corresponding info from the database'''
def sqlCFMID(mass=None,ppm=None,mode=None,formula=None):
    db = mysql.connect(host="host name",
                   user="user name",
                   passwd="password",
                   db="db name")      
                       
    cur = db.cursor()
    accuracy_condition = ''
    if mass:
        if ppm>=1:
            accuracy_condition = """(abs(c.mass-"""+str(mass)+""")/"""+str(mass)+""")*1000000<"""+str(ppm)
        if ppm<1 and ppm>0: 
            accuracy_condition = """(abs(c.mass-"""+str(mass)+""")"""
    if formula:
        accuracy_condition = """c.formula='"""+formula+"""'"""
    query= """select t1.dtxcid as DTXCID, t1.formula as FORMULA,t1.mass as MASS, t1.mz as PMASS_x, (t1.intensity/maxintensity)*100.0 as INTENSITY0C,
t1.energy as ENERGY 
from 

(select c.dtxcid, max(p.intensity) as maxintensity, p.energy from peak p
			Inner Join job j on p.job_id=j.id
			Inner Join spectra s on j.spectra_id = s.id
			Inner Join chemical c on j.dtxcid=c.dtxcid
			#Inner Join fragtool ft on j.fragtool_id=ft.id
     		#inner join fragintensity fi on fi.peak_id = p.id       
		    where """ +accuracy_condition + """ 
       and s.type='""" + mode + """'
			group by c.dtxcid, p.energy) t2
    left join
	(select c.*, p.* from peak p
		Inner Join job j on p.job_id=j.id
		Inner Join spectra s on j.spectra_id = s.id
		Inner Join chemical c on j.dtxcid=c.dtxcid
		#Inner Join fragtool ft on j.fragtool_id=ft.id   
    	#inner join fragintensity fi on fi.peak_id = p.id 
		where """ +accuracy_condition + """ 
       and s.type='""" + mode + """') t1

on t1.dtxcid=t2.dtxcid and t1.energy=t2.energy
order by DTXCID,ENERGY,INTENSITY0C desc;"""
    #Decided to chunk the query results for speed optimization in post porocessing (spectral matching)
    cur.execute(query)
    chunks=list()
    for chunk in pd.read_sql(query,db,chunksize=1000):
        chunks.append(chunk)
    return chunks
            

def list_maker(fpcdl,dfg,mode,pcdl_mode):
    # make a list of the MGF masses corresponding to the PCL Monoisotopic masses 
    dfpc = pd.read_csv(fpcdl)
    dfg['Mass_rounded'] = dfg['MASS'].round(1)
    dfpcdl = dfpc.loc[dfpc['Polarity'] == mode]
    dfpcdl['Mass_rounded'] = dfpcdl['Neutral Monoisotopic Mass'].round(1)
    df = pd.merge(dfpcdl,dfg,how='left',on='Mass_rounded') 
    df['MATCHES'] = np.where((((abs(df['Neutral Monoisotopic Mass']-df['MASS'])/df['MASS'])*1000000)<=15),'1','0') 
    df.drop(df[df['MATCHES'] == '0'].index,inplace=True)
    df=df[np.isfinite(df['MASS'])]
    if pcdl_mode:
        mlist = df['MASS'].unique().tolist()
    else:
        mlist = dfg['MASS'].unique().tolist()
    print(mlist)
    return mlist


def input_parser(file,filtering,energy):
    with open(file) as f:
        RESULT = list()
        for line in f:
            if line.startswith('Mode'):
                result = list()
                tolerance = list()
                mode = line.strip().split('=')[1]
                print(mode)
            if line.startswith('Mass'):
                mass = line.split('=')[1]
            if line.startswith('Energy'):
                Energy = line.split('=')[1].strip('\n')                
            if line.startswith('MTolerance'):
                ppm = line.split('=')[1]
            if line.startswith('PTolerance'):
                ppm_sl = line.split('=')[1]
            if line.startswith('#BEGIN'):
                for line in f:
                    if line.split(' ')[0] == '#END':
                        break
                    a, b  = line.split('\t')
                    result.append([float(mass), str(mode), str(Energy), float(a),float(b)])
                    tolerance.append(float(ppm))
                    tolerance.append(float(ppm_sl))
                RESULT.append(result)
    print(RESULT[0])
    categories = [ "RUN %s" %i for i in range(0,len(RESULT))]
    df = pd.concat([pd.DataFrame(d) for d in RESULT], keys=categories)
    df.columns = ["MASS", "MODE", "ENERGY", "PMASS_y","INTENSITY"]
    df.sort_values(['MASS','PMASS_y'],ascending=[True,True],inplace=True) 
    df.loc[df['PMASS_y']>df['MASS'],'INTENSITY']=None
    df.dropna(inplace=True)
    #df1['INTENSITY0M'] = df1.groupby('MASS')['INTENSITY'].apply(lambda x: (x/x.nlargest(2).min())*100.0)
    df['INTENSITY0M'] = (df['INTENSITY']/df.groupby('MASS')['INTENSITY'].transform(np.max))*100.0
    df.loc[df['INTENSITY0M']>100,'INTENSITY0M']=None
    #df1.loc[df1['INTENSITY0M']<0.1,'INTENSITY0M']=None
    df.reset_index(drop=True, inplace=True) # reset index
    df.to_csv('try.csv',float_format='%.5f',index=False)
    ppm = tolerance[0]
    ppm_sl = tolerance[1]    
    compare_df(df,ppm,ppm_sl,filtering,energy)
    return df    
    
    
def compare_df(dfg,ppm,ppm_sl,filtering,energy):
    mode = dfg.at[0,'MODE'] 
    if mode=='ESI-MSMS-pos':
        dfg['MASS'] = dfg.groupby('MODE')['MASS'].transform(lambda x: (x-1.007825))
        dfg['MASS'] = dfg['MASS'].round(6)
    if mode=='ESI-MSMS-neg':
        dfg['MASS'] = dfg.groupby('MODE')['MASS'].transform(lambda x: (x+1.007825)) 
        dfg['MASS'] = dfg['MASS'].round(6)    
    if mode=='ESI-MSMS-neutral':
        dfg['MASS'] = dfg['MASS'].round(6)   
    mass = dfg.at[0,'MASS']        
    dfcfmid = sqlCFMID(mass,ppm,mode)
    print(dfg)
    #print dfcfmid
    dfmgf = None
    df = None
    dfmgf = dfg[dfg['MASS'] == mass].reset_index()
    dfmgf.sort_values('MASS',ascending=True,inplace=True)
    df = cpd.Score(dfcfmid,dfmgf,mass,ppm_sl,filtering,energy)
    dfAE = df[0] #all energies scores
    dfS = df[1]
    dfAE.to_excel('Energies_Scores_input_spectrum.xlsx',engine='xlsxwriter')
    dfS.to_excel('Sum_Scores_input_spectrum.xlsx',engine='xlsxwriter')



def compare_mgf_df(file,filename,fpcdl,ppm,ppm_sl,POSMODE,filtering,energy,bymass,pcdl_mode):    
    dfg = spectrum_reader(file)
    if POSMODE:
        mode='ESI-MSMS-pos'
        polarity=['ESI+','Esi+']
        #CMass = Mass - Adduct_Mass
        dfg['MASS'] = dfg.groupby('RETENTION TIME')['MASS'].transform(lambda x: (x-1.007825))
        dfg['MASS'] = dfg['MASS'].round(6)
    else:
        mode='ESI-MSMS-neg'
        polarity=['ESI-','Esi-']
        #CMass = Mass + Adduct_Mass
        dfg['MASS'] = dfg.groupby('RETENTION TIME')['MASS'].transform(lambda x: (x+1.007825)) 
        dfg['MASS'] = dfg['MASS'].round(6)
    #dfg.to_csv("dfg_aftercsv.csv",float_format='%.7f',index=False)  
    #mass_list = dfg['MASS'].unique().tolist()
    #mass_list = [312.184525]
    Polarity = polarity[0]
    mass_list = list_maker(fpcdl,dfg,Polarity,pcdl_mode)
    if not mass_list:
        Polarity = polarity[1]
        mass_list = list_maker(fpcdl,dfg,Polarity,pcdl_mode)
    #mass_list = mass_list[1:2]
    print(mass_list)
    print("Number of masses to search: " + str(len(mass_list)))
    dfAE_list=list()
    dfS_list=list()  
    for mass in mass_list:
        index = mass_list.index(mass) + 1
        print("searching mass " + str(mass) + " number " + str(index) + " of " + str(len(mass_list)))
        dfcfmid = sqlCFMID(mass,ppm,mode)
        if not dfcfmid:
            print("No matches for this mass in CFMID library, consider changing the accuracy of the queried mass")
        else:    
            dfmgf = None
            df = None
            dfmgf = dfg[dfg['MASS'] == mass].reset_index()
            dfmgf.sort_values('MASS',ascending=True,inplace=True)
            df = cpd.Score(dfcfmid,dfmgf,mass,ppm_sl,filtering,energy)
            if mode=='ESI-MSMS-pos':
                df[0]['MASS_in_MGF'] = mass + 1.007825
                df[1]['MASS_in_MGF'] = mass + 1.007825
            if mode=='ESI-MSMS-neg':
                df[0]['MASS_in_MGF'] = mass - 1.007825
                df[1]['MASS_in_MGF'] = mass - 1.007825
            dfAE_list.append(df[0]) #all energies scores
            dfS_list.append(df[1]) #sum of all energies score
    if not  dfAE_list:
        print("No matches All Energies found")
    else:
        dfAE_total = pd.concat(dfAE_list) #all energies scores for all matches

    if not dfS_list:
        print ("No matches Single Energies found")
    else:        
        dfS_total = pd.concat(dfS_list) #Sum of scores for all matches
    #th_dtxcid = dfAE_list[0].at[0,'DTXCID'] #top hit dtxcid for plotting
    #print(th_dtxcid) 
    #dfcfm = pd.concat(dfcfmid)[(pd.concat(dfcfmid)['DTXCID'] == th_dtxcid) & (pd.concat(dfcfmid)['ENERGY'] == 'energy2')].reset_index()
    #cpd.plot(dfcfm,dfmgf)
    if pcdl_mode:
        df_resultAE = merge_pcdl(fpcdl,dfAE_total,Polarity,bymass,ppm)
        NFile = fpcdl.rsplit('/',1)[-1]
        firstFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_MultiScores_byformula_DTXCID.xlsx"
        if bymass:
            secondFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_MultiScores_bymass_AllHits.xlsx"    
        else:
            secondFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_MultiScores_byformula_AllHits.xlsx"    
        df_resultAE[0].to_excel(firstFile,engine='xlsxwriter')
        df_resultAE[1].to_excel(secondFile,engine='xlsxwriter')
        df_resultS = merge_pcdl(fpcdl,dfS_total,Polarity,bymass,ppm)

        NFile = fpcdl.rsplit('/',1)[-1]
        firstFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_OneScore_wformula_DTXCID.xlsx"
        if bymass:    
            secondFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_OneScore_bymass_AllHits.xlsx" 
        else:
            secondFile = NFile.rsplit('.',1)[0] + "_compared_with_CFMID_OneScore_byformula_AllHits.xlsx" 
        df_resultS[0].to_excel(firstFile,engine='xlsxwriter')
        df_resultS[1].to_excel(secondFile,engine='xlsxwriter')    
    else:
        #dfAE_total.to_excel('CFMID_Multiscores_AllHits.xlsx',engine='xlsxwriter')
        #dfS_total.to_excel('CFMID_OneScore_AllHits.xlsx',engine='xlsxwriter')
        dfAE_total.to_excel(filename+'_CFMID_Multiscores_AllHits.xlsx',engine='xlsxwriter')
        dfS_total.to_excel(filename+'_CFMID_OneScore_AllHits.xlsx',engine='xlsxwriter')




"""///This function merges the results from the CFMID analysis back to the
PCDL file. It first merges them entirely by the dtxcid, then either by formula matches and in that case only
formula matches are merged back, or by mass in ncase no formula prediction was done
of the pcdl part. I added a filter to give you an On/Off switch to control whether to merge
by formula or by mass. Set bymass=True to m,erge by mass"""
def merge_pcdl(fpcdl,df,polarity,bymass,ppm):
    dfT = list()
    df['Polarity'] = polarity
    #print df
    dfpcdl = pd.read_csv(fpcdl)
    dfpcdl['Predicted/Matched Formula'] = dfpcdl['Predicted/Matched Formula'].str.replace(' ','')
    dfm = pd.merge(dfpcdl,df,how='left',on='DTXCID')
    dfT.append(dfm)
    dfm = None

    if bymass:
        dfpcdl['MASS_rounded'] = dfpcdl['Neutral Monoisotopic Mass'].round(1)
        df['MASS_rounded'] = df['MASS'].round(1)
        dfm = pd.merge(dfpcdl,df,suffixes=['','_CFMID'],how='left',on='MASS_rounded')        
        dfm.drop(['MASS_rounded','PCDL MATCH','RT1 (min)'],axis=1,inplace=True)
        #dfm = dfm.set_index(['Neutral Monoisotopic Mass','Predicted/Matched Formula','Predicted/Matched Compound'	,'DTXSID','CASRN','DTXCID','MASS','FORMULA','DTXCID_CFMID']).head()

    else:
        dfm = pd.merge(dfpcdl,df,suffixes=['','_CFMID'],how='left',left_on='Predicted/Matched Formula',right_on='FORMULA')
        #dfm.drop(['PCDL MATCH','RT1 (min)'],axis=1,inplace=True)
        dfm=dfm.dropna(axis=1,how='all')
        #dfm.set_index(dfm.columns.values.tolist(),inplace=True)
        print(dfm)
        #dfm = dfm.set_index(dfm.columns.values.tolist()[0:11]).head()


    #dfm.sort_values(['Neutral Monoisotopic Mass','MASS','RANK'],ascending=[True,True,True],inplace=True)
    dfT.append(dfm)
    return dfT



def indexing(file):
    df = pd.read_csv(file)
    df.set_index(df.columns.values.tolist(),inplace=True)
    print(df)
    df.to_excel('try_indexing_outcome.xlsx',engine='xlsxwriter')

        


    
#read_NTA_data('/home/hussein/Documents/NTA/Python_alt/ENTACT_DataReporting_EPA_MS2.csv')
#parseCFMID('/home/hussein/Documents/NTA/Python_alt/spectra_ESI-MSMS-neg_mass.dat')
#compare_df(183.057312)
#compare_df(183.058217) 
    
#parseMGF(os.getcwd()+'/20180418_505_CE40.mgf') #<--Convert MGF to CSV

# to read a signle spectrum from text file input spectrum:
#input_parser('input_spectrum.txt') 


# to process MGD files uncomment the following lines
#indexing('try_indexing.csv')
'''
file = os.getcwd()+'/20180418_505_CE40.csv'
fpcdl = os.getcwd()+'/20180419 505 pos CE40_PCDL.csv'
t0=time.clock()
compare_mgf_df(file,fpcdl,10,0.02,POSMODE=True)
t1=time.clock()
print ("time to Process is " + str(t1-t0))
'''








