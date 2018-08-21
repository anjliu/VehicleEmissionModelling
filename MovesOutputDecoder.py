# -*- coding: utf-8 -*-
"""
Created on Sat May 26 12:05:56 2018

@author: Anjie
"""

import pandas as pd
import json, os

# manual input:
scenario='TrajectoryReconstruction-wIdle'


inputpath='MOVESoutputs/'+scenario+'/'
outputpath='MOVESoutputsDecoded/'

# get run specs: periods
runSpecs=json.load(open('runSpecs-limited.json')) # use limited runspecs - not all periods are available
periods=runSpecs['periods']

# columns in moves output to use
cols=['linkID','pollutantID','processID','sourceTypeID','fuelTypeID','modelYearID','emissionQuant']
colsrename={'sourceType':'Vehicle Type'}

# columns in links to use
linkcols=[0,4,5,6]
linkcolsrename={'linkID':'link','linkLength':'length','linkVolume':'volume','linkAvgSpeed':'averageSpeed'}

# moves output ID decoders
decpath='MOVESdecoders'
decoders={}
for file in os.listdir(decpath):
    decoders[file[:-4]]=pd.read_csv(decpath+'/'+file,index_col='ID').iloc[:,0]
# link number decoders
ldecpath='linkDecoders'
ldecoders={}
for file in os.listdir(ldecpath):
    ldecoders[file[:-4]]=pd.read_csv(ldecpath+'/'+file,index_col='ID').iloc[:,0]
# link digits and their meaning
ldigit=pd.read_csv('linkDigits.csv',index_col='digit').iloc[:,0]

mo_decoded=pd.DataFrame()
links_ft=pd.DataFrame()

for p in periods:
    # replace IDs with names    
    movesPath=inputpath+str(p)+'MovesOutput.csv'
    mo=pd.read_csv(movesPath,usecols=cols)
    mo.columns = mo.columns.str.replace('ID','')
    for d in decoders:
        decoder=decoders[d]
        for i in decoder.index:
            mo[d]=mo[d].replace(i,decoder[i])
    mo['hour']=periods[p]['hour']
    
    mo_decoded=mo_decoded.append(mo)
    
    print('Processed ' +str(p))
    
mo_decoded['Emissions [tonnes]']=mo_decoded['emissionQuant']/10e6

    
def parselinkno(df):
    
    for d in ldecoders:
        ldecoder=ldecoders[d]
        n=ldigit[ldigit==d].index[0]
        df[d]=df['link']%10**(n)/10**(n-1)
        df[d]=df[d].astype(int)
        for i in ldecoder.index:
            df[d]=df[d].replace(i,ldecoder[i])
        print('Decoded links for ' +str(d))
    return df

savefiles=[mo_decoded,links_ft]
savenames=['MovesOutputForTableau','LinksForTableau']

mo_decoded=parselinkno(mo_decoded)
mo_decoded['cycle']=mo_decoded['link']%100
mo_decoded.to_csv(outputpath+scenario+'.csv',index=False)
