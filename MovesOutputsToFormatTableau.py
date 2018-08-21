# -*- coding: utf-8 -*-
"""
Created on Sat May 26 12:05:56 2018

@author: Anjie
"""

import pandas as pd
import json, os

# get run specs: periods
runSpecs=json.load(open('runSpecs.json'))
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

# file paths
movesInputPath='MOVESinputs/all/'
savepath='ForTableau/'

mo_ft=pd.DataFrame()
links_ft=pd.DataFrame()

for p in periods:
    # replace IDs with names    
    movesOutputPath='MOVESoutputs/allMovements/'+str(p)+'MovesOutput.csv'
    mo=pd.read_csv(movesOutputPath,usecols=cols)
    mo.columns = mo.columns.str.replace('ID','')
    for d in decoders:
        decoder=decoders[d]
        for i in decoder.index:
            mo[d]=mo[d].replace(i,decoder[i])
    mo['hour']=periods[p]['hour']
    
    mo_ft=mo_ft.append(mo)
    
    # add traffic data from MOVES inputs
    linkspath=movesInputPath+str(p)+'/links.xlsx'
    links=pd.read_excel(linkspath,usecols=linkcols)
    links=links.rename(columns={'linkID':'link'})
    links_ft=links_ft.append(links)
    
    print('Processed ' +str(p))
    
mo_ft['Emissions [tonnes]']=mo_ft['emissionQuant']/10e6

    
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

savefiles=[mo_ft,links_ft]
savenames=['MovesOutputForTableau','LinksForTableau']

file=0
for dfout in savefiles:
    dfout=parselinkno(dfout)
    dfout['cycle']=dfout['link']%100
    dfout.to_csv(savepath+savenames[file]+'.csv',index=False)
    file+=1