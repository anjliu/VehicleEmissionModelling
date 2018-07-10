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
periodsused=['morning1','morning2','morning3','afternoon1','afternoon2']

# columns in moves output to use
cols=['linkID','pollutantID','processID','sourceTypeID','fuelTypeID','modelYearID','emissionQuant']
colsrename={'sourceType':'Vehicle Type'}

# moves output ID decoders
decpath='MOVESdecoders'
decoders={}
for file in os.listdir(decpath):
    decoders[file[:-4]]=pd.read_csv(decpath+'/'+file,index_col='ID').iloc[:,0]
del decoders['sourceType']

# link number decoders
ldecpath='linkDecoders'
ldecoders={}
for file in os.listdir(ldecpath):
    ldecoders[file[:-4]]=pd.read_csv(ldecpath+'/'+file,index_col='ID').iloc[:,0]
# link digits and their meaning
ldigit=pd.read_csv('linkDigits.csv',index_col='digit').iloc[:,0]

mo_ft=pd.DataFrame()

for p in periodsused:
    # replace IDs with names    
    movesOutputPath='MOVESoutputs/thru/'+str(p)+'MovesOutput.csv'
    mo=pd.read_csv(movesOutputPath,usecols=cols)
    mo.columns = mo.columns.str.replace('ID','')
    for d in decoders:
        decoder=decoders[d]
        for i in decoder.index:
            mo[d]=mo[d].replace(i,decoder[i])
    mo['period']=p
    
    mo_ft=mo_ft.append(mo)
    
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

inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']

monew=pd.DataFrame()
moperiod=[]
mointer=[]
modirec=[]
mocycle=[]
motype=[]
moemi=[]
mopol=[]
mocols=[moperiod,mointer,modirec,mocycle,motype,mopol,moemi]
mocolnames=['period','inter','direc','cycle','type','pollutant','emi']

mo_ft=parselinkno(mo_ft)
mo_ft['cycle']=mo_ft['link']%100

pollutants=mo_ft['pollutant'].unique()


for p in periodsused:
    mop=mo_ft.loc[mo_ft['period']==p]
    for i in inters:
        mopi=mop.loc[mop['intersection']==i]
        for d in direcs:
            mopid=mopi.loc[mopi['approachDirection']==d]
            cycles=mopid['cycle'].unique()
            for c in cycles:
                mopidc=mopid.loc[mopid['cycle']==c]
                types=mopidc['sourceType'].unique()
                for t in types:
                    mot=mopidc.loc[mopidc['sourceType']==t]
                    for pol in pollutants:
                        mopol=mot.loc[mot['pollutant']==pol]
                        emi=sum(mopol['emissionQuant'])
                    
                        toappend=[p,i,d,c,t,pol,emi]
                        n=0
                        for x in mocols:
                            x.append(toappend[n])
                            n+=1
n=0
for x in mocols:
    monew[mocolnames[n]]=x
    n+=1
savename='MovesOutputDecoded-thru'
monew.to_csv(savename+'.csv',index=False)