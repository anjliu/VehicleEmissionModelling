# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 18:17:23 2018

@author: Anjie
"""
import json
import pandas as pd

# inputs
rates=pd.read_csv('rates.csv',names=['typeModePolProcess','rate'])
scenario='test-temp'
inputpath='MOVESinputs/'+scenario
opm=pd.read_excel(inputpath+'/OpModeDistribution.xlsx')
ls=pd.read_excel(inputpath+'/links.xlsx')
links=ls['linkID']
lst=pd.read_excel(inputpath+'/linkSourceTypes.xlsx')

# specs
## vehicle types
typeSpecs=pd.read_csv('MOVESdecoders/sourceType.csv',index_col=0)
types=typeSpecs.index
## operating modes
runSpecs=json.load(open('runSpecs.json'))
modes=runSpecs['opModeFormat']['opModeSets']['modeSet']
# pollutant processes
pps=[101,501,515,601,9001,9101]

# variables
totalvolume=ls['linkVolume'][0]

# get volume for each vehicle type
typeVolumes=[]
for t in types:
    ratio=lst.loc[lst['sourceTypeID']==t,'sourceTypeHourFraction'].iloc[0]
    typeVolumes.append(ratio*totalvolume)
typeSpecs['Volume']=typeVolumes

fractions=opm.loc[opm['polProcessID']==9101]

out_emissions=[]
out_type=[]
out_pollutant=[]
out_links=[]

for l in links:
    opm_l=opm.loc[opm['linkID']==l]
    for t in types:
        opm_t=opm_l.loc[opm_l['sourceTypeID']==t]
        fractions_t=fractions.loc[fractions['sourceTypeID']==t]
    
        for p in pps:
            out_pollutant.append(p)
            opm_tp=opm_t.loc[opm_t['polProcessID']==p]
    
            emissionsPerMode=[]
            for m in modes:
                opm_tpm=opm_tp.loc[opm_tp['opModeID']==m]
                fractions_tm=fractions_t.loc[fractions_t['opModeID']==m]
                fraction=fractions_tm['opModeFraction'].iloc[0]
    
                tmp=int(str(t*100+m)+str(p))
                tmp_rate=rates.loc[rates['typeModePolProcess']==tmp,'rate']
                tmp_rate=tmp_rate.iloc[0]
                emissionsPerMode.append(fraction*tmp_rate*typeSpecs['Volume'][t])
                
            out_type.append(t)
            out_emissions.append(sum(emissionsPerMode))
            out_links.append(l)
            
emissions=pd.DataFrame({'link':out_links,
                        'type':out_type,
                        'pollutant':out_pollutant,
                        'emissions':out_emissions})