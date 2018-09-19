# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 18:17:23 2018

@author: Anjie
"""
import json, os
import pandas as pd

from time import time

start_t=time()

# manual inputs
scenario='advDetector100m'

# constant inputs
periods=['morning1']#,'morning2','morning3']#,'afternoon1','afternoon2','afternoon3','evening1','evening2']
rates=pd.read_csv('rates.csv',names=['typeModePolProcess','rate'])
inputpath='MOVESinputs/'+scenario
savepath='RealTimeEstimates/'+scenario
if not os.path.exists(savepath):
    os.makedirs(savepath)

# specs
## vehicle types
typeSpecs=pd.read_csv('MOVESdecoders/sourceType.csv',index_col=0)
types=typeSpecs.index
## operating modes
runSpecs=json.load(open('runSpecs.json'))
modes=runSpecs['opModeFormat']['opModeSets']['modeSet']
# pollutant processes
pps=[101,501,515,601,9001,9101]


for p in periods:
    start_p=time()
    inputpath_p=inputpath+'/'+p
    savepath_p=savepath+'/'+p
    
    opm=pd.read_excel(inputpath_p+'/OpModeDistribution.xlsx')
    ls=pd.read_excel(inputpath_p+'/links.xlsx')
    links=ls['linkID']
    lst=pd.read_excel(inputpath_p+'/linkSourceTypes.xlsx')
    
    out_emissions=[]
    out_type=[]
    out_pollutant=[]
    out_links=[]
    
    for l in links: # for each link
        totalvolume=ls.loc[ls['linkID']==l,'linkVolume'].iloc[0] # link volume
        # get volume for each vehicle type
        typeVolumes=[]
        for t in types:
            ratio=lst.loc[lst['sourceTypeID']==t,'sourceTypeHourFraction'].iloc[0]
            typeVolumes.append(ratio*totalvolume)
        typeSpecs['Volume']=typeVolumes
            
        opm_l=opm.loc[opm['linkID']==l] # op mode distribution for this link
        fractions=opm_l.loc[opm_l['polProcessID']==9101] # filter to fractions for one polprocess (same for all polprocesses)
    
        for t in types:
            opm_t=opm_l.loc[opm_l['sourceTypeID']==t] # op mode distribution for link and vehicle type
            
            if not opm_t.empty: # if there are vehicle types t on this link
                fractions_t=fractions.loc[fractions['sourceTypeID']==t] # fraction of this vehicle type on this link
            
                for pp in pps: # for each pollutant process
                    out_pollutant.append(pp) # append to pollutants column
                    opm_tp=opm_t.loc[opm_t['polProcessID']==pp] # op mode distribution for vehicle type and pollutant
            
                    emissionsPerMode=[] # start list for emissions of each operating mode
                    for m in modes: # for each mode
                        opm_tpm=opm_tp.loc[opm_tp['opModeID']==m] # op mode distribution for type, pollutant process, and mode
                        fractions_tm=fractions_t.loc[fractions_t['opModeID']==m] # op mode fraction for this vehicle type
                        fraction=fractions_tm['opModeFraction'].iloc[0] # get value from cell
            
                        tmp=int(str(t*100+m)+str(pp)) # type-mode-polprocess ID
                        tmp_rate=rates.loc[rates['typeModePolProcess']==tmp,'rate'] # emission rate for this type-mode-polprocess
                        tmp_rate=tmp_rate.iloc[0] # get value from cell
                        emissionsPerMode.append(fraction*tmp_rate*typeSpecs['Volume'][t]) # emissions for this type, polprocess, and mode
                        
                    out_type.append(t)
                    out_emissions.append(sum(emissionsPerMode)) # add up emissions of all modes for this vehicle type and polprocess
                    out_links.append(l)
    
#        print('Processed link '+str(l))
    
    # output dataframe
    emissions=pd.DataFrame({'link':out_links,
                            'type':out_type,
                            'pollutant':out_pollutant,
                            'emissions':out_emissions})
    
    emissions.to_csv(savepath_p+'.csv',index=False)
    time_p=time()-start_p
    print('completed '+p+'in %.fs' %time_p)

rt=time()-start_t
print('runtime: %.f' %rt)