# -*- coding: utf-8 -*-
"""
Created on Mon May 28 18:57:30 2018

@author: Anjie
"""

import pandas as pd
import numpy as np
import os, json, timeit

import toolSpace as ts

t_start=timeit.default_timer()

#%% manual inputs
mov = 'thru' # movements of interest: {'thru','all'}

#%% get specs and properties

# intersections, directions, movements
inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']

# get run specs
runSpecs=json.load(open('runSpecs.json'))
periods=runSpecs['periods']

# get model specs
modelSpecs=json.load(open('modelSpecs.json'))

#%%import vissim data and format to join with movesoutput for tableau

cols=['link','type','Vehicles','avgspd']
cols_rename={'avgspd':'Average Speed [km/h]'}

# list files
fzps=ts.listFiles('fzp')
lsas=ts.listFiles('lsa')

fzp_ft={}

# get data to match links for opmodes

for s in range(len(fzps)):
    # import vissim data
    fzp=ts.import_fzp(fzps[s])
    
    for p in periods:
        
        fzp_ft[p]=pd.DataFrame(columns=cols) # start template
        fzp_p=ts.filterbyperiod(fzp,p) # filter by period
        
        inter_index=1
        for i in inters: # for each intersection
            direc_index=1
            for d in direcs: # for each direction (NB, SB)
    
                # Dunbar has no northbound and Hwy401 has no southbound available (end intersections)
                if not ((i==inters[0] and d==direcs[0]) or (i==inters[-1] and d==direcs[-1])):
                    
                    # model specs
                    interspecs=modelSpecs[i]
                    links=interspecs[d]['links']
                    
                    # list links for intersection and approach
                    allLinks=[]
                    if mov=='all':
                        for section in links:
                            for movement in links[section]:
                                allLinks.extend(links[section][movement])
                    else:
                        for section in links:
                            allLinks.extend(links[section]['thru'])
                            
                    #filter vissim data to intersection and direction of concern
                    fzp=fzp.loc[fzp['link'].isin(allLinks)] # filter by links
                    lsa=lsa.loc[lsa['SC']==sc] # filter by signal controller
                    lsa=lsa.loc[lsa['SG']==sg] # filter by signal group
                    
                    ### sort data
                    
                    # find green phases
                    tg=lsa.loc[lsa['newState'].str.contains('green'),'tsim'] # green start times
                    lsa=lsa.loc[lsa['tsim']>=tg.iloc[0]] # filter to after first green start
                    tr=lsa.loc[lsa['newState'].str.contains('red'),'tsim'] # red start times
                    tg=tg.loc[tg<=tr.iloc[-1]] # remove green starts without an end
                    tg=tg.reset_index(drop=True)
                    tr=tr.reset_index(drop=True)
                    
                    # list vehicle numbers for each green phase
                    
                    for c in range(len(tg)): # for each cycle
                        fzpc=fzp.loc[(fzp['sec']>=tg[c]) & (fzp['sec']<=tr[c])] # fzp for this intersection-approach and this cycle
                        vehc=fzpc['no'].unique() # vehicle numbers in this cycle
                        
                        types=fzpc['vehtype'].unique()
                        typecounts=fzpc['vehtype'].value_counts()
                        for t in types:
                            
                    