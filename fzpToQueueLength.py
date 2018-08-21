# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 19:39:09 2018

@author: Anjie
"""

import pandas as pd
import json
from time import time

import toolSpace as ts

t_start=time()

# capacity of a lane [veh/hr]
cap=2000

# get specs and properties

# intersections, directions
inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']

# get run specs
runSpecs=json.load(open('runSpecs.json'))
periods=runSpecs['periods']

# get model specs
modelSpecs=json.load(open('modelSpecs.json'))

fzps=ts.listFiles('fzp')
lsas=ts.listFiles('lsa')

linkIDlist=[]
qlengthlist=[]
congestionlist=[]
hourlist=[]
nvehicleslist=[]

# for each simulation
for s in range(len(fzps)):
    
    # import vissim data
    fzp=ts.import_fzp(fzps[s])
    lsa=ts.import_lsa(lsas[s])
    tg=lsa.loc[lsa['newState'].str.contains('green')] # filter signal timing to green times
    
    for p in periods:

        hour=periods[p]['hour']
        
        # filter by period
        fzp_p=ts.filterbyperiod(fzp,p)
        lsa_p=ts.filterbyperiod(lsa,p)    
            
        inter_index=1
        for i in inters: # for each intersection
            tgi=tg.loc[tg['SC']==modelSpecs[i]['SignalController']] # filter to signal controller

            direc_index=1
            for d in direcs: # for each direction (NB, SB)
    
                # Dunbar has no northbound and Hwy401 has no southbound available (end intersections)
                if not ((i==inters[0] and d==direcs[0]) or (i==inters[-1] and d==direcs[-1])):

                                        # segment model specs
                    segSpecs=modelSpecs[i][d]
                    
                    # links in the segment
                    links=segSpecs['links']['approach']['thru']
                                        
                    # link lengths
                    lengthMap=ts.getLengthMap(links)
                    segLength=sum(lengthMap)
                    
                    # filter signal timing data to segment
                    tgid=tgi.loc[lsa['SG']==segSpecs['StartingSignalGroup']] # filter signal timing to signal group
                    tgid=tgid.reset_index(drop=True) # reset index for filtering
                    tgidp=ts.filterbyperiod(tgid,p) # filter signal timing to greens starting this period
                    if p!='evening2':
                        tgidp=tgidp.append(tgid.loc[tgidp.index[-1]+1]) # add on the next green start after the last green start in this period to close the last bin
                    tgidp=tgidp['tsim'] # only the timing data is needed
                    if p=='evening2': # for last period, there is no next green after the period (simulation stopped)
                        tgidp=tgidp.append(pd.Series(fzp['sec'].iloc[-1]),ignore_index=True) # append last second of vehicle records
                    tgidp=tgidp.reset_index(drop=True) # reset index for iterating through cycles

                    
                    # vehicle records - get absolute position for vehicles
                    fzppid=fzp_p.loc[fzp_p['link'].isin(links)] # vehicle records for the intersection and approach
                    ## find absolute position of vehicle with end of last major intersection as reference
                    fzpx=0 #cumulative distance from reference
                    # main segment positions
                    for l in links:
                        fzppid.loc[fzppid['link']==l,'abspos']=fzppid['pos']+fzpx
                        fzpx +=lengthMap[l]

                    # iterate through cycles in this period
                    cycles=range(len(tgidp)-1) # the added on green start is not counted - just for closing the bin               
                    
                    for c in cycles:
                        
                        fzpc=fzppid.loc[(fzppid['sec']>=tgidp[c])&(fzppid['sec']<=tgidp[c+1])] # filter vehicle records to cycle
                        fzp_stopped=fzpc.loc[fzpc['speed']==0] # vehicle records for stopped vehicles
                        x_stopped=fzp_stopped['abspos'] # positions of stopped vehicles
                        qlength=0 # qlength is zero unless...
                        if not x_stopped.empty: # ...there was a queue
                            qlength=segLength-min(x_stopped) # then queue length is the segment length - the lowest (last) stopped position
                        nvehicles=len(fzpc['no'].unique())
                        Tcycle=tgidp[c+1]-tgidp[c]

                        if not fzpc['lane'].empty:
                            nlanes=max(fzpc['lane'])
                            congestion=nvehicles/nlanes/Tcycle*3600/2000
                        else:
                            nlanes=0
                            congestion=0

                        linkID=c+1+s*100+inter_index*10000+direc_index*1000 # linkID denoting intersection, direction and cycle
                        linkIDlist.append(linkID)
                        hourlist.append(hour)
                        qlengthlist.append(qlength)
                        congestionlist.append(congestion)
                        nvehicleslist.append(nvehicles)

                        
                direc_index+=1
            inter_index+=1

output=pd.DataFrame({'link':linkIDlist,'hour':hourlist,'qlength':qlengthlist,'congestion':congestionlist,'nvehicles':nvehicleslist})
output.to_csv('metrics.csv') # save as csv

print('runtime: %.1fs' %(time()-t_start))