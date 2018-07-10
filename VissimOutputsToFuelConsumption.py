# -*- coding: utf-8 -*-
"""
Created on Wed May 30 17:14:07 2018

@author: Anjie
"""

import pandas as pd
import json, timeit

import toolSpace as ts

t0=timeit.default_timer()

#%% set up

# list fzp files
fzps=ts.listFiles('fzp')
lsas=ts.listFiles('lsa')

periods=['morning1','morning2','morning3','afternoon1','afternoon2'] # currently only morning periods are available
mov='thru' # use only thru movements as trajectory reconstruction is for thru movements, so that errors for actual and reconstructed trajectories can be compared

# intersections, directions
inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']
# get run specs
runSpecs=json.load(open('runSpecs.json'))
periodspecs=runSpecs['periods']
# get model specs
modelSpecs=json.load(open('modelSpecs.json'))

# vehicle masses
massdf=pd.read_excel('PhysicsTable.xlsx')
massmap=pd.Series(index=massdf['ID'],data=massdf['mass'].values)

# decoders
typemap=ts.getVehTypeMap()
typedecoder=pd.read_csv('MOVESdecoders/sourceType.csv')
typedecoder=pd.Series(typedecoder['sourceType'],index=typedecoder['ID'])

#%% get fuel consumption per cycle by approach and direction
fuel=pd.DataFrame()

dfperiod=[]
dfinter=[]
dfdirec=[]
dfcycle=[]
dfsim=[]
dftype=[]
dffuel=[]

s=1
# for each simulation
for s in range(len(fzps)):
    
    # import vissim data
    fzp=ts.import_fzp(fzps[s])
    lsa=ts.import_lsa(lsas[s])

    for p in periods:
        fzpp=fzp.loc[(fzp['sec']>=periodspecs[p]['start']) & (fzp['sec']<=periodspecs[p]['end'])]
        lsap=lsa.loc[(lsa['tsim']>=periodspecs[p]['start']) & (lsa['tsim']<=periodspecs[p]['end'])]
        
        for i in inters:
            
            for d in direcs:
                
                # Dunbar has no northbound and Hwy401 has no southbound available (end intersections)
                if not ((i==inters[0] and d==direcs[0]) or (i==inters[-1] and d==direcs[-1])):
                    t1=timeit.default_timer()
                    
                    # model specs
                    interspecs=modelSpecs[i]
                    allLinks=ts.getInterLinks(interspecs,d,mov) # get links for this intersection and approach direction
                    
                    sc=interspecs['SignalController'] # signal controller of intersection
                    sg=interspecs[d]['SignalGroup']['thru'] # signal group of interest
                    
                    #filterdata to intersection and direction of concern
                    fzp0=fzpp.loc[fzpp['link'].isin(allLinks)] # filter by links
                    lsa0=lsap.loc[(lsap['SC']==sc) & (lsap['SG']==sg)] # filter by signal controller and group
                    
                    ### sort data
                    
                    # find green phases
                    tg,tr=ts.findGreenPhases(lsa0)
                    
                    # filter by cycle 
                    cycles=range(len(tg))
                    for c in cycles:

                        fzpc=fzp0.loc[(fzp0['sec']>=tg[c]) & (fzp0['sec']<=tr[c])] # fzp for this intersection-approach and this cycle
                        types=fzpc['vehtype'].unique() # vehicle types in this cycle
                        for t in types:
                            
                            fzpct=fzpc.loc[fzpc['vehtype']==t]
                            mt=typemap[t] # moves type ID
                            mass=massmap[mt] # get mass from physics table
                            fuelct=ts.getFuel(fzpct,mass)
                            
                            dffuel.append(fuelct)
                            dfcycle.append(c+1)
                            dftype.append(mt)
                            dfsim.append(s)
                            dfperiod.append(p)
                            dfinter.append(i)
                            dfdirec.append(d)
                            
                    t1=timeit.default_timer()-t1
                            
                    print('Processed %i %s %s %s in %.1fs' %(s,p,i,d,t1))
                            
    s+=1

# populate dataframe and export
fuel['sim']=dfsim
fuel['period']=dfperiod
fuel['inter']=dfinter
fuel['direc']=dfdirec
fuel['cycle']=dfcycle
fuel['type']=dftype
fuel['fuel']=dffuel

fuel.to_csv('FuelEstimate-VissimOutput-'+str(mov)+'.csv',index=False)

runtime=timeit.default_timer()-t0
print('Runtime: %.2f' %runtime)