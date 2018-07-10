# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 14:50:59 2018

@author: Anjie
"""

import pandas as pd
import json, timeit
import numpy as np

import toolSpace as ts

t0=timeit.default_timer()

#%% set up

# conversion variables
m2k=1.60934 # km per mile
kph2mps=1/3.6 # m/s per km/h
mph2mps=m2k*kph2mps
grade=0

spdcut=25*m2k

# list fzp files
fzps=ts.listFiles('fzp')
lsas=ts.listFiles('lsa')
rsrs=ts.listFiles('rsr')

periods=['morning1','morning2','morning3','afternoon1','afternoon2'] # currently available
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
massmap=pd.Series(massdf['mass'].values,index=massdf['ID'])

# decoders
typemap=ts.getVehTypeMap()
typedecoder=pd.read_csv('MOVESdecoders/sourceType.csv')
alltypes=typedecoder['ID']
typedecoder=pd.Series(typedecoder['sourceType'],index=alltypes)

# load VSP parameters. Source: Population and Activity of On-Road Vehicles in MOVES2014, page 103
vspParams=pd.read_excel('PhysicsTable.xlsx',sheet_name='Revised',dtype={'sourceTypeID':int,'cA':float,'cB':float,'cC':float,'mass':float})
params={}
for i in alltypes:
    params[i]=ts.getVSPparams(i,vspParams)

#%% get fuel consumption per cycle by approach and direction
vspDF=pd.DataFrame()

dfperiod=[]
dfinter=[]
dfdirec=[]
dfcycle=[]
dfsim=[]
dftype=[]
dfvspfast=[]
dfvspslow=[]
dfvspmass=[]
dfvehslow=[]
dfvehfast=[]
dfavgspeed=[]
dfavgaccel=[]
dfavgTT=[]

def periodfilter(df,timeheader):
    dfout=df.loc[(df[timeheader]>=periodspecs[p]['start']) & (df[timeheader]<=periodspecs[p]['end'])]
    return dfout

s=1
# for each simulation
for s in range(len(fzps)):
    
    # import vissim data
    fzp=ts.import_fzp(fzps[s])
    lsa=ts.import_lsa(lsas[s])
    rsr=ts.import_rsr(rsrs[s])

    for p in periods:
        fzpp=periodfilter(fzp,'sec')
        lsap=periodfilter(lsa,'tsim')
        rsrp=periodfilter(rsr,'Time')
        
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
                    
                    TTmeas=modelSpecs[i][d]['TravelTimeMeas'] # travel time measurer number
                    
                    #filterdata to intersection and direction of concern
                    fzp0=fzpp.loc[fzpp['link'].isin(allLinks)] # filter by links
                    lsa0=lsap.loc[(lsap['SC']==sc) & (lsap['SG']==sg)] # filter by signal controller and group
                    rsr0=rsrp.loc[rsrp['No.']==TTmeas]
                    
                    ### sort data
                    
                    # find green phases
                    tg,tr=ts.findGreenPhases(lsa0)
                    
                    # filter by cycle 
                    cycles=range(len(tg))
                    for c in cycles:

                        fzpc=fzp0.loc[(fzp0['sec']>=tg[c]) & (fzp0['sec']<=tr[c])] # fzp for this intersection-approach and this cycle
                        types=fzpc['vehtype'].unique() # vehicle types in this cycle
                        
                        rsrc=rsr0.loc[(rsr0['Time']>=tg[c]) & (rsr0['Time']<=tr[c])]
                        
                        for t in types:
                            
                            # filter to vehicle type
                            fzpct=fzpc.loc[fzpc['vehtype']==t]
                            rsrct=rsrc.loc[rsrc['VehType']==t]
                            
                            fzpslow=fzpct.loc[fzpct['speed']<spdcut]
                            fzpfast=fzpct.loc[fzpct['speed']>=spdcut]
                            
                            mt=typemap[t] # moves type ID
                            mass=massmap[mt] # get mass from physics table
                            vehslow=len(fzpslow['no'].unique()) # number of vehicles
                            vehfast=len(fzpfast['no'].unique()) # number of vehicles
                            
                            # calculate VSP and get other parameters
                            vspfast=ts.getVSP(fzpfast['speed']*kph2mps,fzpfast['accel'],grade,params[mt],mt)
                            vspslow=ts.getVSP(fzpslow['speed']*kph2mps,fzpslow['accel'],grade,params[mt],mt)
                            totalvspfast=sum(vspfast)
                            totalvspslow=sum(vspslow)
                            totalvsp=sum([totalvspfast,totalvspslow])
                            avgspeed=np.mean(fzpct['speed'])
                            avgaccel=np.mean(fzpct['accel'])
                            avgTT=np.mean(rsrct['Trav'])
                            
                            dfvspfast.append(totalvspfast)
                            dfvspslow.append(totalvspslow)
                            dfcycle.append(c+1)
                            dftype.append(mt)
                            dfsim.append(s)
                            dfperiod.append(p)
                            dfinter.append(i)
                            dfdirec.append(d)
                            dfvspmass.append(totalvsp*mass)
                            dfvehslow.append(vehslow)
                            dfvehfast.append(vehfast)
                            dfavgspeed.append(avgspeed)
                            dfavgaccel.append(avgaccel)
                            dfavgTT.append(avgTT)
                            
                    t1=timeit.default_timer()-t1
                            
                    print('Processed %i %s %s %s in %.1fs' %(s,p,i,d,t1))
                            
    s+=1

# populate dataframe and export
vspDF['sim']=dfsim
vspDF['period']=dfperiod
vspDF['inter']=dfinter
vspDF['direc']=dfdirec
vspDF['cycle']=dfcycle
vspDF['type']=dftype
vspDF['vspfast']=dfvspfast
vspDF['vspslow']=dfvspslow
vspDF['vspmass']=dfvspmass
vspDF['vehicleslow']=dfvehslow
vspDF['vehiclefast']=dfvehfast
vspDF['avgspeed']=dfavgspeed
vspDF['avgaccel']=dfavgaccel
vspDF['avgTT']=dfavgTT

vspDF.to_csv('VSP-VissimOutput-'+str(mov)+'.csv',index=False)

runtime=timeit.default_timer()-t0
print('Runtime: %.2f' %runtime)