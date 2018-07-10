# -*- coding: utf-8 -*-
"""
Created on Wed May  2 11:49:06 2018

@author: Anjie
"""

import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np

from toolSpace import startFig, getLengthMap

#manual inputs
lanes=1

plt.style.use('seaborn-whitegrid')

#%%import data

#import simulation data
visOutputPath='C:/Users/anjie/OneDrive - University of Waterloo/MASc Research/VISSIMmodel/WholeCorridor/SimOutputs-wDataCollectionPoints/Hespeler--Dunbar-Hwy401_AllVehTypes_wDataCollectionPoints_001'
#import simulated detector (simDetector) data
drow=6
hrow=157 #row number for headers where detector data starts
mer=pd.read_csv(visOutputPath+'.mer',header=hrow,sep=';')
mer=mer.rename(columns=lambda x: x.strip()) #strip whitespace from headers
mer=mer.astype({'Measurem.':int,'t(Entry)':float,'t(Exit)':float,'VehNo':int,'Vehicle type':int})
merDetail=pd.read_csv(visOutputPath+'.mer',nrows=hrow-drow-1,skiprows=9,header=None,delim_whitespace=True)
merDetail[3]=merDetail[3].replace({':':''}, regex=True)
merDetail=merDetail.astype({3:int})

#import signal phasing
lsa=pd.read_csv(visOutputPath+'.lsa',index_col=False,skiprows=126,sep=';',names=['tsim','tcycle','SC','SG','newState','lastDuration','type','SGcause'])
lsa=lsa.astype({'tsim':float,'newState':str})

#import travel times
rsr=pd.read_csv(visOutputPath+'.rsr',sep=';',index_col=False,header=5)
rsr=rsr.rename(columns=lambda x: x.strip())

#%% get simulation specs

segment0='DunbarBishop'#segment of interest
segment1='BishopSheldon'

#load simulation and model specifications
simSpecs=json.load(open('simSpecs.json'))
segment0=simSpecs['segments'][segment0]
segment1=simSpecs['segments'][segment1]
seglinks0=segment0['links'] #IDs of links in Vissim model categorized by appproach and movement
wt=simSpecs['warmUpTime'] #warm-up time: time at beginning of simulation to exclude [s]
pt=simSpecs['durationOfInterest'] #period to plot or analyse
mainthrulinks=seglinks0['main']['thru']
#mainLTlinks=seglinks0['main']['LT']

#simulation detector numbers
simDet={}
simDet['adv']=segment0['detectors']['advance']
simDet['stop']=segment0['detectors']['stopbar']
simDet['end']=segment1['detectors']['upstream']

x={} #positions of detectors
for sd in simDet:
    x[sd]=int(merDetail.loc[merDetail[3]==simDet[sd][0],9].iloc[0])
x['end']=x['stop']+40 #add on intersection width
lengthMap=getLengthMap(mainthrulinks)
print(pd.Series(mainthrulinks))
posOffset=sum(lengthMap.loc[pd.Series(mainthrulinks).iloc[0:-1]])
for sd in x:
    x[sd]=x[sd]+posOffset

#TODO: don't hardcode intersection width
#%% filter data

#filter simDetector data

sd_all=simDet['adv']+simDet['stop']+simDet['end'] #all detectors
if lanes==1:
    sd_all=[simDet['adv'][0],simDet['stop'][0],simDet['end'][0]] #detectors of one lane
mer=mer.loc[mer['Measurem.'].isin(sd_all)] #by simDetectors of interest
merEntry=mer.loc[(mer['t(Entry)']>wt)&(mer['t(Entry)']<wt+pt)] #by detector on times
merExit=mer.loc[(mer['t(Exit)']>wt)&(mer['t(Exit)']<wt+pt)] #by detector off times

#filter signal phasing to time and Signal Group of interest
lsa=lsa.loc[lsa['SC']==segment0['SC']]
lsa=lsa.loc[lsa['SG']==segment0['SG']['thru']]
lsa=lsa.loc[lsa['tsim']>wt]
lsa=lsa.loc[lsa['tsim']<wt+pt]
lsa=lsa.reset_index(drop=True)

#filter travel time data
rsr=rsr.loc[rsr['Time']>=wt]
rsr=rsr.loc[rsr['Time']<=wt+pt]
rsr=rsr.loc[rsr['No.']==10]

#%% compose inputs 
#compose simDetector data for plotting
vehicles=merEntry['VehNo'].unique()
t1=[]
t2=[]
vehn=[]
dpos=[]
for n in vehicles:
    merEntry_n=merEntry.loc[merEntry['VehNo']==n]
    merExit_n=merExit.loc[merExit['VehNo']==n]
    for sd in simDet:
        for d in simDet[sd]:
            merEntry_nd=merEntry_n.loc[merEntry['Measurem.']==d]
            merExit_nd=merExit_n.loc[merExit['Measurem.']==d]
            if not (merEntry_nd['t(Entry)'].empty or merExit_nd['t(Exit)'].empty):
                vehn.append(n)
                dpos.append(x[sd])
                t1.append(merEntry_nd['t(Entry)'].iloc[0])
                t2.append(merExit_nd['t(Exit)'].iloc[0])
detTraj=pd.DataFrame({'vehn':vehn,'x':dpos,'t1':t1,'t2':t2})

#compose signal data for plotting
tred=lsa.loc[lsa['newState'].str.contains('red'),'tsim']
redindex=tred.index.values
tgreen=lsa.loc[lsa['newState'].str.contains('green'),'tsim']

#%% plot inputs

#plot simDetector data
ps=json.load(open('plotSpecs.json'))
f=startFig('Trajectory Reconstruction','Time [s]','Position [m]',1)

#plot simdetection
for n in range(len(vehn)):
    pltx=[t1[n],t2[n]]
    plty=[dpos[n]]*2
    detplt=plt.plot(pltx,plty,linewidth='4',color='orange')
#connect detections by vehicle
for n in vehicles:
    detTrajn=detTraj.loc[detTraj['vehn']==n]
    detTrajn_adv=detTrajn.loc[detTrajn['x']==x['adv']]
    detTrajn_stop=detTrajn.loc[detTrajn['x']==x['stop']]
    detTrajn_end=detTrajn.loc[detTrajn['x']==x['end']]
    if not (detTrajn_stop.empty or detTrajn_end.empty): #plot thru intersection interpolations
        pltx=[detTrajn_stop['t2'].iloc[0],detTrajn_end['t2'].iloc[0]]
        plty=[x['stop'],x['end']]
#        plt.plot(pltx,plty,alpha=.3)
    if not (detTrajn_stop.empty or detTrajn_adv.empty): #plot approach interpolations
        pltx=[detTrajn_adv['t1'].iloc[0],detTrajn_stop['t1'].iloc[0]]
        plty=[x['adv'],x['stop']]
#        plt.plot(pltx,plty,alpha=.3)

#plot signal phasing
for i in redindex:
    if i+1 < len(lsa):
        sigplt=plt.plot([lsa.loc[i,'tsim'],lsa.loc[i+1,'tsim']],[x['stop']+2]*2,alpha=.9,linewidth=6,color='r')
        
#%% compose inputs for lead vehicle trajectories

#find lead vehicles based on wait time at stopline
detTrajStop=detTraj.loc[detTraj['x']==x['stop']]
detTrajStop['wait']=detTrajStop['t2']-detTrajStop['t1']
maxwaits=detTrajStop['wait'].nlargest(lanes)
leadvehs=detTrajStop.loc[detTrajStop['wait'].isin(maxwaits)]
leadvehns=leadvehs['vehn'].tolist()

detTrajAdv=detTraj.loc[detTraj['x']==x['adv']]
tadv=detTrajAdv.loc[detTrajAdv['vehn'].isin(leadvehns)]['t1'] #time when lead vehicle is at advance detector

detTrajEnd=detTraj.loc[detTraj['x']==x['end']]

#inputs for deceleration/acceleration curve construction
ttcruise=np.percentile(rsr['Trav'],0) #nth percentile travel time to find cruise speed
#TODO: length of segment
v0=521.75/ttcruise #starting speed
dt=.2 #time step [s]
tdelta=9999
p=.1 #parameter adjustment step for search

#%% trajectory of lead vehicles
leadTraj={}
for i in range(len(leadvehns)): #for each lead vehicle
    leadTraj[i]={}
    
    #deceleration curve
    #original values from literature
    k3o=0.005 
    k4o=0.154
    k5o=0.493 
    #start with original values
    k3=k3o
    k4=k4o
    k5=k5o
    while tdelta>.2:#time difference between advance detection and extrapolated cruise trajectory at advance
        v=[v0]
        dc=[0]
        xdec=[0]
        tdec=[0]
        while v[-1]>0:
            dc.append(-k3*v[-1]**2+k4*v[-1]+k5) #deceleration curve from literature
            v.append(v[-1]-(dc[-1]+dc[-2])/2*dt)
            xdec.append(xdec[-1]+v[-1]*dt+(dc[-1]+dc[-2])/4*dt**2)
            tdec.append(tdec[-1]+dt)
            
        #shift deceleration trajectory to match end with beginning of stopped trajectory
        tdeci=tdec+leadvehs['t1'].iloc[i]-tdec[-1]
        xdeci=xdec+leadvehs['x'].iloc[i]-xdec[-1]
        tcruise=(xdeci[0]-x['adv'])/v[0]
        textra=tdeci[0]-tcruise #extrapolated time of reaching advance detector
        tdelta=tadv.iloc[i]-textra
        if xdeci[0]<x['adv']:
            tdelta=tadv.iloc[i]-tdeci[abs(xdeci-x['adv']).argmin()]
        if tdelta>0:
            k3=k3-k3o*p 
            k4=k4+k4o*p
            k5=k5+k5o*p
        else:
            k3=k3+k3o*p 
            k4=k4-k4o*p
            k5=k5-k5o*p
            
    #stopped curve
    tStopEnd=detTrajStop.loc[detTrajStop['vehn']==leadvehns[i],'t2'].iloc[0]
    
    #acceleration curve
    v=[0]    
    tdelta=9999
    a=[0]
    tdwell=1 #time vehicle dwells in the detector position while accelerating
    #calibrated parameters from literature
    beta0=3.369
    beta1=0.072
    
    xacc=[0]
    tacc=[0]
    while v[-1]<v0:
        a.append(beta0-beta1*v[-1])
        v.append(v[-1]+(a[-1]+a[-2])/2*dt)
        xacc.append(xacc[-1]+v[-1]*dt+(dc[-1]+dc[-2])/4*dt**2)
        tacc.append(tacc[-1]+dt)
    tacc+=tStopEnd-tdwell
    xacc=[xx+x['stop'] for xx in xacc]
    leadTraj[i]['t']=[textra]+list(tdeci)+list(tacc)
    leadTraj[i]['x']=[x['adv']]+list(xdeci)+list(xacc) 
    
    #put together lead vehicle trajectory

    trajplt=plt.plot(leadTraj[i]['t'],leadTraj[i]['x'],'b--')
    
#%% construct follow wait trajectories

traj={0:{'x':leadTraj[0]['x'],'t':leadTraj[0]['t']}}
#get car-following parameters time interval and effective vehicle length
leff=5.5 #effective length [m]
hDis=2 #maximum headway for queue discharge
#TODO: effective length depends on vehicle type

folvehs=detTraj.loc[~detTraj['vehn'].isin(leadvehns)]
folvehns=folvehs.loc[folvehs['x']==x['adv'],'vehn'].unique()
folvehns2=folvehs.loc[folvehs['x']==x['end'],'vehn'].unique()
folvehns=np.intersect1d(folvehns,folvehns2)
detTrajFol=detTraj.loc[detTraj['vehn'].isin(folvehs['vehn'])]
vehns=leadvehns+list(folvehns)
tadvfol=detTrajAdv.loc[detTrajAdv['vehn'].isin(folvehns)]['t1'].tolist()
tstopfol=detTrajStop.loc[detTrajStop['vehn'].isin(folvehns)]['t1'].tolist()
tendfol=detTrajEnd.loc[detTrajEnd['vehn'].isin(folvehns)]['t1'].tolist()
tendfol2=detTrajEnd.loc[detTrajEnd['vehn'].isin(folvehns)]['t2'].tolist()

tstop=detTrajStop.loc[detTrajStop['vehn']==vehns[0],'t1'].iloc[0]
tgo=detTrajStop.loc[detTrajStop['vehn']==vehns[0],'t2'].iloc[0]
xstop=x['stop']

labels=[]

for j in range(len(folvehns)):
    traj[j+1]={}
    vapproach=(x['stop']-x['adv'])/(tstopfol[j]-tadvfol[j])
    vthru=(x['end']-x['stop'])/(tendfol[j]-tstopfol[j])
    h=tstopfol[j]-tstopfol[j-1]
    
    if vapproach<.7*vthru: #check if vehicle stopped in queue and followed closely
        #find T (time interval) for Newell's model
        detTrajEnd=detTraj.loc[detTraj['x']==x['end']]
        t2=detTrajEnd.loc[detTrajEnd['vehn'].isin(folvehns),'t2'].iloc[j]
        xlead=x['end']+leff
        index=abs(np.array(traj[j]['x'])-xlead).argmin()
        t1=traj[j]['t'][index]
    
        if xlead>traj[j]['x'][-1]:
            t1=traj[j]['t'][-1]+(x['end']+leff-traj[j]['x'][-1])/v0
        tlag=t2-t1
        
        tgo+=tlag
        xstop-=leff
        
        xdeci=pd.Series(xdec).add(xstop-xdec[-1])
    
        if xdeci[0]>x['adv']:
            tdeci=pd.Series(tdec).add(tadvfol[j]+(xdeci[0]-x['adv'])/v0)
            tdeci=pd.Series(tadvfol[j]).append(tdeci)
            xdeci=pd.Series(x['adv']).append(xdeci)
        else:
            tdelta=tdec[abs(np.array(xdeci-x['adv'])).argmin()]-tdec[0]
            tdeci=pd.Series(tdec).add(tadvfol[j]-tdelta)
        
        #acceleration curve
        v=[0]    
        tdelta=9999
        a=[0]
        tdwell=1 #time vehicle dwells in the detector position while accelerating
        #calibrated parameters from literature
        beta0=3.369
        beta1=0.072
        
        xacc=[0]
        tacc=[0]
        while v[-1]<v0:
            a.append(beta0-beta1*v[-1])
            v.append(v[-1]+(a[-1]+a[-2])/2*dt)
            xacc.append(xacc[-1]+v[-1]*dt+(a[-1]+dc[-2])/4*dt**2)
            tacc.append(tacc[-1]+dt)
        tacc+=tgo-tdwell
        xacc=[xx+xstop for xx in xacc]
        
        #put together vehicle trajectory
        traj[j+1]['t']=list(tdeci)+list(tacc)
        traj[j+1]['x']=list(xdeci)+list(xacc)
        
        labels.append(str(folvehns[j])+' stopped')
    
    elif h<hDis: #didn't stop but followed closely
        t0=tadvfol[j]
        traj[j+1]={'t':[t0],'x':[x['adv']]}
        for i in range(1,len(traj[j-1]['t'])):
            #build trajectory
            tfol=traj[j]['t'][i]+tlag
            xfol=traj[j]['x'][i]-leff
            traj[j+1]['x'].append(xfol)
            tfolcruise=(traj[j+1]['x'][i]-traj[j+1]['x'][i-1])/v0+traj[j+1]['t'][i-1]
            if tfolcruise>tfol:
                tfol=tfolcruise
            traj[j+1]['t'].append(tfol)
            
        labels.append(str(folvehns[j])+' following')
            
    else: #free flow vehicles
        traj[j+1]['t']=[tadvfol[j],tstopfol[j],tendfol[j]]
        traj[j+1]['x']=[x['adv'],x['stop'],x['end']]
        
        labels.append(str(folvehns[j])+' free flow')

for j in range(len(folvehns)):
    if traj[j+1]['x'][-1]<x['end']:
        traj[j+1]['x']=traj[j+1]['x']+[x['end']]
        traj[j+1]['t']=traj[j+1]['t']+[tendfol2[j]]
    plt.plot(traj[j+1]['t'],traj[j+1]['x'],'b--')

#%% Simulated Trajectories
    
segment='DunbarBishop' #segment of interest

#load simulation and model specifications
simSpecs=json.load(open('simSpecs.json'))

seglinks=simSpecs['segments'][segment]['links'] #IDs of links in Vissim model categorized by appproach and movement
ltLanes=simSpecs['LTlaneConnections'] #links and positions where left-turn lanes connect
wt=simSpecs['warmUpTime'] #warm-up time: time at beginning of simulation to exclude [s]
pt=simSpecs['durationOfInterest'] #period to plot or analyse

#file path of Vissim outputs
visOutputPath='C:/Users/anjie/OneDrive - University of Waterloo/MASc Research/VISSIMmodel/WholeCorridor/SimOutputs-wDataCollectionPoints/Hespeler--Dunbar-Hwy401_AllVehTypes_wDataCollectionPoints_001'
fzpPath=visOutputPath +'.fzp'

#Load vehicle records to dataframe
fzp=pd.read_csv(fzpPath,sep=';',skiprows=range(0,22),
                names=['sec','no','link','lane','dist','pos','speed','accel','routeDec','route','vehType','sim'])
#                dtype={'sec':float,'no':int,'link':int,'lane':int,'dist':float,'pos':float,'speed':float,
#                       'accel':float,'routeDec':int,'route':int,'vehType':int,'sim':int})
print('Vehicle records imported.')

#filter vehicle records by links of interest

links0=[]
for i in seglinks:
    for j in seglinks[i]:
        links0.extend(seglinks[i][j])
fzp=fzp.loc[fzp['link'].isin(links0)]
print('Filtered to links of interest.')

fzp=fzp.loc[fzp['lane']==3] #filter to lane 
#fzp=fzp.loc[fzp['vehType']==200] #filter by vehicle type

#remove vehicle records from warm-up period and period after plot duration
fzp=fzp.loc[fzp['sec']>wt]
fzp=fzp.loc[fzp['sec']<wt+pt]

#import link lengths
linkProperties=pd.read_excel('LinkProperties.xlsx',header=1)
#create table of links and their lengths
linklengths=[]
for l in links0:
    linklengths.append(linkProperties.loc[linkProperties['$LINK:NO']==l,'LENGTH2D'].values[0])
lengthMap=pd.Series(linklengths,index=links0)

##find absolute position of vehicle with end of last major intersection as reference
x=0 #cumulative distance from reference
#main segment
for l in seglinks['main']['thru']:
    fzp.loc[fzp['link']==l,'abspos']=fzp['pos']+x
    x +=lengthMap[l]
#output intersection
for m in seglinks['downstream']:
    l=seglinks['downstream'][m][0]
    fzp.loc[fzp['link']==l,'abspos']=fzp['pos']+x
#left turn storage lane
ltLanes0=ltLanes[segment]
for l in seglinks['main']['LT']:
    delta=lengthMap[ltLanes0['fromLane']]-ltLanes0['at']
    fzp.loc[fzp['link']==l,'abspos']=fzp['pos']+x-delta
    x +=lengthMap[l]
#input intersection
for m in seglinks['upstream']:
    l=seglinks['upstream'][m][0]
    fzp.loc[fzp['link']==l,'abspos']=fzp['pos']-lengthMap[l]
print('Positions calculated.')

#fzp=fzp.loc[fzp['sec']<wt+120] #filter to first cycle
vehicles=fzp['no'].unique() #all vehicle numbers that appear
for i in vehicles:
    veh=fzp.loc[fzp['no']==i]
    simplt=plt.plot(veh['sec'],veh['abspos'],'black')
#%%

#plt.legend((sigplt[0],detplt[0]),('Red Signal Phase','Simulated Detector Data'))
plt.legend((sigplt[0],detplt[0],trajplt[0],simplt[0]),('Red Signal Phase','Simulated Detector Data','Reconstructed Trajectory','Simulated Trajectories'))

#%% misc. settings
plt.gca().set_ylim(420,540)
plt.gca().set_xlim(2320,2440)
