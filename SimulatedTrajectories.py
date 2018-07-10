# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 13:27:47 2018

@author: Anjie
"""

import pandas as pd
import matplotlib.pyplot as plt
import json

import toolSpace as ts

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
f=ts.startFig('Simulated Trajectories','Time [s]','Position [m]',1)
for i in vehicles:
    veh=fzp.loc[fzp['no']==i]
    simplt=plt.plot(veh['sec'],veh['abspos'],'black')
#plt.legend()

#export plot
plt.savefig('Vehicle Trajectories of ' +segment+ '.pdf')

#%% misc. settings
plt.gca().set_ylim(420,540)
plt.gca().set_xlim(2320,2440)
#leadtraj=fzp.loc[fzp['no']==4633]
#leadtraj.to_csv('C:/Users/Anjie/OneDrive - University of Waterloo/MASc Research/EmissionModelling/LeadTrajectory.csv',sep=',',index=False)
