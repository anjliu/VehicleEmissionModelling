# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 14:01:26 2018

@author: Anjie
"""

import json

import toolSpace as ts

# conversion variables
m2k=1.60934 # km per mile
kph2mps=1/3.6 # m/s per km/h
mph2mps=m2k*kph2mps
grade=0

# get run specs
runSpecs=json.load(open('runSpecs.json'))
periods=runSpecs['periods']

# get model specs
modelSpecs=json.load(open('modelSpecs.json'))

# specs for MOVES inputs
county=runSpecs['county']
zone=runSpecs['zone']
roadType=runSpecs['roadType']
hourDayID=runSpecs['hourDayID']

opmFormat=runSpecs['opModeFormat']
opmModeSetProcesses=opmFormat['processes']['modeSet']
opmAllRunningProcesses=opmFormat['processes']['allRunning']
opmModeSet=opmFormat['opModeSets']['modeSet']
nmodes=len(opmModeSet)
opmAllRunningMode=opmFormat['opModeSets']['allRunning'][0]
spdbins=opmFormat['speedbins']
for i in spdbins: # convert speeds from miles/hr to km/hr
    i=i*m2k
    
vspbins={}
for i in opmFormat['vspbins']:  
    vspbins[int(i)]=opmFormat['vspbins'][i]
    


# input combinations
types=[21,42,52,61]

# lists for MOVES inputs

ls_linkID=[]
ls_linkVolume=[]

lst_linkID=[]
lst_sourceTypeID=[]
lst_fraction=[]

opm_type=[]
opm_link=[]
opm_process=[]
opm_mode=[]
opm_fraction=[]

# create MOVES inputs and links for each combination of vehicle type, operating mode, pollutant process

for t in types: # for all vehicle types
    for m in opmModeSet: # for running operating operating modes categorized by speed and VSP
        linkid=t*100+m
        
        ls_linkID.append(linkid)
        ls_linkVolume.append(1)
        
        lst_linkID.append(linkid)
        lst_sourceTypeID.append(t)
        lst_fraction.append(1)
        
        for pr in opmModeSetProcesses: # for pollutant processes for categorized operating modes
            
            opm_type.append(t)
            opm_link.append(linkid)
            opm_process.append(pr)
            opm_mode.append(m)
            opm_fraction.append(1)
        
        for pr in opmAllRunningProcesses: # for pollutant processes for the all running operating mode
        
            opm_type.append(t)
            opm_link.append(linkid)
            opm_process.append(pr)
            opm_mode.append(opmAllRunningMode)
            opm_fraction.append(1)

# Set up for result saving - create dataframes for MOVES inputs
ls=ts.makeLinksTemplate()
lst=ts.makeLinkSourceTypeTemplate()
opm=ts.makeOpModeTemplate()

# links
ls['linkID']=ls_linkID
lsSize=len(ls)
ls['linkLength']=[1]*lsSize
ls['linkVolume']=ls_linkVolume
ls['countyID']=[county]*lsSize
ls['zoneID']=[zone]*lsSize
ls['roadTypeID']=[roadType]*lsSize
ls['linkAvgSpeed']=[0]*lsSize

# link source types
lst['linkID']=lst_linkID
lst['sourceTypeID']=lst_sourceTypeID
lst['sourceTypeHourFraction']=lst_fraction

# operating mode distributions
opmSize=len(opm_link)
opm['sourceTypeID']=opm_type
opm['hourDayID']=[hourDayID]*opmSize
opm['linkID']=opm_link
opm['polProcessID']=opm_process
opm['opModeID']=opm_mode
opm['opModeFraction']=opm_fraction

savepath='MOVESinputs/rates'
# export to csv
ls.to_excel(savepath+'/links.xlsx',index=False)
lst.to_excel(savepath+'/linkSourceTypes.xlsx',index=False)
opm.to_excel(savepath+'/OpModeDistribution.xlsx',index=False)
