# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:42:22 2018

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

# intersections, directions
inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']

# conversion variables
m2k=1.60934 # km per mile
kph2mps=1/3.6 # m/s per km/h
mph2mps=m2k*kph2mps
grade=0

# get run specs
runSpecs=json.load(open('runSpecs.json'))
periods=runSpecs['periods']

# get model specs
modelSpecs=json.load(open('modelSpecs-MajorOnly.json')) # approach segments to major intersections excluding minor intersections

# specs for MOVES
county=runSpecs['county']
zone=runSpecs['zone']
roadType=runSpecs['roadType']
hourDayID=runSpecs['hourDayID']

opmFormat=runSpecs['opModeFormat']
opmModeSetProcesses=opmFormat['processes']['modeSet']
opmAllRunningProcesses=opmFormat['processes']['allRunning']
opmModeSet=opmFormat['opModeSets']['modeSet']
modes=len(opmModeSet)
opmAllRunning=opmFormat['opModeSets']['allRunning'][0]
spdbins=opmFormat['speedbins']
for i in spdbins: # convert speeds from miles/hr to km/hr
    i=i*m2k
    
vspbins={}
for i in opmFormat['vspbins']:  
    vspbins[int(i)]=opmFormat['vspbins'][i]
    
brake_threshold=-2*mph2mps

# load VSP parameters. Source: Population and Activity of On-Road Vehicles in MOVES2014, page 103
vspParams=pd.read_excel('PhysicsTable.xlsx',sheet_name='Revised',dtype={'sourceTypeID':int,'cA':float,'cB':float,'cC':float,'mass':float})

# vehicle type ID mapping from Vissim to MOVES
typemap=ts.getVehTypeMap()

def findCycleVehicles(fzp,tstart,tend,links):
    fzp=fzp.loc[fzp['link'].isin(links)] # filter to links used
    fzp=fzp.loc[(fzp['sec']>=tstart) & (fzp['sec']<=tend)]
    vehcycle=fzp['no'].unique()
    return(vehcycle)

#%% function getting operating mode distribution 

def getOpMode(inter,direc,movements,fzp,lsa):

    # model specs
    interspecs=modelSpecs[inter]
    
    allLinks=ts.getSegmentLinks(interspecs,direc,movements)
    intLinks=ts.getIntLinks(interspecs,direc,movements)
    
    sc=interspecs['SignalController'] # signal controller of intersection
    sg=interspecs[direc]['StartingSignalGroup'] # signal group of interest
    
    # get link lengths
    lengthMap=ts.getLengthMap(allLinks)
    lengthtotal=sum(lengthMap)
    
    # Set up for result saving - create dataframes for MOVES inputs
    ls=ts.makeLinksTemplate()
    lst=ts.makeLinkSourceTypeTemplate()
    opm=ts.makeOpModeTemplate()
 
    # filter vissim data to intersection and direction of concern
    fzp=fzp.loc[fzp['link'].isin(allLinks)] # filter by links
    lsa=lsa.loc[(lsa['SC']==sc) & (lsa['SG']==sg)] # filter by signal controller and group

    ### sort data
    
    # find green phases
    tg=lsa.loc[lsa['newState'].str.contains('green'),'tsim']
    tg=tg.reset_index(drop=True)
    
    # list vehicle numbers and fzp for each green phase
    vehpercycle={}
    fzppercycle={}
    cycles=range(len(tg))
    tg=tg.append(pd.Series([99999999]),ignore_index=True) # add large value to end to close the last bin
    
    for c in cycles:
        vehpercycle[c]=findCycleVehicles(fzp,tg[c],tg[c+1],intLinks) # vehicle numbers in this cycle
        fzppercycle[c]=fzp.loc[fzp['no'].isin(vehpercycle[c])] # fzp for vehicles crossing in this cycle
        
    ### calculate traffic flow characteristics
    
    ls_linkID=[]
    ls_linkVolume=[]
    ls_linkAvgSpeed=[]
    
    lst_linkID=[]
    lst_sourceTypeID=[]
    lst_fraction=[]
    
    opm_type=[]
    opm_link=[]
    opm_process=[]
    opm_mode=[]
    opm_fraction=[]
    
    opmSeconds=pd.DataFrame(index=opmModeSet,columns=cycles)
    opmFraction=pd.DataFrame(index=opmModeSet,columns=cycles)
        
    for c in cycles:
        fzpc=fzppercycle[c]
        vehc=vehpercycle[c]
        if len(vehc)==0:
            print('no vehicle records in '+str(inter)+str(direc)+' cycle'+str(c))
           
        # calculate links inputs
        v_avg=np.mean(fzpc['speed'])
        q=len(vehc)
        # built links inputs
        ls_linkID.append(c+1)
        ls_linkVolume.append(q)
        ls_linkAvgSpeed.append(v_avg)
        
        # calculate link source type inputs and op mode
        types=fzpc['vehtype'].unique()
        typecounts=fzpc['vehtype'].value_counts()
        typefraction=pd.Series(index=types)
        for t in types:
            # link source type inputs
            typefraction[t]=typecounts[t]/sum(typecounts)
            lst_linkID.append(c+1)
            lst_sourceTypeID.append(typemap[t])
            lst_fraction.append(typefraction[t])
            
            # op mode inputs
            
            # get parameters for calculating VSP
            sourceID=typemap[t]
            params=ts.getVSPparams(sourceID,vspParams)
            
            fzpct=fzpc.loc[fzpc['vehtype']==t] # filter by type
            
            ### assign vsp and speed bins
            
            # separate by brake/go
            fzp_br=fzpct.loc[fzpct['accel']<=brake_threshold]
            # TODO: add 3-second brake rate condition
            opmSeconds[c][0]=len(fzp_br.index)
            fzp_go=fzpct.loc[fzpct['accel']>brake_threshold]
    
            # count idle mode vehicle seconds
            fzp_idle=fzp_go.loc[(fzp_go['speed']>=-m2k) & (fzp_go['speed']<=m2k)]
            opmSeconds[c][1]=len(fzp_idle.index)
    
            ## accelerating/cruising modes/coasting modes
    
            fzp_spd_vsp={}
            for i in range(len(spdbins)-1):
                # separate by speed bins
                fzp_spd=fzp_go.loc[(fzp_go['speed']>=spdbins[i]) & (fzp_go['speed']<spdbins[i+1])]
                
                # get vsp for each second
                fzp_spd['vsp']=ts.getVSP(fzp_spd['speed']*kph2mps,fzp_spd['accel'],grade,params,sourceID)
                
                # separate by vsp bins
                fzp_spd_vsp[i]={}
                for j in range(1,len(vspbins[i])):
                    fzp_spd_vsp[i][j]=fzp_spd.loc[(fzp_spd['vsp']>=vspbins[i][j-1]) & (fzp_spd['vsp']<vspbins[i][j])]
                    
                # first and last vsp bin is open
                fzp_spd_vsp[i][0]=fzp_spd.loc[fzp_spd['vsp']<vspbins[i][0]]
                fzp_spd_vsp[i][len(vspbins[i])]=fzp_spd.loc[fzp_spd['vsp']>=vspbins[i][-1]]
            
            # last speed bin is open
            
            fzp_spd=fzp_go.loc[fzp_go['speed']>=spdbins[-1]]
            # get vsp for each second
            fzp_spd['vsp']=ts.getVSP(fzp_spd['speed']*kph2mps,fzp_spd['accel'],grade,params,sourceID)
            
            # separate by vsp bins
            lastspd=len(spdbins)-1
            fzp_spd_vsp[lastspd]={}
            for j in range(1,len(vspbins[lastspd])):
                fzp_spd_vsp[lastspd][j]=fzp_spd.loc[(fzp_spd['vsp']>=vspbins[lastspd][j-1]) & (fzp_spd['vsp']<vspbins[lastspd][j])]
                
            # first and last vsp bin is open
            fzp_spd_vsp[lastspd][0]=fzp_spd.loc[fzp_spd['vsp']<vspbins[lastspd][0]]       
            fzp_spd_vsp[lastspd][len(vspbins[lastspd])]=fzp_spd.loc[fzp_spd['vsp']>=vspbins[lastspd][-1]]
            
            # fill in opm seconds with each op mode bin
            i=0
            mi=2
            while i<=lastspd:
                lastvsp=len(vspbins[i])
                j=0
                while j<=lastvsp:
                    m=opmModeSet[mi]
                    if not fzp_spd_vsp[i][j].empty:
                        opmSeconds[c][m]=len(fzp_spd_vsp[i][j].index)
                    else:
                        opmSeconds[c][m]=0
                    j=j+1
                    mi=mi+1
                i=i+1
            
            # get op mode fractions using op mode seconds
            for m in opmModeSet:
                if sum(opmSeconds[c])==0:
                    print(str(inter)+str(d)+' cycle'+str(c)+' type'+str(t)+'doesn\'t have any operating mode seconds.')
                opmFraction[c][m]=opmSeconds[c][m]/sum(opmSeconds[c])

            
            # populate op mode counts
            for p in opmModeSetProcesses:
                opm_type.extend([typemap[t]]*modes)
                opm_process.extend([p]*modes)
                opm_link.extend([c+1]*modes)
                opm_mode.extend(opmModeSet)
                opm_fraction.extend(list(opmFraction[c]))
            for p in opmAllRunningProcesses:
                opm_type.append(typemap[t])
                opm_process.append(p)
                opm_link.append(c+1)
                opm_mode.append(opmAllRunning)
                opm_fraction.append(1)
    
    ### populate MOVES input dataframes and export
            
    # links
    ls['linkID']=ls_linkID
    lsSize=len(ls)
    ls['linkLength']=[lengthtotal]*lsSize
    ls['linkVolume']=ls_linkVolume
    ls['linkAvgSpeed']=ls_linkAvgSpeed
    ls['countyID']=[county]*lsSize
    ls['zoneID']=[zone]*lsSize
    ls['roadTypeID']=[roadType]*lsSize
    
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

    return ls, lst, opm

#%%import vissim data and get operating modes

fzps=ts.listFiles('fzp')
lsas=ts.listFiles('lsa')

ls={}
lst={}
opm={}

# for each simulation
for s in range(len(fzps)):
    
    # import vissim data
    fzp=ts.import_fzp(fzps[s])
    lsa=ts.import_lsa(lsas[s])

    for p in periods:
        
        # start templates for links, link source types, and operating mode distributions
        ls[p]=ts.makeLinksTemplate()
        lst[p]=ts.makeLinkSourceTypeTemplate()
        opm[p]=ts.makeOpModeTemplate()
        
        lstest=ts.makeLinksTemplate()

        # filter by period
        fzp_p=ts.filterbyperiod(fzp,p)
        lsa_p=ts.filterbyperiod(lsa,p)    
            
        inter_index=1
        for i in inters: # for each intersection
            direc_index=1
            for d in direcs: # for each direction (NB, SB)
    
                # Dunbar has no northbound and Hwy401 has no southbound available (end intersections)
                if not ((i==inters[0] and d==direcs[0]) or (i==inters[-1] and d==direcs[-1])):
                    
                    ts0=timeit.default_timer()

                    # get operating mode distributions
                    ls0,lst0,opm0=getOpMode(i,d,mov,fzp_p,lsa_p)
                    
                    ts3=timeit.default_timer()
                    
                    # change link number to reflect simulation, intersection, and direction
                    change=s*100+inter_index*10000+direc_index*1000
                    ls0['linkID']+=change
                    lst0['linkID']+=change
                    opm0['linkID']+=change
                    
                    # append to existing set of operating modes
                    ls[p]=ls[p].append(ls0)
                    lst[p]=lst[p].append(lst0)
                    opm[p]=opm[p].append(opm0)
                    
                    tsn=timeit.default_timer()
                    print('Generated MOVES inputs for '+p+', '+i+d+', simulation'+str(s)+' in %.fs' %(tsn-ts0))
                                
                direc_index+=1
            inter_index+=1

# export to csv for all periods
for p in periods:
    
    # save path
    savepath='MOVESinputs/'+str(mov)+'/'+str(p)
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        
    # export to csv
    ls[p].to_excel(savepath+'/links.xlsx',index=False)
    lst[p].to_excel(savepath+'/linkSourceTypes.xlsx',index=False)
    opm[p].to_excel(savepath+'/OpModeDistribution.xlsx',index=False)

#%% timer
runtime=timeit.default_timer()-t_start
print('Total runtime: %.2f seconds' %runtime)




