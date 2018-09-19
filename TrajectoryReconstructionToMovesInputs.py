# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 12:00:35 2018

@author: Anjie
"""

import pandas as pd
import numpy as np
import json, os
import matplotlib.pyplot as plt
from time import time

import toolSpace as ts 

start_t=time() 

# manual inputs
scenario='advDetector25m'
advpos=2 # position number of advance detector [2:25m, 3:50m, 4:75m, 5:100m]
majorOnly=True
realPeriods=True

#%%

# file directories
outpath='MOVESinputs/'+scenario+'/'
figoutpath='TrajectoryReconstructionPlots/'+scenario
if not os.path.exists(figoutpath):
    os.makedirs(figoutpath)

# intersections, directions
inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
direcs=['NB','SB']

# list files for simulation data
lsas=ts.listFiles('lsa')
mers=ts.listFiles('mer')
rsrs=ts.listFiles('rsr')
fzps=ts.listFiles('fzp')

# runSpecs
runSpecs=json.load(open('runSpecs.json'))
periods=json.load(open('periods-variations.json'))
if realPeriods:
    periods=json.load(open('periods-real.json'))

# modelSpecs
modelSpecs=json.load(open('modelSpecs.json'))
if majorOnly:
    modelSpecs=json.load(open('modelSpecs-MajorOnly.json'))


#%% parameters for trajectory reconstruction
cruisepercentile=80
dt=1 # time step [s]
adj=.1 # parameter adjustment step for curve fitting
stoplag=1.5 # time between vehicle detection at stopbar and vehicle coming to a complete stop
startlag=2 # time between vehicle starting to accelerate and vehicle no longer being detected at stopped position
stoptime=2 # minimum time of detection for a vehicle that stopped at the detector position
pstopcruise=.5 # if a vehicle fully stops, its average speed is at most this fraction of the cruise speed

t_tol=1 # acceptable difference between detected time and fitted curve time at the advance position
dec_max=-10 # maximum deceleration [m/s^2]
max_iter=10 # maximum iterations for curve fitting
max_accel_dist=100 # maximum distance for building acceleration curve [m]
max_decel_dist=300 # maximum distance for building deceleration curve [m]
# original deceleration curve parameters from literature
k3o=0.005 
k4o=0.154
k5o=0.493
# acceleration curve parameters from literature
beta0o=3.369
beta1o=0.072
# car-following parameters: time interval and effective vehicle length
leffs=pd.Series([5.5,13,11.5,13.5],index=[100,150,200,210]) # effective length [m]
hDis=2 # maximum headway for queue discharge

#%% for generating MOVES inputs

# load VSP parameters. Source: Population and Activity of On-Road Vehicles in MOVES2014, page 103
vspParams=pd.read_excel('PhysicsTable.xlsx',sheet_name='Revised',dtype={'sourceTypeID':int,'cA':float,'cB':float,'cC':float,'mass':float})
# vehicle type ID mapping from Vissim to MOVES
typemap=ts.getVehTypeMap()
# parameters for calculating VSP by vehicle type
paramdict={}
for source in typemap:
    paramdict[source]=ts.getVSPparams(source,vspParams)

# conversion variables
m2k=1.60934 # km per mile
kph2mps=1/3.6 # m/s per km/h
mph2mps=m2k*kph2mps
grade=0

# specs for MOVES
county=runSpecs['county']
zone=runSpecs['zone']
roadType=runSpecs['roadType']
hourDayID=runSpecs['hourDayID']

opmFormat=json.load(open('opModeFormat.json'))
opmModeSetProcesses=opmFormat['processes']['modeSet']
opmAllRunningProcesses=opmFormat['processes']['allRunning']
opmModeSet=opmFormat['opModeSets']['modeSet']
modes=len(opmModeSet)
opmAllRunning=opmFormat['opModeSets']['allRunning'][0]
spdbins=opmFormat['speedbins']
for n in spdbins: # convert speeds from miles/hr to km/hr
    n=n*m2k
    
vspbins={}
for n in opmFormat['vspbins']:  
    vspbins[int(n)]=opmFormat['vspbins'][n]
    
brake_threshold=-2*mph2mps

#%% filter data, reconstruct trajectories, and append to trajectory records
for s in range(len(lsas)): # for each simulation
    
    # import vissim data
    lsa=ts.import_lsa(lsas[s])
    mer,merDetail=ts.import_mer(mers[s])
    rsr=ts.import_rsr(rsrs[s])
    fzp=ts.import_fzp(fzps[s])
    merEn=mer.loc[mer['t(Entry)']>0] # divide detector data into dectector on and detector off data
    merEx=mer.loc[mer['t(Exit)']>0]
    tg=lsa.loc[lsa['newState'].str.contains('green')] # filter signal timing to green times
    
    detn_all=merDetail[3]
    merDetail['intersection']=detn_all/1000
    merDetail['direction']=detn_all%1000/100
    merDetail['lane']=detn_all%100/10
    merDetail['position']=detn_all%10
    
    merDetail=merDetail.astype({'intersection':int,'direction':int,'lane':int,'position':int})

    for p in periods: # for each period (hour)
#        print('Starting '+p)
        rsrp=ts.filterbyperiod(rsr,p)
        fzpp=ts.filterbyperiod(fzp,p)
        
        # templates for outputs - MOVES inputs
        ls_linkID=[]
        ls_linkVolume=[]
        ls_linkAvgSpeed=[]
        ls_linkLength=[]
        lst_linkID=[]
        lst_sourceTypeID=[]
        lst_fraction=[]
        opm_type=[]
        opm_link=[]
        opm_process=[]
        opm_mode=[]
        opm_fraction=[]
        
        inter_index=1 # intersection index (for link ID, detector identification)
        for i in inters: # for each intersection
            tgi=tg.loc[tg['SC']==modelSpecs[i]['SignalController']] # filter to signal controller
            merDetail_i=merDetail.loc[merDetail['intersection']==inter_index] # filter to intersection
            
            direc_index=1
            for d in direcs: # for each direction
                                
                # Dunbar northbound and Hwy401 southbound are out of scope (end intersections)
                if not ((i=='Dunbar' and d=='NB') or (i=='Hwy401' and d=='SB')):
                    time_id0=time()
#                    print('Starting next segment')
                    
                    merDetail_id=merDetail_i.loc[merDetail_i['direction']==direc_index] # filter to direction
                    
                    f=ts.startFig('Lead Trajectories of '+str(p)+str(i)+str(d),'Time','Position',4)
                    
                    # segment model specs
                    segSpecs=modelSpecs[i][d]
                    
                    # links in the segment
                    intlink=segSpecs['links']['intersection']['thru']
                    approachlinks=segSpecs['links']['approach']['thru']
                    links=approachlinks+intlink
                                        
                    # link lengths
                    lengthMap=ts.getLengthMap(links)
                    segLength=sum(lengthMap)
                    
                    # simulated trajectories from vehicle records
                    fzppid=fzpp.loc[fzpp['link'].isin(links)] # vehicle records for the intersection and approach
                    ## find absolute position of vehicle with end of last major intersection as reference
                    fzpx=0 #cumulative distance from reference
                    # main segment positions
                    for l in approachlinks:
                        fzppid.loc[fzppid['link']==l,'abspos']=fzppid['pos']+fzpx
                        fzpx +=lengthMap[l]
                    # intersection positions
                    fzppid.loc[fzppid['link']==intlink,'abspos']=fzppid['pos']+fzpx
                    # plot trajectories vehicle by vehicle
                    traj_vehicles=fzppid['no'].unique() # all vehicles numbers
                    for n in traj_vehicles:
                        veh=fzppid.loc[fzp['no']==n]
                        plt.plot(veh['sec'],veh['abspos'],'grey')

                    # filter signal timing data to segment
                    tgid=tgi.loc[lsa['SG']==segSpecs['StartingSignalGroup']] # filter signal timing to signal group
                    tgid=tgid.reset_index(drop=True) # reset index for filtering
                    tgidp=ts.filterbyperiod(tgid,p) # filter signal timing to greens starting this period
                    if p!='evening2':
                        tgidp=tgidp.append(tgid.loc[tgidp.index[-1]+1]) # add on the next green start after the last green start in this period to close the last bin
                    tgidp=tgidp['tsim'] # only the timing data is needed
                    if p=='evening2': # for last period, there is no next green after the period (simulation stopped)
                        tgidp=tgidp.append(pd.Series(merEx['t(Exit)'].iloc[-1]),ignore_index=True)
                    tgidp=tgidp.reset_index(drop=True) # reset index for iterating through cycles
                    
                    rsrpid=rsrp.loc[rsrp['No.']==segSpecs['TravelTimeMeas']] # filter travel time measurements 
                    
                    # get detector numbers
                    detn={}                
                    lanes=max(merDetail_id['lane'])
                    detn_start=inter_index*1000+direc_index*100
                    detn['adv']=list(detn_start+advpos+pd.Series(range(1,lanes+1))*10)
                    detn['stop']=list(detn_start+1+pd.Series(range(1,lanes+1))*10)
                    detn['end']=list(detn_start+pd.Series(range(1,lanes+1))*10)
                    
                    xdet={} # get positions of detectors
                    for det in detn:
                        xdet[det]=int(merDetail.loc[merDetail[3]==detn[det][0],9].iloc[0])
                    # for stopbar detector
                    posOffset=sum(lengthMap.loc[pd.Series(approachlinks).iloc[0:-1]]) # offset distance is length of upstream links
                    xdet['stop']=xdet['stop']+posOffset
                    # for advance detector
                    detlink=merDetail.loc[merDetail[3]==detn['adv'][1],5].iloc[0] # link detector is on
                    detlinkindex=links.index(detlink)
                    posOffset=sum(lengthMap.loc[pd.Series(approachlinks).iloc[0:detlinkindex]]) # offset distance is length of upstream links up to link of detector
                    xdet['adv']=xdet['adv']+posOffset
                    # end of intersection
                    xdet['end']=xdet['stop']+float(lengthMap[intlink]+2) #add on intersection span

                    # get cruising speed
                    merExid_stop=merEx.loc[merEx['Measurem.'].isin(detn['stop'])] # stopbar exit detections
                    vehicles=merExid_stop['VehNo'].unique()
                    merAdv=merEx.loc[merEx['Measurem.'].isin(detn['adv'])]
                    merStop=merEx.loc[merEx['Measurem.'].isin(detn['stop'])]
                    timeStop=[]
                    timeAdv=[]
                    for n in vehicles:
                        timeStop.append(merStop.loc[merStop['VehNo']==n,'t(Exit)'].iloc[0])
                        timeAdv.append(merAdv.loc[merAdv['VehNo']==n,'t(Exit)'].iloc[0])
                    timeStop=pd.Series(timeStop)
                    timeAdv=pd.Series(timeAdv)
    
                    time_adv2stop=timeStop-timeAdv # times to travel from advance to stopbar
                    distance_adv2stop=pd.Series([xdet['stop']-xdet['adv']]*len(time_adv2stop)) # distance from advance to stop as an array
                    speeds_adv2stop=(distance_adv2stop)/(time_adv2stop) # speeds based on distance and times
                    v0=np.percentile(speeds_adv2stop,cruisepercentile) # cruise speed as percentile of speeds
                    print('Cruise speed: %.1f' %v0)
                    
                    # iterate through cycles in this period
                    cycles=range(len(tgidp)-1) # the added on green start is not counted - just for closing the bin               
                    
                    for c in cycles:
#                        print('Processing '+str(p)+' '+i+' '+d+' cycle '+str(c))
                        cyclevehn=[]
                        link=[]
                        speed=[]
                        accel=[]
                        vehtype=[]
                        
#                        lanes=len(detn['stop'])
                        for l in range(lanes):
                            
                            detectorStop=detn['stop'][l]

                            # find vehicles detected in this lane during this cycle at the stopbar
                            merExid_stopbar=merEx.loc[merEx['Measurem.']==detectorStop]
                            merExid_stopbar_cycle=merExid_stopbar.loc[(merExid_stopbar['t(Exit)']>=tgidp[c])&(merExid_stopbar['t(Exit)']<=tgidp[c+1])]
                            vehicles_cycle=merExid_stopbar_cycle['VehNo'].unique()
                            # find vehicles that went thru end of intersection
                            merExid_end=merEx.loc[merEx['Measurem.'].isin(detn['end'])]
                            vehicles_thru=merExid_end['VehNo'].unique()
                            
                            # vehicles in this lane and going thru
                            vehicles=np.intersect1d(vehicles_cycle,vehicles_thru)
                            
                            if len(vehicles)>0:
                            
                                # filter detections of thru vehicles of this cycle in this lane
                                # (so that adv detections before the cycle starts are included as long as the vehicle passes the stopbar during the cycle.)
                                merEnc=merEn.loc[merEn['VehNo'].isin(vehicles)]
                                merExc=merEx.loc[merEx['VehNo'].isin(vehicles)]
                                
                                # data for stopbar detection of this lanes
                                detectionStopEn=merEnc.loc[merEnc['Measurem.']==detectorStop]
                                detectionStopEx=merExc.loc[merExc['Measurem.']==detectorStop]
                                # data for advance detection of all lanes
                                detectionAdvEn=merEnc.loc[merEnc['Measurem.'].isin(detn['adv'])]
                                detectionAdvEx=merExc.loc[merExc['Measurem.'].isin(detn['adv'])]
                                # data for end of intersection detection of all lanes
                                detectionEndEn=merEnc.loc[merEnc['Measurem.'].isin(detn['end'])]
                                detectionEndEx=merExc.loc[merExc['Measurem.'].isin(detn['end'])]
                                
                                # compose detector data to combine entry and exit times into one data entry for each vehicle
                                stopt1=[]
                                stopt2=[]
                                advt1=[]
                                advt2=[]
                                endt1=[]
                                endt2=[]
                                vehn=[]
                                vtype=[]
                                for n in vehicles:
                                    # compile detector data
                                    vehn.append(n)
                                    vtype.append(detectionStopEn.loc[detectionStopEn['VehNo']==n,'Vehicle type'].iloc[0])
                                    stopt1n=detectionStopEn.loc[detectionStopEn['VehNo']==n,'t(Entry)'].iloc[0]
                                    stopt2n=detectionStopEx.loc[detectionStopEx['VehNo']==n,'t(Exit)'].iloc[0]
                                    advt1n=detectionAdvEn.loc[detectionAdvEn['VehNo']==n,'t(Entry)'].iloc[0]
                                    advt2n=detectionAdvEx.loc[detectionAdvEx['VehNo']==n,'t(Exit)'].iloc[0]
                                    endt1n=detectionEndEn.loc[detectionEndEn['VehNo']==n,'t(Entry)'].iloc[0]
                                    endt2n=detectionEndEx.loc[detectionEndEx['VehNo']==n,'t(Exit)'].iloc[0]
                                    stopt1.append(stopt1n)
                                    stopt2.append(stopt2n)
                                    advt1.append(advt1n)
                                    advt2.append(advt2n)
                                    endt1.append(endt1n)
                                    endt2.append(endt2n)

                                    plt.plot([stopt1n,stopt2n],[xdet['stop']]*2,'orange')
                                    plt.plot([advt1n,advt2n],[xdet['adv']]*2,'orange')
                                    plt.plot([endt1n,endt2n],[xdet['end']]*2,'orange')

                                # dataframe of detector data
                                detection=pd.DataFrame({'vehn':vehn,'advt1':advt1,'advt2':advt2,'stopt1':stopt1,'stopt2':stopt2,'endt1':endt1,'endt2':endt2,'vtype':vtype})

                                detection=detection.sort_values(by='stopt1') # sort by time of reaching stopbar
                                detection=detection.reset_index(drop=True)
                                
                                leadTraj={} # trajectory of lead vehicle
                                
                                leaddetection=detection.iloc[0]
                                stopbarWait=leaddetection['stopt2']-leaddetection['stopt1']
                                
                                # check if there was a queue - how long the first vehicle stopped at the stopbar
                                if stopbarWait>stoptime: # if the vehicle had stopped, construct deceleration and acceleration curves for the lead vehicle
#                                    print('lead vehicle stopped')
                                    followcheck=True
                                    
                                    tdelta=9999 # time difference for fitting curve to detector data (to be minimized)
                                    
                                    # deceleration curve
                                    
                                    # start with original values
                                    k3=k3o
                                    k4=k4o
                                    k5=k5o
                                    itercount=0
                                    while abs(tdelta)>t_tol and itercount<max_iter: # time difference between advance detection and extrapolated cruise trajectory at advance
                                        vdec=[v0] # speed profile
                                        adec=[0] # deceleration profile
                                        xdec=[0] # positions 
                                        tdec=[0] # time
                                        while vdec[-1]>0 and xdec[-1]<max_decel_dist: # until speed has reached zero 
                                            dec_rate=-k3*vdec[-1]**2+k4*vdec[-1]+k5
                                            if dec_rate<dec_max:
                                                dec_rate=dec_max
                                            adec.append(dec_rate) # using deceleration curve from literature
                                            vdec.append(vdec[-1]-(adec[-1]+adec[-2])/2*dt) # get speed for each time step
                                            xdec.append(xdec[-1]+(vdec[-1]+vdec[-2])/2*dt+(adec[-1]+adec[-2])/4*dt**2) # append to positions
                                            tdec.append(tdec[-1]+dt) # append to times
                                
                                        # shift deceleration trajectory to match end of deceleration with beginning of stopped trajectory
                                        tdeci=pd.Series(tdec)+leaddetection['stopt1']-tdec[-1]+stoplag
                                        xdeci=pd.Series(xdec)+xdet['stop']-xdec[-1]
                                        tdelta=leaddetection['advt1']-tdeci[abs(xdeci-xdet['adv']).idxmin()]
                                        if xdeci[0]<xdet['adv']: # if deceleration started after advance detection
                                            tcruise=(xdeci[0]-xdet['adv'])/vdec[0] # time spent cruising after advance detection before deceleration
                                            textra=tdeci[0]-tcruise # extrapolated time of reaching advance detector from beginning of deceleration
                                            tdelta=leaddetection['advt1']-textra # difference between advance detection time and extrapolated time
                                        if tdelta>0: # if curve is earlier than detection, decelerate faster
                                            dec_change='faster'
                                            k3=k3-k3o*adj 
                                            k4=k4+k4o*adj
                                            k5=k5+k5o*adj
                                        else: # if curve is later than detection, decelerate slower
                                            dec_change='slower'
                                            k3=k3+k3o*adj 
                                            k4=k4-k4o*adj
                                            k5=k5-k5o*adj
                                        
                                        itercount+=1
#                                        print('itercount='+str(itercount)+', decelerated '+dec_change)
                                    
                                    # stopped curve
                                    tStopEnd=leaddetection['stopt2']
                                    
                                    # acceleration curve
                                    tdelta=9999 # to be minimized
                                    
                                    # start with original values
                                    beta0=beta0o
                                    beta1=beta1o
                                    
                                    itercount=0
                                    while abs(tdelta)>t_tol and itercount<max_iter:
                                        #initial values at the beginning of the curve
                                        vacc=[0]
                                        aacc=[0]
                                        xacc=[0]
                                        tacc=[0]
                                        
                                        while vacc[-1]<v0 and xacc[-1]<max_accel_dist:
                                            aacc.append(beta0-beta1*vacc[-1]) # acceleration function from literature
                                            vacc.append(vacc[-1]+(aacc[-1]+aacc[-2])/2*dt)
                                            xacc.append(xacc[-1]+(vacc[-1]+vacc[-2])/2*dt+(aacc[-1]+aacc[-2])/4*dt**2)
                                            tacc.append(tacc[-1]+dt)
                                                                                      
                                        # shift acceleration curve to start at the end of the stopped curve
                                        tacci=pd.Series(tacc)+tStopEnd-startlag
                                        xacci=pd.Series(xacc)+xdet['stop']

                                        tdelta=leaddetection['endt2']-tacci[abs(pd.Series(xacci)-xdet['end']).idxmin()] # shortest distance between acceleration curve and end detection
                                        if xacci.iloc[-1]>xdet['end']: # if acceleration ended after end detection
                                            tcruise=(xdet['end']-xacci.iloc[-1])/vacc[-1] # time spent cruising after acceleration before reaching end detector
                                            textra=tacci.iloc[-1]+tcruise # extrapolated time of reaching end detector from end of acceleration
                                            tdelta=leaddetection['endt2']-textra # difference between end detection time and extrapolated time
                                        if tdelta>0: # if curve is earlier than detection, accelerate slower
                                            beta0-=beta0o*adj
                                            beta1+=beta1o*adj

                                        else: # if curve is later than detection, accelerate faster
                                            beta0+=beta0o*adj
                                            beta1-=beta1o*adj
                                        
                                        itercount+=1
#                                        print('itercount='+str(itercount))
                                        
                                        
                                    # points for idle trajectory
                                    v_idle=[0]*int(stopbarWait)
                                    a_idle=[0]*int(stopbarWait)
                                    
                                    # compile curves into lead trajectory
                                    leadTraj['t']=list(tdeci)+list(tacci)
                                    leadTraj['x']=list(xdeci)+list(xacci)
                                    leadTraj['a']=list(adec)+v_idle+list(aacc)
                                    leadTraj['v']=list(vdec)+a_idle+list(vacc)
                                    leadTraj['vehn']=detection.loc[0,'vehn']
                                    leadTraj['vtype']=detection.loc[0,'vtype']
                                    
                                    ### following trajectories
                                    tstop=detection.loc[0,'stopt1'] # time preceding vehicle stops at stopbar
                                    tgo=detection.loc[0,'stopt2'] # time preceding vehicle leaves from stopbar
                                    xstop=xdet['stop']
                                    
                                    traj={0:{'x':leadTraj['x'],'t':leadTraj['t'],'v':leadTraj['v'],'a':leadTraj['a'],'vehn':leadTraj['vehn'],'vtype':leadTraj['vtype']}}
                                    for j in detection[1:].index:
                                        traj[j]={'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']} # fill in vehicle number and type first
                                        vapproach=(xdet['stop']-xdet['adv'])/(detection.loc[j,'stopt1']-detection.loc[j,'advt1']) # average speed approaching the intersection
                                        h=detection.loc[j,'stopt2']-detection.loc[j-1,'stopt2'] # headway wrt preceding vehicle at stopbar 
                                        leff=leffs[detection.loc[j-1,'vtype']] # get effective length of preceding vehicle based on vehicle type
                                        
                                        # find T (time interval) for Newell's model - time lag of trajectories
                                        t2=detection.loc[j,'stopt2']
                                        xlead=xdet['stop']+leff # effective position of preceding vehicle corresponding the following vehicle at the stopbar
                                        index=abs(np.array(traj[j-1]['x'])-xlead).argmin() # find point in lead closest to stopbar
                                        t1=traj[j-1]['t'][index] # time when preceding vehicle is at this position
                                        if xlead>traj[j-1]['x'][-1]: # if lead acceleration curve ends before effective lead position,
                                            t1=traj[j-1]['t'][-1]+(xdet['stop']+leff-traj[j-1]['x'][-1])/v0 # then extrapolate using cruise speed
                                        tlag=t2-t1 # time interval for Newell's model
                                        
                                        if followcheck:
                                        
                                            if vapproach<pstopcruise*v0: #check if vehicle stopped in queue
#                                                print('vehicle stopped in queue')
                                                # offset end of stopped curve of preceding vehicle to get end of stopped curve for current following vehicle
                                                tgo+=tlag
                                                xstop-=leff
                                                # shift deceleration curve spatially so that it ends where the stopped curve starts
                                                xdeci=pd.Series(xdec).add(xstop-xdec[-1])
                                                # shift deceleration curve temporally
                                                if xdeci[0]>xdet['adv']: # if deceleration curve starts after the advance location
                                                    tdeci=pd.Series(tdec).add(detection.loc[j,'advt1']+(xdeci[0]-xdet['adv'])/v0) # match extrapolated upstream position to advance detector
                                                else: # if deceleration curve starts before the advance location
                                                    tdelta=tdec[abs(np.array(xdeci-xdet['adv'])).argmin()]-tdec[0] # find the point closest to the advance detector
                                                    tdeci=pd.Series(tdec).add(detection.loc[j,'advt1']-tdelta) # shift deceleration curve to match this point to advance detector
                                                
                                                #acceleration curve

                                                # shift acceleration curve to start at the end of the stopped curve
                                                tacci=[tt+tgo-startlag for tt in tacc]
                                                xacci=[xx+xstop for xx in xacc]
                                                
                                                # idling 
                                                t_idle=tgo-tdeci.iloc[-1]
                                                v_idle=[0]*int(t_idle)
                                                a_idle=[0]*int(t_idle)
                                                
                                                # put together deceleration and acceleration curves
                                                traj[j]['t']=list(tdeci)+list(tacci)
                                                traj[j]['x']=list(xdeci)+list(xacci)
                                                traj[j]['v']=list(vdec)+v_idle+list(vacc)
                                                traj[j]['a']=list(adec)+a_idle+list(aacc)
                                                
                                                                                        
                                            elif h<hDis: # didn't stop but followed closely
#                                                print('Vehicle followed closely')
                                                # values at start of trajectory
                                                t0=detection.loc[j,'advt1']
                                                traj[j]={'t':[t0],'x':[xdet['adv']],'a':[0],'v':[v0],'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']}
                                                for k in range(1,len(traj[j-1]['t'])): # offset point by point
                                                    # build following trajectory as time-space offset of preceding trajectory (Newell's model)
                                                    tfol=traj[j-1]['t'][k]+tlag # time offset
                                                    xfol=traj[j-1]['x'][k]-leff # position offset
                                                    traj[j]['x'].append(xfol)
                                                    tfolcruise=(traj[j]['x'][k]-traj[j]['x'][k-1])/v0+traj[j]['t'][k-1]
                                                    if tfolcruise>tfol: # if cruising gets the vehicle to the follow position later than the time offset, then assume cruising instead
                                                        tfol=tfolcruise
                                                    traj[j]['t'].append(tfol)
                                                    # calculate speed and acceleration per point increment
                                                    traj[j]['v'].append(traj[j]['x'][-1]-traj[j]['x'][-2]) # speed as change in position
                                                    traj[j]['a'].append(traj[j]['v'][-1]-traj[j]['v'][-2]) # acceleration as change in speed
                                            
                                            else:
#                                                print('vehicle in free flow')
                                                followcheck=False # vehicle are no longer following - go to free flow reconstruction
                                                # straight trajectories from advance detector to stop detector
                                                traj[j]={'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']}
                                                ti=detection.loc[j,'advt2']
                                                tf=detection.loc[j,'stopt2']
                                                tflow=tf-ti
                                                xi=xdet['adv']
                                                xf=xdet['stop']
                                                xflow=xf-xi
                                                vflow=xflow/tflow
                                                traj[j]['t']=np.arange(ti,ti+int(tflow)+dt,dt)
                                                points=len(traj[j]['t'])
                                                traj[j]['x']=np.linspace(xi,xf,points)
                                                traj[j]['v']=[v0]*points
                                                traj[j]['a']=[0]*points
                                                                                                   
                                        else: #free flow vehicles
                                            # straight trajectories from advance detector to stop detector
#                                            print('vehicle in free-flow')
                                            traj[j]={'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']}
                                            ti=detection.loc[j,'advt2']
                                            tf=detection.loc[j,'stopt2']
                                            tflow=tf-ti
                                            xi=xdet['adv']
                                            xf=xdet['stop']
                                            xflow=xf-xi
                                            vflow=xflow/tflow
                                            traj[j]['t']=np.arange(ti,ti+int(tflow)+dt,dt)
                                            points=len(traj[j]['t'])
                                            traj[j]['x']=np.linspace(xi,xf,points)
                                            traj[j]['v']=[vflow]*points
                                            traj[j]['a']=[0]*points
                                            
                                        plt.plot(traj[j]['t'],traj[j]['x'],'g--')
                               
                                else: # if first vehicle didn't stop, then trajectory is assumed to be straight from advance to stopbar
                                    # lead trajectory
#                                    print('lead vehicle didn\'t stop')
                                    ti=detection.loc[0,'advt2']
                                    tf=detection.loc[0,'stopt2']
                                    tflow=tf-ti
                                    xi=xdet['adv']
                                    xf=xdet['stop']
                                    xflow=xf-xi
                                    vflow=xflow/tflow
                                    leadTraj['t']=np.arange(ti,ti+int(tflow)+dt,dt)
                                    points=len(leadTraj['t'])
                                    leadTraj['x']=np.linspace(xi,xf,points)
                                    leadTraj['v']=[vflow]*points
                                    leadTraj['a']=[0]*points
                                    leadTraj['vehn']=detection.loc[0,'vehn']
                                    leadTraj['vtype']=detection.loc[0,'vtype']
                                    
                                    traj={0:{'t':leadTraj['t'],'x':leadTraj['x'],'a':leadTraj['a'],'v':leadTraj['v'],'vehn':leadTraj['vehn'],'vtype':leadTraj['vtype']}}
                                    for j in detection.index[1:]:
                                        traj[j]={'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']} # fill in vehicle number and type first
                                        ti=detection.loc[j,'advt2']
                                        tf=detection.loc[j,'stopt2']
                                        tflow=tf-ti
                                        xi=xdet['adv']
                                        xf=xdet['stop']
                                        xflow=xf-xi
                                        vflow=xflow/tflow
                                        traj[j]['t']=np.arange(ti,ti+int(tflow)+dt,dt)
                                        points=len(traj[j]['t'])
                                        traj[j]['x']=np.linspace(xi,xf,points)
                                        traj[j]['v']=[vflow]*points
                                        traj[j]['a']=[0]*points
                                        plt.plot(traj[j]['t'],traj[j]['x'],'g--')
                                
                                
                                
                                # plot lead trajectories
                                
                                plt.plot(leadTraj['t'],leadTraj['x'],'r--')
                                
                                # complete reconstructed trajectories for lane and add to cycle collection
                                for j in traj:
#                                    plt.plot(traj[j]['t'],traj[j]['x'],'bo')

                                    # filling for before the reconstructed trajectory (constant speed)
                                    xfill=traj[j]['x'][0]
                                    tfill=xfill/v0
                                    fillpoints=int(tfill) # number points, i.e. seconds, to fill in before the trajectory
                                    vfillList1=[v0]*fillpoints # array of speeds to fill in
                                    afillList1=[0]*fillpoints # array of accelerations to fill in 
                                    
                                    # after the acceleration, check if there's a gap or overshoot
                                    gap=xdet['end']-traj[j]['x'][-1]                          
                                    if gap>0:
                                        # filling for after the reconstructed trajectory (constant speed)
                                        xfill=gap
                                        tfill=xfill/v0
                                        fillpoints=int(tfill) # number points, i.e. seconds, to fill in before the trajectory
                                        vfillList2=[v0]*fillpoints # array of speeds to fill in
                                        afillList2=[0]*fillpoints # array of accelerations to fill in
                                    else:
                                        
                                        # remove acceleration curve after end of intersection
                                        acc_remove=len(np.where(np.array(traj[j]['x'])>xdet['end']))
                                        traj[j]['v']=traj[j]['v'][:len(traj[j]['v'])-acc_remove]
                                        traj[j]['a']=traj[j]['a'][:len(traj[j]['a'])-acc_remove]
                                        vfillList2=[]
                                        afillList2=[]
                                    
                                    totalpoints=len(vfillList1)+len(traj[j]['v'])+len(vfillList2) # number of points in whole trajectory
                                    
                                    # append to table of trajectory records for the cycle
                                    speed.extend(vfillList1+traj[j]['v']+vfillList2) # join trajectories before, during, and after reconstructed portion
                                    accel.extend(afillList1+traj[j]['v']+afillList2) # join trajectories before, during, and after reconstructed portion
                                    cyclevehn.append(traj[j]['vehn'])
                                    vehtype.extend([traj[j]['vtype']]*totalpoints)
                             
                        # create data frame of cycle vehicle trajectories
                        trajc=pd.DataFrame({'speed':speed,'accel':accel,'vehtype':vehtype})
                                                
                        #%% generate MOVES inputs for this cycle                
                        linkID=c+1+s*100+inter_index*10000+direc_index*1000
                        opmSeconds=pd.DataFrame(index=opmModeSet,columns=cycles)
                        opmFraction=pd.DataFrame(index=opmModeSet,columns=cycles)
                        
                        ## links inputs
                        v_avg=np.mean(trajc['speed'])
                        q=len(cyclevehn)
                        # built links inputs                    
                        ls_linkID.append(linkID)
                        ls_linkVolume.append(q)
                        ls_linkAvgSpeed.append(v_avg)
                        ls_linkLength.append(segLength)
                        
                        #%% calculate link source type inputs and op mode
                        types=trajc['vehtype'].unique()
                        typecounts=trajc['vehtype'].value_counts()
                        typefraction=pd.Series(index=types)
                        for t in types:
                            # link source type inputs
                            typefraction[t]=typecounts[t]/sum(typecounts)
                            lst_linkID.append(linkID)
                            lst_sourceTypeID.append(typemap[t])
                            lst_fraction.append(typefraction[t])
                            
                            # op mode inputs
                            
                            # get parameters for calculating VSP
                            
                            sourceID=typemap[t]
                            params=paramdict[sourceID]
                            
                            trajct=trajc.loc[trajc['vehtype']==t] # filter by type
                            
                            ### assign vsp and speed bins
                            
                            # separate by brake/go
                            traj_br=trajct.loc[trajct['accel']<=brake_threshold]
                            # TODO: add 3-second brake rate condition
                            opmSeconds[c][0]=len(traj_br.index)
                            traj_go=trajct.loc[trajct['accel']>brake_threshold]
                    
                            # count idle mode vehicle seconds
                            traj_idle=traj_go.loc[(traj_go['speed']>=-m2k) & (traj_go['speed']<=m2k)]
                            opmSeconds[c][1]=len(traj_idle.index)
                    
                            ## accelerating/cruising modes/coasting modes
                    
                            traj_spd_vsp={}
                            for n in range(len(spdbins)-1):
                                # separate by speed bins
                                traj_spd=traj_go.loc[(traj_go['speed']>=spdbins[n]) & (traj_go['speed']<spdbins[n+1])]
                                
                                # get vsp for each second
                                traj_spd['vsp']=ts.getVSP(traj_spd['speed']*kph2mps,traj_spd['accel'],grade,params,sourceID)
                                
                                # separate by vsp bins
                                traj_spd_vsp[n]={}
                                for j in range(1,len(vspbins[n])):
                                    traj_spd_vsp[n][j]=traj_spd.loc[(traj_spd['vsp']>=vspbins[n][j-1]) & (traj_spd['vsp']<vspbins[n][j])]
                                    
                                # first and last vsp bin is open
                                traj_spd_vsp[n][0]=traj_spd.loc[traj_spd['vsp']<vspbins[n][0]]
                                traj_spd_vsp[n][len(vspbins[n])]=traj_spd.loc[traj_spd['vsp']>=vspbins[n][-1]]
                            
                            # last speed bin is open
                            
                            traj_spd=traj_go.loc[traj_go['speed']>=spdbins[-1]]
                            # get vsp for each second
                            traj_spd['vsp']=ts.getVSP(traj_spd['speed']*kph2mps,traj_spd['accel'],grade,params,sourceID)
                            
                            # separate by vsp bins
                            lastspd=len(spdbins)-1
                            traj_spd_vsp[lastspd]={}
                            for j in range(1,len(vspbins[lastspd])):
                                traj_spd_vsp[lastspd][j]=traj_spd.loc[(traj_spd['vsp']>=vspbins[lastspd][j-1]) & (traj_spd['vsp']<vspbins[lastspd][j])]
                                
                            # first and last vsp bin is open
                            traj_spd_vsp[lastspd][0]=traj_spd.loc[traj_spd['vsp']<vspbins[lastspd][0]]       
                            traj_spd_vsp[lastspd][len(vspbins[lastspd])]=traj_spd.loc[traj_spd['vsp']>=vspbins[lastspd][-1]]
                            
                            # fill in opm seconds with each op mode bin
                            n=0
                            mi=2
                            while n<=lastspd:
                                lastvsp=len(vspbins[n])
                                j=0
                                while j<=lastvsp:
                                    m=opmModeSet[mi]
                                    if not traj_spd_vsp[n][j].empty:
                                        opmSeconds[c][m]=len(traj_spd_vsp[n][j].index)
                                    else:
                                        opmSeconds[c][m]=0
                                    j+=1
                                    mi+=1
                                n+=1
                            
                            # get op mode fractions using op mode seconds
                            for m in opmModeSet:
                                if sum(opmSeconds[c])==0:
                                    print(str(i)+str(d)+' cycle'+str(c)+' type'+str(t)+'doesn\'t have any operating mode seconds.')
                                opmFraction[c][m]=opmSeconds[c][m]/sum(opmSeconds[c])
                
                            
                            # populate op mode counts
                            for pr in opmModeSetProcesses:
                                opm_type.extend([typemap[t]]*modes)
                                opm_process.extend([pr]*modes)
                                opm_link.extend([linkID]*modes)
                                opm_mode.extend(opmModeSet)
                                opm_fraction.extend(list(opmFraction[c]))
                            for pr in opmAllRunningProcesses:
                                opm_type.append(typemap[t])
                                opm_process.append(pr)
                                opm_link.append(linkID)
                                opm_mode.append(opmAllRunning)
                                opm_fraction.append(1)

#%%                                
                    time_id=time()-time_id0
                    print('Processed '+str(p)+' '+i+' '+d+' in %.fs' %time_id)
                    f.savefig(figoutpath+'/'+str(p)+str(i)+str(d)+'.png', bbox_inches='tight')
                    plt.show()
                direc_index+=1
            inter_index+=1
            
        ## put together outputs for exporting
    
        # Set up for result saving - create dataframes for MOVES inputs
        ls=ts.makeLinksTemplate()
        lst=ts.makeLinkSourceTypeTemplate()
        opm=ts.makeOpModeTemplate()
        
        # links
        ls['linkID']=ls_linkID
        lsSize=len(ls)
        ls['linkLength']=ls_linkLength
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
        
        # save path
        savepath=outpath+str(p)
        if not os.path.exists(savepath):
            os.makedirs(savepath)
            
        # export to csv
        ls.to_excel(savepath+'/links.xlsx',index=False)
        lst.to_excel(savepath+'/linkSourceTypes.xlsx',index=False)
        opm.to_excel(savepath+'/OpModeDistribution.xlsx',index=False)
    
print('runtime = ' +str(time()-start_t))
                        
                        
                        
                        
                        
                        