# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 12:00:35 2018

@author: Anjie
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from time import time

import toolSpace as ts

start_t=time()

#temp
inters=['Pinebush']
direcs=['NB']

# intersections, directions
#inters=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
#direcs=['NB','SB']

## start template for output file 
## (similar to fzp files so that the VissimOutputToMovesInput module can be used)
#sec=[]
#no=[]
#link=[]
#speed=[]
#accel=[]
#vehtype=[]

traj=pd.DataFrame(columns=['sec','no','link','lane','dist','pos','speed','accel','vehtype','sim'])

#import simulation data
lsas=ts.listFiles('lsa')
mers=ts.listFiles('mer')
rsrs=ts.listFiles('rsr')
fzps=ts.listFiles('fzp')

# runSpecs
runSpecs=json.load(open('runSpecs.json'))
periods=runSpecs['periods']

# modelSpecs
modelSpecs=json.load(open('modelSpecs-MajorOnly.json'))

# parameters for trajectory reconstruction
dt=1 # time step [s]
adj=.1 # parameter adjustment step for curve fitting
stoplag=1 # time between vehicle detection at stopbar and vehicle coming to a complete stop
startlag=1 # time between vehicle starting to accelerate and vehicle no longer being detected at stopped position
stoptime=1 # minimum time of detection for a vehicle that stopped at the detector position
t_tol=1 # acceptable difference between detected time and fitted curve time at the advance position
# original deceleration curve parameters from literature
k3o=0.005 
k4o=0.154
k5o=0.493 
# acceleration curve parameters from literature
beta0=3.369
beta1=0.072
# car-following parameters: time interval and effective vehicle length
leffs=pd.Series([5.5,13,11.5,13.5],index=[100,150,200,210]) # effective length [m]
hDis=2 # maximum headway for queue discharge
pstopcruise=.75 # if a vehicle stops, its speed is at most this fraction of its cruise speed


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

    for p in periods: # for each period (hour)
        rsrp=ts.filterbyperiod(rsr,p)
        fzpp=ts.filterbyperiod(fzp,p)
        
        inter_index=1 # for naming linkID
        for i in inters: # for each intersection
            tgi=tg.loc[tg['SC']==modelSpecs[i]['SignalController']] # filter to signal controller
            
            direc_index=1
            for d in direcs: # for each direction
                                
                # Dunbar northbound and Hwy401 southbound are out of scope (end intersections)
                if not ((i=='Dunbar' and d=='NB') or (i=='Hwy401' and d=='SB')):
                    time_id0=time()
                    f=ts.startFig('Lead Trajectories of '+str(p)+str(i)+str(d),'Time','Position',4)
                    
                    # segment model specs
                    segSpecs=modelSpecs[i][d]
                    
                    # links in the segment
                    intlink=segSpecs['links']['intersection']['thru']
                    approachlinks=segSpecs['links']['approach']['thru']
                    links=approachlinks+intlink
                                        
                    # link lengths
                    lengthMap=ts.getLengthMap(links)
                    
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
                    fzp_vehicles=fzppid['no'].unique() # all vehicles numbers
                    for n in fzp_vehicles:
                        veh=fzppid.loc[fzp['no']==n]
                        plt.plot(veh['sec'],veh['abspos'],'grey')

                    # filter data to segment
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
                    
                    # get cruising speed from travel times
                    ttcruise=np.percentile(rsrpid['Trav'],5) # use 5th percentile of travel times
                    segLength=sum(lengthMap) # length of segment
                    v0=segLength/ttcruise

                    # get detector numbers
                    detn={}
                    segDetSpecs=segSpecs['detectors']
                    detn['adv']=segDetSpecs['advance']
                    detn['stop']=segDetSpecs['stopbar']
                    detn['end']=segDetSpecs['end']

                    
                    xdet={} # get positions of detectors
                    for sd in detn:
                        xdet[sd]=int(merDetail.loc[merDetail[3]==detn[sd][0],9].iloc[0])
                    # for stopbar detector
                    posOffset=sum(lengthMap.loc[pd.Series(approachlinks).iloc[0:-1]]) # offset distance is length of upstream links
                    xdet['stop']=xdet['stop']+posOffset
                    # for advance detector
                    detlink=merDetail.loc[merDetail[3]==detn['adv'][0],5].iloc[0] # link detector is on
                    detlinkindex=links.index(detlink)
                    posOffset=sum(lengthMap.loc[pd.Series(approachlinks).iloc[0:detlinkindex]]) # offset distance is length of upstream links up to link of detector
                    xdet['adv']=xdet['adv']+posOffset
                    # end of intersection
                    xdet['end']=xdet['stop']+float(lengthMap[intlink]) #add on intersection span

                    
                    # iterate through cycles in this period
                    cycles=range(len(tgidp)-1) # the added on green start is not counted - just for closing the bin               
                    
                    for c in cycles:
                        no=[]
                        link=[]
                        speed=[]
                        accel=[]
                        vehtype=[]
                        
                        lanes=len(detn['stop'])
                        for l in range(lanes):
                            
                            # detectors for this lane
                            detectorStop=detn['stop'][l]

                            # fine vehicles detected in this lane during this cycle at the stopbar
                            merExid_stopbar=merEx.loc[merEx['Measurem.']==detectorStop] # stopbar exit detections
                            merExid_stopbar_cycle=merExid_stopbar.loc[(merExid_stopbar['t(Exit)']>=tgidp[c])&(merExid_stopbar['t(Exit)']<=tgidp[c+1])]
                            vehicles=merExid_stopbar_cycle['VehNo'].unique()
                            
                            if len(vehicles)>0:
                            
                                # filter detections to vehicles 
                                # (so that adv detections before the cycle starts are included as long as the vehicle passes the stopbar during the cycle.)
                                merEnc=merEn.loc[merEn['VehNo'].isin(vehicles)]
                                merExc=merEx.loc[merEx['VehNo'].isin(vehicles)]
                                
                                # data for stopbar detection of this lane
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
                                    stopt1.append(stopt1n)
                                    stopt2.append(stopt2n)
                                    advt1.append(advt1n)
                                    advt2.append(advt2n)
                                    plt.plot([stopt1n,stopt2n],[xdet['stop']]*2,'orange')
                                    plt.plot([advt1n,advt2n],[xdet['adv']]*2,'orange')
                                    

                                # dataframe of detector data
                                detection=pd.DataFrame({'vehn':vehn,'advt1':advt1,'advt2':advt2,'stopt1':stopt1,'stopt2':stopt2,'vtype':vtype})

                                detection=detection.sort_values(by='stopt1') # sort by time of reaching stopbar
                                detection=detection.reset_index(drop=True)
                                
                                leadTraj={} # trajectory of lead vehicle
                                
                                leaddetection=detection.iloc[0]
                                stopbarWait=leaddetection['stopt2']-leaddetection['stopt1']
                                
                                # check if there was a queue - how long the first vehicle stopped at the stopbar
                                if stopbarWait>stoptime: # if the vehicle had stopped, construct deceleration and acceleration curves for the lead vehicle
                                    
                                    tdelta=9999 # time difference for fitting curve to detector data (to be minimized)
                                    
                                    # deceleration curve
                                    
                                    # start with original values
                                    k3=k3o
                                    k4=k4o
                                    k5=k5o
                                    while tdelta>t_tol: # time difference between advance detection and extrapolated cruise trajectory at advance
                                        v=[v0] # speed profile
                                        dc=[0] # deceleration profile
                                        xdec=[0] # positions 
                                        tdec=[0] # time
                                        while v[-1]>0: # until speed has reached zero,
                                            dc.append(-k3*v[-1]**2+k4*v[-1]+k5) # using deceleration curve from literature
                                            v.append(v[-1]-(dc[-1]+dc[-2])/2*dt) # get speed for each time step
                                            xdec.append(xdec[-1]+v[-1]*dt+(dc[-1]+dc[-2])/4*dt**2) # append to positions
                                            tdec.append(tdec[-1]+dt) # append to times
                                
                                        # shift deceleration trajectory to match end of deceleration with beginning of stopped trajectory
                                        tdeci=pd.Series(tdec)+leaddetection['stopt1']-tdec[-1]+stoplag
                                        xdeci=pd.Series(xdec)+xdet['stop']-xdec[-1]
                                        tcruise=(xdeci[0]-xdet['adv'])/v[0] # time spent cruising after advance detection before deceleration
                                        textra=tdeci[0]-tcruise # extrapolated time of reaching advance detector from beginning of deceleration
                                        tdelta=leaddetection['advt1']-textra # difference between advance detection time and extrapolated time
                                        if xdeci[0]<xdet['adv']: # if deceleration started before advance detection
                                            tdelta=leaddetection['advt1']-tdeci[abs(xdeci-xdet['adv']).idxmin()]
                                        if tdelta>0: # if curve is earlier than detection
                                            k3=k3-k3o*adj 
                                            k4=k4+k4o*adj
                                            k5=k5+k5o*adj
                                        else: # if curve is later than detection
                                            k3=k3+k3o*adj 
                                            k4=k4-k4o*adj
                                            k5=k5-k5o*adj
                                    
                                    vdec=v
                                    adec=dc
                                    
                                    # stopped curve
                                    tStopEnd=leaddetection['stopt2']
                                    
                                    # acceleration curve
                                    v=[0] # start at zero speed
                                    tdelta=9999 # to be minimized
                                    a=[0] # acceleration starts at 0
                                    
                                    xacc=[0]
                                    tacc=[0]
                                    while v[-1]<v0:
                                        a.append(beta0-beta1*v[-1])
                                        v.append(v[-1]+(a[-1]+a[-2])/2*dt)
                                        xacc.append(xacc[-1]+v[-1]*dt+(dc[-1]+dc[-2])/4*dt**2)
                                        tacc.append(tacc[-1]+dt)
                                    
                                    # shift acceleration curve to start at the end of the stopped curve
                                    tacc+=tStopEnd-startlag
                                    xacc=[xx+xdet['stop'] for xx in xacc]
                                    aacc=a
                                    vacc=v
                                    
                                    
                                    # compile curves into lead trajectory
                                    leadTraj['t']=list(tdeci)+list(tacc)
                                    leadTraj['x']=list(xdeci)+list(xacc)
                                    leadTraj['a']=list(adec)+list(aacc)
                                    leadTraj['v']=list(vdec)+list(vacc)
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
                                        h=detection.loc[j,'stopt1']-detection.loc[j-1,'stopt1'] # headway wrt preceding vehicle at stopbar 
                                        leff=leffs[detection.loc[j-1,'vtype']] # get effective length of preceding vehicle based on vehicle type
                                        
                                        # find T (time interval) for Newell's model - time lag of trajectories
                                        t2=detection.loc[j,'stopt2']
                                        xlead=xdet['stop']+leff # effective position of preceding vehicle corresponding the following vehicle at the stopbar
                                        index=abs(np.array(traj[j-1]['x'])-xlead).argmin() # find point in lead closest to stopbar
                                        t1=traj[j-1]['t'][index] # time when preceding vehicle is at this position
                                        if xlead>traj[j-1]['x'][-1]: # if lead acceleration curve ends before effective lead position,
                                            t1=traj[j-1]['t'][-1]+(xdet['stop']+leff-traj[j-1]['x'][-1])/v0 # then extrapolate using cruise speed
                                        tlag=t2-t1 # time interval for Newell's model
                                        
                                        if vapproach<pstopcruise*v0: #check if vehicle stopped in queue
                                            # offset end of stopped curve of preceding vehicle to get end of stopped curve for current following vehicle
                                            tgo+=tlag
                                            xstop-=leff
                                            # shift deceleration curve spatially so that it ends where the stopped curve starts
                                            xdeci=pd.Series(xdec).add(xstop-xdec[-1])
                                            # shift deceleration curve temporally
                                            if xdeci[0]>xdet['adv']: # if deceleration curve starts after the advance location
                                                tdeci=pd.Series(tdec).add(detection.loc[j,'advt1']+(xdeci[0]-xdet['adv'])/v0) # match extrapolated upstream position to advance detector
#                                                tdeci=pd.Series(detection.loc[j,'advt1']).append(tdeci)
#                                                xdeci=pd.Series(xdet['adv']).append(xdeci)
                                            else: # if deceleration curve starts before the advance location
                                                tdelta=tdec[abs(np.array(xdeci-xdet['adv'])).argmin()]-tdec[0] # find the point closest to the advance detector
                                                tdeci=pd.Series(tdec).add(detection.loc[j,'advt1']-tdelta) # shift deceleration curve to match this point to advance detector
                                            
                                            #acceleration curve
                                            v=[0] # speed vector
                                            a=[0] # acceleration vector
                                          
                                            xacc=[0] # position vector
                                            tacc=[0] # time vector
                                            tdelta=9999 # to be minimized

                                            while v[-1]<v0: # until speed reaches cruise speed
                                                a.append(beta0-beta1*v[-1]) # accleration of this time step, from literature
                                                v.append(v[-1]+(a[-1]+a[-2])/2*dt) # append speed for each time step
                                                xacc.append(xacc[-1]+v[-1]*dt+(a[-1]+dc[-2])/4*dt**2) # append to position vector
                                                tacc.append(tacc[-1]+dt) # append to time vector
                                            vacc=v
                                            aacc=a
                                            tacc+=tgo-startlag # shift start of acceleration to end of stopped curve, account for start lag
                                            xacc=[xx+xstop for xx in xacc]
                                            
                                            # put together decelerationg and acceleration curves
                                            traj[j]['t']=list(tdeci)+list(tacc)
                                            traj[j]['x']=list(xdeci)+list(xacc)
                                            traj[j]['v']=list(vdec)+list(vacc)
                                            traj[j]['a']=list(adec)+list(aacc)
                                            
                                                                                    
                                        elif h<hDis: #didn't stop but followed closely
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
                                                                                               
                                        else: #free flow vehicles
                                            # straight trajectories from advance detector to stop detector
                                            traj[j]={'t':[t0],'x':[xdet['adv']],'a':[0],'v':[v0],'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']}
                                            ti=detection.loc[j,'advt1']
                                            tf=detection.loc[j,'stopt1']
                                            tflow=tf-ti
                                            xi=xdet['adv']
                                            xf=xdet['stop']
                                            xflow=xf-xi
                                            vflow=xflow/tflow
                                            traj[j]['t']=np.arange(ti,ti+int(tflow),1)
                                            points=len(traj[j]['t'])
                                            traj[j]['x']=np.linspace(xi,xf,points)
                                            traj[j]['v']=[v0]*points
                                            traj[j]['a']=[0]*points
                                            
                                        plt.plot(traj[j]['t'],traj[j]['x'],'g--')
                               
                                else: # if first vehicle didn't stop, then trajectory is assumed to be straight from advance to stopbar
                                    # lead trajectory
                                    ti=detection.loc[0,'advt1']
                                    tf=detection.loc[0,'stopt1']
                                    tflow=tf-ti
                                    xi=xdet['adv']
                                    xf=xdet['stop']
                                    xflow=xf-xi
                                    vflow=xflow/tflow
                                    leadTraj['t']=np.arange(ti,ti+int(tflow),1)
                                    points=len(leadTraj['t'])
                                    leadTraj['x']=np.linspace(xi,xf,points)
                                    leadTraj['v']=[v0]*points
                                    leadTraj['a']=[0]*points
                                    leadTraj['vehn']=detection.loc[0,'vehn']
                                    leadTraj['vtype']=detection.loc[0,'vtype']
                                    
                                    traj={0:{'t':leadTraj['t'],'x':leadTraj['x'],'a':leadTraj['a'],'v':leadTraj['v'],'vehn':leadTraj['vehn'],'vtype':leadTraj['vtype']}}
                                    for j in detection.index[1:]:
                                        traj[j]={'vehn':detection.loc[j,'vehn'],'vtype':detection.loc[j,'vtype']} # fill in vehicle number and type first
                                        ti=detection.loc[j,'advt1']
                                        tf=detection.loc[j,'stopt1']
                                        tflow=tf-ti
                                        xi=xdet['adv']
                                        xf=xdet['stop']
                                        xflow=xf-xi
                                        vflow=xflow/tflow
                                        traj[j]['t']=np.arange(ti,ti+int(tflow),1)
                                        points=len(traj[j]['t'])
                                        traj[j]['x']=np.linspace(xi,xf,points)
                                        traj[j]['v']=[v0]*points
                                        traj[j]['a']=[0]*points
                                
                                
                                
                                # plot lead trajectories
                                
                                plt.plot(leadTraj['t'],leadTraj['x'],'r--')
                                
                                # complete reconstructed trajectories for lane and add to cycle collection
                                for j in traj:
                                    # filling for before the reconstructed trajectory (constant speed)
                                    xfill=traj[j]['x'][0]
                                    tfill=xfill/v0
                                    fillpoints=int(tfill) # number points, i.e. seconds, to fill in before the trajectory
                                    vfillList1=[v0]*fillpoints # array of speeds to fill in
                                    afillList1=[0]*fillpoints # array of accelerations to fill in 

                                    # filling for after the reconstructed trajectory (constant speed)
                                    xfill=xdet['end']-traj[j]['x'][-1]
                                    tfill=xfill/v0
                                    fillpoints=int(tfill) # number points, i.e. seconds, to fill in before the trajectory
                                    vfillList2=[v0]*fillpoints # array of speeds to fill in
                                    afillList2=[0]*fillpoints # array of accelerations to fill in 
                                    
                                    totalpoints=len(vfillList1)+len(traj[j]['v'])+len(vfillList2) # number of points in whole trajectory
                                    
                                    # append to table of trajectory records for the cycle
                                    speed.extend(vfillList1+traj[j]['v']+vfillList2) # join trajectories before, during, and after reconstructed portion
                                    accel.extend(afillList1+traj[j]['v']+afillList2) # join trajectories before, during, and after reconstructed portion
                                    no.append(traj[j]['vehn'])
                                    vehtype.extend([traj[j]['vtype']]*totalpoints)
                                    
                                    
                                    
                                    
#                                    linkID=c+1+s*100+inter_index*10000+direc_index*1000
#                                    link.extend([linkID]*)
                                
                    time_id=time()-time_id0
                    print('Processed '+str(p)+' '+i+' '+d+' in %.fs' %time_id)
                    f.savefig('TrajReconstruction-TestPlots/'+str(p)+str(i)+str(d)+'.png', bbox_inches='tight')
                direc_index+=1
            inter_index+=1
print('runtime = ' +str(time()-start_t))
                        
                        
                        
                        
                        
                        