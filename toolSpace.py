# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 13:42:42 2018

@author: Anjie
"""
import json, os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.style.use('seaborn-whitegrid')

g=9.81 # gravitational constant [m/s^2]
# get run specs
runSpecs=json.load(open('runSpecs.json'))
dt=runSpecs['timeInterval'] # time interval for aggregation [s]
wt=runSpecs['warmUpTime'] # warm-up time: time at beginning of simulation to exclude [s]
durt=runSpecs['durationOfInterest'] # duration of data to include starting from warm-up time [s]

# fuel formula constants
fi=1200
b1=.16
b2=.000518
beta=.1
pmax=120
gr=0 # tentative

periods=runSpecs['periods']


def startFig(title,xlabel,ylabel,size): # start a standard figure. Size:{1,2}
    pltSpec=json.load(open('plotSpecs.json'))
    fL=pltSpec['figL'+str(size)]
    fW=pltSpec['figW'+str(size)]
    fs1=pltSpec['fontsize']['main']
    fs2=pltSpec['fontsize']['sub']
    f=plt.figure(figsize=(fL,fW))
    f.suptitle(title,fontsize=fs1)
    plt.xlabel(xlabel,fontsize=fs2)
    plt.ylabel(ylabel,fontsize=fs2)
    plt.tick_params(labelsize=16)
    return f

def getLengthMap(links): # get a map of link lengths for given links
    #import link lengths
    linkProperties=pd.read_excel('LinkProperties.xlsx',header=1)
    #create table of links and their lengths
    linklengths=[]
    for l in links:
        linklengths.append(linkProperties.loc[linkProperties['$LINK:NO']==l,'LENGTH2D'].values[0])
    lengthMap=pd.Series(linklengths,index=links)
    return lengthMap

def makeLinksTemplate(): # start a Links data frame
    ls=pd.DataFrame(columns=['linkID','countyID','zoneID','roadTypeID','linkLength','linkVolume','linkAvgSpeed','linkDescription','linkAvgGrade'])
    return ls

def makeLinkSourceTypeTemplate(): # start a Link Source Type data frame
    lst=pd.DataFrame(columns=['linkID','sourceTypeID','sourceTypeHourFraction'])
    return lst

def makeOpModeTemplate(): # start a Operating Mode Distribution data frame
    opm=pd.DataFrame(columns=['sourceTypeID','hourDayID','linkID','polProcessID','opModeID','opModeFraction'])
    return opm

def getVSPparams(sourceID,vspParams): # get VSP parameters given vehicle type
    cA=vspParams.loc[vspParams['ID']==sourceID,'cA'].iloc[0]
    cB=vspParams.loc[vspParams['ID']==sourceID,'cB'].iloc[0]
    cC=vspParams.loc[vspParams['ID']==sourceID,'cC'].iloc[0]
    m=vspParams.loc[vspParams['ID']==sourceID,'mass'].iloc[0]
    f=vspParams.loc[vspParams['ID']==sourceID,'massFactor'].iloc[0]
    params=[cA,cB,cC,m,f]
    return params
    
def getVSP(v,a,gr,params,sourceID): # calculate VSP given speed, acceleration, grade, parameters
    cA,cB,cC,m,f=params
    if sourceID==21:
        vsp=cA/m*v+cB/m*v**2+cC/m*v**3+v*(a+g*np.sin(gr))
    else:
        vsp=(cA*v+cB*v**2+cC*v**3+m*v*(a+g*np.sin(gr)))/f
    for x in vsp:
        if x<0:
            x=0
    return vsp

def filterbyperiod(df,period): # filter to period of day
    df=df.loc[(df.iloc[:,0]>=periods[period]['start']) & (df.iloc[:,0]<=periods[period]['end'])]
    return df

#%% importing Vissim output files
def import_fzp(path): # import vehicle records file
    fzp=pd.read_csv(path,sep=';',skiprows=range(0,23),usecols=list(range(8))+[10,11],
                    names=['sec','no','link','lane','dist','pos','speed','accel','vehtype','sim'],
                    dtype={'sec':float,'no':int,'link':int,'lane':int,'dist':float,'pos':float,'speed':float,
                           'accel':float,'vehtype':int,'sim':int})
    return fzp

def import_lsa(path): # import signal phasing file
    lsa=pd.read_csv(path,index_col=False,skiprows=127,sep=';',
                    names=['tsim','tcycle','SC','SG','newState','lastDuration','type','SGcause'],
                    dtype={'tsim':float,'newState':str})
    return lsa

def import_rsr(path): # import travel time measurements
    rsr=pd.read_csv(path,sep=';',index_col=False,header=5)
    rsr=rsr.rename(columns=lambda x: x.strip())
    return rsr

def import_mer(path): # import data point measurements, i.e. simulated detector data
    drow=6
    hrow=171 #row number for headers where detector data starts
    
    # reads the first section of data: data collection point number, lane & position
    merDetail=pd.read_csv(path,nrows=hrow-drow-1,skiprows=9,header=None,delim_whitespace=True)
    merDetail[3]=merDetail[3].replace({':':''}, regex=True)
    merDetail=merDetail.astype({3:int})
    
    # reads the remaining section: the data collection point records
    mer=pd.read_csv(path,header=hrow,sep=';')
    mer=mer.rename(columns=lambda x: x.strip()) #strip whitespace from headers
    mer=mer.astype({'Measurem.':int,'t(Entry)':float,'t(Exit)':float,'VehNo':int,'Vehicle type':int})
    return mer, merDetail
#%%

def remapLinkNo(links): # replace link numbers with values starting from 1
    linksunique=np.array(links).unique()
    newmap=pd.Series(range(1,len(linksunique)+1,index=linksunique))
    return newmap
 
def listFiles(extension):
    curdir=os.getcwd()
    vispath=curdir+'/VissimOutputs'
    filelist=[]
    for file in os.listdir(vispath):
        if file.endswith('.'+extension):
            filelist.append(vispath+'/'+file)
    return filelist

def importMovesOutput(path):
    mo=pd.read_excel(path)
    return mo
    
def getSegmentLinks(interspecs,direc,mov):
    links=interspecs[direc]['links']  
    allLinks=[] # list links for intersection and approach
    if mov=='all':
        for section in links:
            for movement in links[section]:
                allLinks.extend(links[section][movement])
    else:
        for section in links:
            allLinks.extend(links[section]['thru'])
    return allLinks

def getIntLinks(interspecs,direc,mov):
    links=interspecs[direc]['links']['intersection']
    intLinks=[] # list links for intersection and approach
    if mov=='all':
        for movement in links:
            intLinks.extend(links[movement])
    else:
        intLinks=links['thru']
    return intLinks  
    
# superceded
#def findGreenPhases(lsa):
#    tg=lsa.loc[lsa['newState'].str.contains('green'),'tsim'] # green start times
#    lsa=lsa.loc[lsa['tsim']>=tg.iloc[0]] # filter to after first green start
#    tr=lsa.loc[lsa['newState'].str.contains('red'),'tsim'] # red start times
#    tg=tg.loc[tg<=tr.iloc[-1]] # remove green starts without an end
#    tg=tg.reset_index(drop=True)
#    tr=tr.reset_index(drop=True)
#    return tg,tr
    
def getFuel(fzp,m): # requires fzp data and vehicle mass
    v=fzp['speed']/3.6
    a=fzp['accel']
    pc=b1*v+b2*v**3
    pii=m*a*v/1000
    pg=9.81*m*(gr)*v/1000
    alpha=fi/3600
    pt2=pc+pii+pg
    pt=[min(pmax,x) for x in pt2]
    ft=[alpha+beta*x if x>0 else alpha for x in pt]
    f=sum(ft)
    return f

def getVehTypeMap(): # vehicle type ID mapping from Vissim to MOVES
    vehtypefile=pd.read_csv('VehicleTypeMapping.csv')
    typemap=pd.Series(vehtypefile['MOVEStype'].values,index=vehtypefile['VissimType'])
    return typemap
