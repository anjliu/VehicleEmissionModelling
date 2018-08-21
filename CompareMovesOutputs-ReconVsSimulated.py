# -*- coding: utf-8 -*-
"""
Created on Sat Jul 28 13:44:47 2018

@author: Anjie
"""

import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import matplotlib.pyplot as plt
import timeit
import numpy as np

import toolSpace as ts

t0=timeit.default_timer()

# decoded MOVES output scenarios
scenarios=['TrajectoryReconstruction-wIdle','thruMovements']

# decoded MOVES output path
mo_decodedpath='MOVESoutputsDecoded/'

# data of interest
hours=['7-8']#,'8-9']
intersections=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
types=['Combination Truck','Single Unit Truck','Transit Bus','Passenger Car']

pollutants=['Atmospheric CO2','Nitrous Oxide','Methane','Total Gaseous Hydrocarbons','Total Energy']
approachDirections=['NB','SB']
fuelTypes=['Diesel','Gasoline']
modelYears=range(1986,2017,1)

# apply desired filters
filters={
        'hour':hours,
        'intersection':intersections,
        'approachDirection':approachDirections,
        'cycle':range(99),
        'sourceType':[types[-1]],
        'pollutant':[pollutants[0]],
        'fuelType':fuelTypes,
        'modelYear':[2013]
        }

# get data and filter
scenariosdf=[]
for n in scenarios:
    scen_df=pd.read_csv(mo_decodedpath+n+'.csv')
    scen_df=scen_df.dropna()
    for f in filters:
        scen_df=scen_df.loc[scen_df[f].isin(filters[f])]
    scenariosdf.append(scen_df)

mergecols=['link','pollutant','process','sourceType','fuelType','modelYear','hour','approachDirection','intersection','cycle']

# load qlengths
qlengthdf=pd.read_csv('metrics.csv')

# merge data sets
mergeddf=pd.merge(scenariosdf[0],scenariosdf[1],on=mergecols,suffixes=['-recon','-simul'])
mergeddf=mergeddf.loc[mergeddf['emissionQuant-simul']!=0] # remove zero emissions records as they interfere with relative error calculations
mergeddf=pd.merge(mergeddf,qlengthdf,on=['link','hour'])

# remove first and last cycles
newdf=pd.DataFrame()
for i in intersections:
    for d in approachDirections:
        dfid=mergeddf.loc[(mergeddf['intersection']==i)&(mergeddf['approachDirection']==d)]
        if not dfid.empty:
            lastcycle=max(dfid['cycle'])
            firstcycle=min(dfid['cycle'])
            dfid=dfid.loc[(dfid['cycle']!=lastcycle)&(dfid['cycle']!=firstcycle)]
            newdf=newdf.append(dfid)
mergeddf=newdf

f=ts.startFig('Comparison: '+str(filters['intersection'])+' '+str(filters['pollutant'])+' '+str(filters['sourceType']),
              'Emission Estimates from Simulated Trajectories [g]',
              'Emission Estimates from Reconstructed Trajectories [g]',
              2)

# for all cycles:
emi_recon=mergeddf['emissionQuant-recon']
emi_simul=mergeddf['emissionQuant-simul']

plt.plot(emi_simul,emi_recon,'og',label='Cycles with queue length > advance detector setback')

difference=emi_recon-emi_simul
errors=difference/emi_simul
errors_abs=[abs(x) for x in errors]
mape=np.mean(errors_abs)
print('MAPE = %.1f' %(mape*100))

# for cycles with qlength < x:
ql_cutoff=51 # [m] - exclude cycles where qlength is greater than the cutoff
mergeddf2=mergeddf.loc[mergeddf['qlength']<ql_cutoff]
emi_recon2=mergeddf2['emissionQuant-recon']
emi_simul2=mergeddf2['emissionQuant-simul']

print('Short queues only:')
plt.plot(emi_simul2,emi_recon2,'sk',label='Cycles with queue length < advance detector setback')

difference2=emi_recon2-emi_simul2
errors2=difference2/emi_simul2
errors_abs2=[abs(x) for x in errors]
mape2=np.mean(errors_abs)
print('MAPE = %.1f' %(mape2*100))

mergeddf['errors_relative']=errors
plt.legend(prop={'size': 15})


# errors vs...
f2=ts.startFig('Relative Errors for Different Levels of Congestion',
               'Level of Congestion (Volume/Capacity) [-]',
               'Relative Error',
               2)
plt.plot(mergeddf['congestion'],errors,'o')

# errors vs...
f3=ts.startFig('Relative Errors for Different Levels of Emissions',
               'Emissions [g]',
               'Relative Error',
               2)
plt.plot(emi_simul,errors,'o')

# errors vs...
f4=ts.startFig('Relative Errors for Different Queue Lengths',
               'Emissions [g]',
               'Relative Error',
               2)
plt.plot(mergeddf['qlength'],errors,'o')

# errors vs...
f4=ts.startFig('Relative Errors for Different Volumes',
               'Number of Vehicles in the Cycle',
               'Relative Error',
               2)
plt.plot(mergeddf['nvehicles'],errors,'o')

# Absolute error vs...
f4=ts.startFig('Absolute Error for Levels of Congestion',
               'Level of Congestion (Volume/Capacity)',
               'Absolute Error [g]',
               2)
plt.plot(mergeddf['congestion'],difference,'o')




