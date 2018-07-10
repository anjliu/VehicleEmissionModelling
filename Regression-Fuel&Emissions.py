# -*- coding: utf-8 -*-
"""
Created on Thu May 31 01:37:59 2018

@author: Anjie
"""

import pandas as pd
import seaborn as sns

import toolSpace as ts

# data of interest
periods=['morning1','morning2','morning3']
intersections=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
types=[21,42,52,61]
pollutants=['Atmospheric CO2','Nitrous Oxide','Methane','Total Gaseous Hydrocarbons','Total Energy']
pol=pollutants[4]

filters={
        'period':periods,
        'inter':intersections,
        'direc':['SB','NB'],
        'cycle':range(999),
        'type':[types[0]]
        }

# get data and filter
fueldf=pd.read_csv('FuelEstimate-VissimOutput-thru.csv')
emissiondf=pd.read_csv('MovesOutputDecoded-thru.csv')
emissionpol=emissiondf.loc[emissiondf['pollutant']==pol]

cols=['period','inter','direc','cycle','type']
# merge data sets
fuelemission=pd.merge(emissionpol,fueldf,on=cols)

# filter further
for x in cols:    
    fuelemission=fuelemission.loc[fuelemission[x].isin(filters[x])]

x=fuelemission['fuel']
y=fuelemission['emi']
f=ts.startFig('Fuel Consumption and Emissions Regression','Fuel Consumed by Cycle and Vehicle Type [mL]','Emissions [g]',1)
sns.regplot(x,y).set(xlabel='Fuel Consumed by Cycle and Vehicle Type [mL]',ylabel='Emissions [g]')


fuelemission.to_csv('fuel&emissions.csv')