# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 12:15:35 2018

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

# data of interest
periods=['morning1','morning2','morning3','afternoon1','afternoon2']
intersections=['Dunbar','Bishop','Sheldon','Pinebush','Hwy401']
types=[21,42,52,61]
pollutants=['Atmospheric CO2','Nitrous Oxide','Methane','Total Gaseous Hydrocarbons','Total Energy']

# apply desired filters
pol=pollutants[2]
filters={
        'period':periods,
        'inter':intersections,
        'direc':['NB','SB'],
        'cycle':range(999),
        'type':[types[0]]
        }

# get data and filter
vspdf=pd.read_csv('VSP-VissimOutput-thru.csv')
vspdf=vspdf.dropna()
emissiondf=pd.read_csv('MovesOutputDecoded-thru.csv')
emissionpol=emissiondf.loc[emissiondf['pollutant']==pol]

mergecols=['period','inter','direc','cycle','type']
# merge data sets
vspemission=pd.merge(emissionpol,vspdf,on=mergecols)

# filter further
for x in mergecols:    
    vspemission=vspemission.loc[vspemission[x].isin(filters[x])]

vspemission['label']=''
for i in vspemission.index:
    row=vspemission.loc[i]
    vspemission['label'][i]=str(row['period'])+' '+str(row['inter'])+' '+str(row['direc'])+' cycle'+str(row['cycle'])

# variables for regression
x11=vspemission['vspslow']
x12=vspemission['vspfast']
x1=x11+x12
x2=vspemission['vspmass']
x3=vspemission['avgspeed']
x4=vspemission['avgaccel']
x51=vspemission['vehicleslow']
x52=vspemission['vehiclefast']
x5=x51+x52
x6=vspemission['avgTT']
y=vspemission['emi']

# regression
xvars=[x11,x12,x51,x52] # select a set of independent variables
X=pd.DataFrame(xvars).transpose() 
X=sm.add_constant(X) # add a constant
regres=sm.OLS(y,X).fit()
results=regres.summary()
coef=regres.params

print(results)

length=len(y)

def predict(coef,xvars,length):
    y=np.array([coef[0]]*length)
    i=0
    for c in coef[1:]:
        y+=c*xvars[i]
        i+=1
    return y

# predicted values and errors
vspemission['predicted']=predict(coef,xvars,length)
predicted=vspemission['predicted']
error=vspemission['predicted']-y
errorPercent=abs(error)/y
mape=sum(errorPercent)/len(y)

vspemission['error']=error
vspemission=vspemission.sort_values('predicted')
vspemission=vspemission.reset_index(drop=True)

# plot figures
f1=ts.startFig('Predicted vs Simulated Emissions - '+str(pol),'Actual','Predicted',1)
ax=sns.regplot(vspemission['emi'],vspemission['predicted'])
ax.set(xlabel='Simulated Emissions [g]', ylabel='Predicted Emissions [g]')
plt.show()
#
#f2=ts.startFig('Predicted vs Simulated Emissions','Sample','Emission',1)
#plt.plot(vspemission['predicted'],'o',label='Predicted Emissions')
#plt.plot(vspemission['emi'],'o',label='Actual Emissions')
#plt.bar(vspemission.index,vspemission['error'],bottom=vspemission['emi'],width=.05)
#plt.legend()
#
#f3=ts.startFig('Predicted vs Simulated Emissions','','Percent Errors',1)
#plt.plot(errorPercent.reset_index(drop=True),'o')

print('MAPE=%.2f' %mape)

vspemission.to_csv('VSP&emissions.csv')

runtime=timeit.default_timer()-t0
print('Runtime: %.2fs' %runtime)

