# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 18:17:23 2018

@author: Anjie
"""

import pandas as pd

ratesdf=pd.read_csv('rates-fromMOVES.csv')
ratesdf=ratesdf.loc[ratesdf['matchlink']!='FALSE']
ratesdf['typeModePolProcess']=ratesdf['typeModePolProcess'].astype(int)

rates=pd.Series(ratesdf['emissionQuant'])
rates=rates.rename(index=ratesdf['typeModePolProcess'])

rates.to_csv('rates.csv')