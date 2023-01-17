# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 16:37:15 2022

@author: daniel
"""
import numpy as np
import pandas as pd
import os

# set working directory
os.chdir('C:/Users/danie/OneDrive - Durham University/Daniel_PhD_projects/CAA_II')

means = pd.read_csv("weather_data/weather_data.csv")

days = [31,28,31,30,31,30,31,31,30,31,30,31]

years = [i for i in range(1997, 2021)]

dailies = []

for i in range(len(years)):   
    df = pd.DataFrame({'id':means.id, 'year':years[i], 'month':means.month, 'day':means.day, 'daytime':means.daytime, 'dir':means.dir, 'speed':means.speed, 'cloud':means.cloud, 'temp':means.temp})
    dailies.append(df)

daily_full = pd.concat(dailies).set_index(['id','year','month','day'])

ghcnd = []

for i in range(461):
   
    try:       
       
        df = pd.read_csv('ghcnd_us/' + np.unique(means.id)[i] + '.csv', usecols = ['STATION','LATITUDE','LONGITUDE','DATE','TMAX','TMIN'])
        df.rename({'STATION':'id'}, inplace=True, axis=1)
        df['year'] = df['DATE'].str[0:4].astype(int)
        df['month'] = df['DATE'].str[5:7].astype(int)
        df['day'] = df['DATE'].str[8:10].astype(int)
        df = df[(df.year >= 1997) & (df.year <= 2020)]
        df['TMAX'] = df.TMAX/10 + 273.
        df['TMIN'] = df.TMIN/10 + 273.
        df['tday'] = df.TMAX *.75 + df.TMIN *.25
        df['tnight'] = df.TMAX *.25 + df.TMIN *.75
               
        df = df.loc[(df.tday.notnull()) & (df.tnight.notnull()), ['id','year','month','day','tday','tnight']]
        
        ghcnd.append(df)
                    
    except:
       
        pass
 
ghcnd = pd.concat(ghcnd).set_index(['id','year','month','day'])

daily_full.loc[daily_full.index.isin(ghcnd.index) & (daily_full.daytime == 1), 'temp'] = ghcnd.tday
daily_full.loc[daily_full.index.isin(ghcnd.index) & (daily_full.daytime == 0), 'temp'] = ghcnd.tnight


daily_full.to_csv('./weather_data/daily_weather.csv', index=False)

print('section completed')