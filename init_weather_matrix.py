# -*- coding: utf-8 -*-
"""
Created on Sun May 15 22:38:43 2022

@author: danie
"""

import os
import numpy as np
import pandas as pd

files=os.listdir('./hourly_normals/')
for i in range(len(files)):
    files[i] = './hourly_normals/'+files[i]

weather_data = []
for i in range(len(files)):
    # load hourly US wind data from NOAA:
    df = pd.read_csv(files[i], header=0)
    df = df[['STATION', 'LATITUDE', 'LONGITUDE', 'month', 'day', 'hour', 'HLY-WIND-VCTDIR', 'HLY-WIND-AVGSPD', 'HLY-CLOD-PCTCLR', 'HLY-CLOD-PCTFEW', 'HLY-TEMP-NORMAL']]
    # variable names: id | latitude | longitude | month | day | hour | wind direction | wind speed | % clear sky | % few clouds | temperature (degrees F)
    df.columns = ['id', 'lat', 'lon', 'month', 'day', 'hour', 'dir', 'speed', 'clear', 'few', 'temp']
    # missing data removed coded -9999
    df.loc[df['dir'] == -9999, 'dir']       = np.mean(df.loc[df['dir'] != -9999, 'dir'])
    df.loc[df['speed'] == -9999, 'speed']   = np.mean(df.loc[df['speed'] != -9999, 'speed'])
    df.loc[df['temp'] == -9999, 'temp']     = np.mean(df.loc[df['temp'] != -9999, 'temp'])
    df.loc[df['clear'] == -9999, 'clear']   = np.mean(df.loc[df['clear'] != -9999, 'clear'])
    df.loc[df['few'] == -9999, 'few']       = np.mean(df.loc[df['few'] != -9999, 'few'])
    # generate variable % cloud cover:
    df['cloud'] = 1 - (df.clear + df.few)/100
    df.loc[df.cloud.isnull(), 'cloud'] = np.mean(df.cloud.dropna())
    # generate day indicator:
    df['daytime'] = 0
    df.loc[(df.hour > 6) & (df.hour < 19), 'daytime'] = 1
    # aggregate daily day-night cycle means
    df = df.groupby(['id','month','day','daytime','lat','lon'])[['dir', 'speed', 'cloud', 'temp']].mean()
    df.reset_index(level=[0,1,2,3,4,5], inplace=True)
    # convert temperature from Fahrenheit to Kelvin:
    df['temp'] = (df.temp - 32)*(5/9) + 273
    # populate weather data matrix:
    weather_data.append(df)

weather_data = pd.concat(weather_data)

weather_data.loc[weather_data.cloud.isnull(), 'cloud'] = np.mean(weather_data.cloud.dropna())

weather_data.to_csv('weather_data.csv', index=False)