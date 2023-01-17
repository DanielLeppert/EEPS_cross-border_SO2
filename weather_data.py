# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:18:31 2022

@author: danie
"""

# This file 

import pandas as pd
import numpy as np
import os
from time import sleep

os.chdir("C:/Users/danie/Documents/")

# load inventory of GHCN-Monthly weather stations:
inv = pd.read_csv("./ghcnd_inventory.txt", sep = "\s+", header = None, error_bad_lines=False)
inv = inv.rename(columns = {0:'id', 1:'lat', 2:'lon', 3:'code', 4:'begin', 5:'end'})

# keep only stations collecting wind data between 1995 and 2020:
inv_wnd = inv[(inv.begin <= 1995) & (inv.end >= 2020) & (inv.code == 'WDF2')]

stations_wnd = []
for station in inv_wnd.id:
    filedir = './ghcnm/'+station+'.csv'
    stations_wnd.append(filedir)
    
wind_list = []
print("---------------------------------")
print("                                 ")
print('Preparing to load weather data...')
sleep(5)
for i in range(len(stations_wnd)):
    # load daily US wind data from NOAA:
    try:
        df = pd.read_csv(stations_wnd[i], header=0)
        df = df[['STATION', 'LATITUDE', 'LONGITUDE', 'DATE','WDF2', 'WSF2','TAVG']]
        df.columns = ['id', 'lat', 'lon', 'date','dir','speed', 'tavg']
        df['year'] = df.date.str[0:4].astype(int)
        df['month'] = df.date.str[5:7].astype(int)
        df = df[(df.year >= 1995) & (df.year <= 2020)]
        # make sure station contains full time series
        if len(df) == 312:
            wind_list.append(df)
        else:
            inv_wnd = inv_wnd.drop(inv_wnd[inv_wnd.id.isin(df.id)].index)
            continue
    # print warning message if file missing in folder:
    except:
        print('missing file no ', i)
        pass
del df

#%%

usw = inv['id'].str.contains('USW')
inv_hrly = inv[usw]

stations = []
for station in np.unique(inv_hrly.id):
    filedir = './us-hourly/'+station+'.csv'
    stations.append(filedir)
    
wind_hourly = []
print("---------------------------------")
print("                                 ")
print('Preparing to load weather data...')
sleep(5)
for i in range(len(stations)):
    # load daily US wind data from NOAA:
    try:
        df = pd.read_csv(stations[i], header=0)
        df = df[['STATION', 'LATITUDE', 'LONGITUDE', 'm', 'day', 'hour', 'HLY-WIND-VCTDIR', 'HLY-WIND-AVGSPD', 'HLY-CLOD-PCTCLR', 'HLY-CLOD-PCTFEW', 'HLY-TEMP-NORMAL']]
        df.columns = ['id', 'lat', 'lon', 'month', 'day', 'hour', 'dir', 'speed', 'clear', 'few', 'temp']
        df.loc[df['dir'] == -9999, 'dir']       = np.mean(df.loc[df['dir'] != -9999, 'dir'])
        df.loc[df['speed'] == -9999, 'speed']   = np.mean(df.loc[df['speed'] != -9999, 'speed'])
        df.loc[df['temp'] == -9999, 'temp']     = np.mean(df.loc[df['temp'] != -9999, 'temp'])
        df.loc[df['clear'] == -9999, 'clear']   = np.mean(df.loc[df['clear'] != -9999, 'clear'])
        df.loc[df['few'] == -9999, 'few']       = np.mean(df.loc[df['few'] != -9999, 'few'])
        wind_hourly.append(df)
    # print warning message if station missing in file:
    except:
        print('missing file no ', i)
        pass

del df 


