# -*- coding: utf-8 -*-
"""
Created on Tue May 17 19:41:29 2022

@author: danie
"""

def annual_dispersion(lat, lon, i, j, k, plant_input, weather_input, plant_m, us_cont_m, data):
    
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    from gaussian_model import gauss_mod
   
    # number of day/night cycles in month k
    cycles = len(weather_input)
    
    # set height above ground (m) where pollution should be estimated
    z = 1.5
    
    # set spatial resolution for GaussMod:
    dxy = 1000.
    x=np.mgrid[-25000.:25000.+dxy:dxy]
    y=np.mgrid[-25000.:25000.+dxy:dxy]
    
    # create grid
    [x,y]=np.meshgrid(x,y)
      
    conc = np.zeros((51,51))
    # for each hour:
    for h in range(cycles):
       # initialize a spatial grid
       c = np.ones((51, 51))            
       # and populate grid with GaussMod output:           
       
       c = gauss_mod(dxy,    # spatial resolution
                       x,    # x coordinates
                       y,    # y coordinates
                       z,    # z coordinate
      plant_input[0],        # emissions rate
      plant_input[1],        # stack height
      plant_input[2],        # stack radius at top
      plant_input[3],        # gas exit temperature
      plant_input[4],        # gas exit velocity
      weather_input[h,0],  # wind direction
      weather_input[h,1],  # wind speed
      weather_input[h,2],  # cloud cover
      weather_input[h,3],  # air temperature                     
      weather_input[h,4])  # daytime indicator
       
       # update concentration matrix
       conc = conc + c
    
    # monthly average concentration at each pixel
    conc = conc/cycles
    
    x = x + lon
    y = y + lat
        
    # transform concentration matrix to spatial pixels
    pixels = pd.DataFrame({'lat': np.repeat(y[:,0], 51), 'lon': x.flatten(order='C'), 'conc': conc.flatten()})    
    pixels = gpd.GeoDataFrame(pixels, geometry=gpd.points_from_xy(pixels.lon, pixels.lat), crs=us_cont_m.crs)   
    pixels = gpd.tools.sjoin(pixels, us_cont_m, op='within', how='left')
            
    # total SO2 dispersion
    total_sulfur = np.sum(np.array(pixels.conc))
    # cross-state concentrations 
    external = np.array(pixels[(pixels.STUSPS != plant_m.state.iloc[i]) & (pixels.STUSPS.notnull())].conc)
    
    data.loc[(data.id==np.unique(data.id)[i]) & \
             (data.year==np.unique(data.year)[j]) & \
             (data.month==np.unique(data.month)[k]), 'total_sulfur'] = total_sulfur
       
    if len(external) != 0:
        # cross-state SO2 sum
        external_sum = np.sum(external)
        # average externality
        external_mean = np.mean(external)
        # highest cross-state concentration
        external_max = np.max(external)
        # populate data
        data.loc[(data.id==np.unique(data.id)[i]) & \
                 (data.year==np.unique(data.year)[j]) & \
                 (data.month==np.unique(data.month)[k]), \
                 ['ext_sulfur_sum','ext_sulfur_mean','ext_sulfur_max']] =  [external_sum, external_mean, external_max]
    
    return data
    
   
