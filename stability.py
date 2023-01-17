# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 22:22:03 2021

@author: danie
"""

# this function determines the atmospheric stability class from Pasquill (1961)

def atmospheric_stability(u, cloud, daytime):
    
    import pandas as pd
    import numpy as np
    
    
    stability = 0
      
    # classify cloudiness
    if cloud < .25:
        CLOUD = 'LOW'
    
    elif cloud >= .25 and cloud <= .75:
        CLOUD = 'MID'
    
    else:
        CLOUD = 'HIGH'
        
    
    # classify wind speed
    if u < 2:
         WIND = 'LOW'
     
    elif u >= 2 and u < 6:
         WIND = 'MID'
     
    else:
         WIND = 'HIGH'
         
    # is it daytime?
    day = False
    
    if daytime == 1:
        day = True
    
    
    for i in range(6):
        
        if CLOUD == 'LOW' and WIND == 'LOW' and day:
            # extremely unstable
            stability = 0
                                                   
        elif (CLOUD == 'LOW' and WIND == 'MID' and day) or (CLOUD == 'MID' and WIND == 'LOW' and day):
            # moderately unstable
            stability = 1
        
        elif (CLOUD == 'MID' and WIND == 'MID' and day) or (CLOUD == 'LOW' and WIND == 'HIGH' and day):
            # slightly unstable
            stability = 2
       
        elif (CLOUD == 'MID' and WIND == 'HIGH' and day) or (CLOUD == 'HIGH' and WIND == 'LOW' and day):
            # neutral
            stability = 3
            
            # slightly stable
        elif (CLOUD == 'HIGH' and WIND == 'HIGH' and not day): 
            stability = 4
        
        elif (CLOUD != 'MID' and WIND != 'HIGH' and not day):
            # stable
            stability = 5
            
    
    return stability
            