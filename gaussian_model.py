# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 15:31:26 2021

@author: daniel
"""
# =============================================================================
# VARIABLE DICTIONARY:
# 
# dxy      = spatial resolution (m^2), default 1,000 square meters
# x        = x-coordinates (m)
# y        = y-coordinates (m)
# z        = distance above ground (m)
# e        = pollutant gas emission rate (g/s)
# H_s      = height of smoke stack (m)
# r        = radius of smoke stack exit (m)
# T_gas    = gas exit temperature (degrees K)
# u_gas    = gas exit velocity (m/s)
# v        = wind direction (degrees, v=0 implies wind from south to north)
# u        = 2-minute peak hourly wind speed (m/s)
# cloud    = cloud cover (%)
# T_air    = ambient average temperature (degrees K)
# C        = pollutant concentration (g/m^3) in location (x,y,z) 
# F        = plume bouyancy flux (m^4/s^3) determined by the exit force of the gas and the force exerted by the air upon it as it rises
# dH       = plume rise (m)
# 
# =============================================================================

def gauss_mod(dxy, x, y, z, e, H_s, r, T_gas, u_gas, v, u, cloud, T_air, daytime):
    
    import numpy as np
   
    # this module computes atmospheric stability from Pasquill (1961)
    from stability import atmospheric_stability
    from calc_sigmas import calc_sigmas
    
    # components of wind speed 'u' in x and y directions
    wx=u*np.sin((v)*np.pi/180.);
    wy=u*np.cos((v)*np.pi/180.);    
    # need angle between point (x, y) and the wind direction, so use scalar product:
    dot_product=wx*x+wy*y;
    # product of magnitude of vectors: 
    magnitudes=u*np.sqrt(x**2 + y**2); 
    # angle between wind direction vector and point (x,y)
    subtended=np.arccos(np.clip(dot_product/(magnitudes+1e-15), -1, 1));
    # distance to point (x,y) from stack
    hypotenuse=np.sqrt(x**2+y**2);   
    # distance along the wind direction to perpendicular line that intesects (x,y)
    downwind=np.cos(subtended)*hypotenuse;

    # now calculate distance cross wind.
    crosswind=np.sin(subtended)*hypotenuse;
    ind=np.where(downwind > 0);
    C=np.zeros((len(x), len(y)));

    # determine stability class:
    stability = atmospheric_stability(u, cloud, daytime)

    # calculate plume rise:

    # gravity constant (m/s^2):
    gravity = 9.81

    # bouyancy factor:
    F = gravity * u_gas * (r**2) * (T_gas - T_air)/T_air

    # computes the peak height above ground for the plume centerline
    if F < 55:
        x_max = 49*F**(5/8)
    else:
        x_max = 119*F**(2/5)
                
    # computes the plume rise above stack height at every pixel downwind
    dH = np.zeros((len(x), len(y)))
    dH[:,:] = 1.6*(F**(1/3))*(x_max**(2/3))/u
        
    # total elevation of plume centerline above ground level
    H = H_s + dH
    
    # calculate sigmas
    (sig_y,sig_z)=calc_sigmas(stability, downwind);
    
    # horizontal dispersion parameter
    f = np.exp(-crosswind[ind]**2./(2.*sig_y[ind]**2.))
    # vertical dispersion parameter
    g = np.exp(-(z-H[ind])**2./(2.*sig_z[ind]**2.))
    g_inv = np.exp(-(z+H[ind])**2./(2.*sig_z[ind]**2.))
    
    C[ind]=(e/u)*f/(np.sqrt(2.*np.pi)*sig_y[ind])*(g + g_inv)/(np.sqrt(2*np.pi)*sig_z[ind])
    
    return C
