# -*- coding: utf-8 -*-
"""
Created on Sat May 21 16:25:57 2022

@author: danie
"""

import numpy as np

def calc_sigmas(stab, downwind):

    dist=np.abs(downwind)
    
    # From Pasquill (1961):
    Iy = np.array([-1.104,-1.634,-2.054,-2.555, -2.754,-3.143])
    Jy = np.array([0.9878,1.0350,1.0231,1.0423,1.0106,1.0148])
    Ky = np.array([-0.0076,-0.0096,-0.0076,-0.0087,-0.0064,-0.0070])
    Iz = np.array([4.679,-1.999,-2.341,-3.186,-3.783,-4.490])
    Jz = np.array([-1.7172,0.8752,0.9477,1.1737,1.3010,1.4024])
    Kz = np.array([0.2770,0.0136,-0.0020,-0.0316,-0.0450,-0.0540])

    Iy = Iy[stab]
    Jy = Jy[stab]
    Ky = Ky[stab]
    Iz = Iz[stab]
    Jz = Jz[stab]
    Kz = Kz[stab]
 
    sig_y = np.exp(Iy + Jy*np.log(dist+1e-15) + Ky*(np.log(dist+1e-15)**2))
    sig_z = np.exp(Iz + Jz*np.log(dist+1e-15) + Kz*(np.log(dist+1e-15)**2))

    return (sig_y, sig_z)
