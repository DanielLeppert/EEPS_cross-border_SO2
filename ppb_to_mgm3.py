# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:48:40 2022

@author: danie
"""

# PPM to mgs / cubic meter converter
def ppb_to_mgm3(ppb, mmass):
  
  mgm3 = ppb*mmass*.0409/1000
  
  return mgm3