# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 18:39:08 2022

@author: Daniel Leppert
"""
# =============================================================================
# This file generates the data panel for the working paper "Is Geography Destiny? 
# Evidence on the impact of spatially non-targeted emission caps on cross-border 
# externalities from the Clean Air Interstate Rule". The dataset is made up of
# publicly available survey- and monitoring products by the U.S. Environmental 
# Protection Agency, the U.S. Energy Information Administration, and National 
# Oceanic and Atmospheric Administration. Complete citations to these data are 
# listed in the working paper bibliography. 
# 
# Before running this program, ensure that all associated folders and files are
# included as specified in the documentation and that no changes have been made 
# to file names or data. Running the file directly from the master folder 
# [Leppert2020Geog] available to download from [https://github.com/DanielLeppert] 
# ensures this. REPRODUCTION FOLDER NOT YET PUBLIC
# -----------------------------------------------------------------------------
# Section 1: Loads data on plant characteristics, sulfur dioxide emissions and
# sulfur emissions control installations. Converts relevant variables to SI units. 
# Restricts the data to coal-fired electric utilities 1997-2020.
# -----------------------------------------------------------------------------
# Section 2: Loads weather data from the Global Historical Climatology Network 
# which includes typical hourly observations from 460 U.S. weather stations. 
# Computes bi-daily weighted averages for daytime and nightime.
# -----------------------------------------------------------------------------
# Section 3: Links each electric utility to its most proximate weather station 
# using a Haversine function of distance between geographic point coordinates.
# -----------------------------------------------------------------------------
# Section 4: Initializes a GaussMod input matrices init_sulfur and init_weather
# and populates with monthly emissions-, weather- and utility variables.
# See inputs_readme.txt for input data documentation.
# -----------------------------------------------------------------------------
# Section 5: Creates a spatial dataframe and coordinate reference system so that
# dispersion model output can be mapped onto geographic space centered on the plant.
# -----------------------------------------------------------------------------
# Section 6: Runs GaussMod (Leppert 2022) using input matrices and computes the 
# cross-state pollution from each utility in each month using the GaussMod
# output. GaussMod is a Gaussian air pollution dispersion model adapted from 
# Abdel-Rahman (2008) and EPA. GaussMod is imported from 'gaussian_model.py' and relies
# on 'stability.py' to calculate atmospheric stability. The user can change the
# spatial resolution of GaussMod. The default is 1,000 meters. A higher resolution
# will result in a slower runtime.
# -----------------------------------------------------------------------------
# Section 7: Adds annual variables and treatment indicators.
#
# =============================================================================

#%% Section 1 -----------------------------------------------------------------

# import libraries
import pandas as pd
import geopandas as gpd
import numpy as np
import os
import sys
    
try:
    # set working directory
    path = input('ENTER WORKING DIRECTORY:_')
    os.chdir(path)
    print('----------------------------------------\nPreparing to load data...\n----------------------------------------')
    
    # load plant characteristics: EIA id, stack height, stack exit area, gas exit temperature, gas exit velocity
    stack_data = pd.read_csv("plants/stack_data.csv", sep=",")
    stack_data = stack_data.groupby(['id'])[['stack_height','area','exit_temp','exit_vol']].mean().reset_index()
    
    # load emissions and sulfur control data
    emit_files=os.listdir('./sulfur/')
    emissions = []
    for i in range(len(emit_files)):
        df = pd.read_csv('./sulfur/'+emit_files[i], index_col=False)
        # generate a dummy variable indicating sulfur control installed
        df['sulfur_control'] = df['sulfur_control'].fillna(0)
        df.loc[df['sulfur_control'] != 0, 'sulfur_control'] = 1
        # aggregate units by plant
        df = df.groupby(['id','year','month']).agg(
            {
                 # emissions summed over units
                 'sulfur' : sum,   
                 # heat input summed over units
                 'heat_input' : sum,  
                 # operating time summed over units
                 'op_time' : sum,
                 # share of units equipped with desulfurization tech as percentage
                 'sulfur_control' : 'mean'
            }
        )
        df = df.reset_index()
        emissions.append(df)
    
    # generate unbalanced panel
    sulfur_data = pd.concat(emissions).sort_values(['id','year','month'])
    
    # load plant id's and coordinates
    plants_xy = pd.read_csv("plants/2020_1995_facilities.txt", sep=",", index_col=False)
    plants_xy['op_year'] = plants_xy.op_date.str[6:10].astype(float)
    plants_xy = plants_xy[    # only keep plants in operation before 1997
                              (plants_xy.op_year < 1997) &
                              # only keep electric utilities
                              (plants_xy.category == 'Electric Utility') &
                              # only keep plants currently in operation
                              (plants_xy.op_status == 'Operating')]
    
    # only keep Acid Rain Program Phase I-affected units
    plants_xy['include'] = 0
    plants_xy.loc[(plants_xy.year < 2000) & (plants_xy.programs.str.contains('ARP', na=False)), 'include'] = 1
     
    # generate dummy variable: covered by CAIR = 1, otherwise 0
    plants_xy['treated'] = 0
    plants_xy.loc[plants_xy['programs'].str.contains("CAIRSO2", na=False),'treated'] = 1
    
    # generate fuel type dummy 
    plants_xy['fuel'] = 2
    plants_xy.loc[plants_xy['fuel_type'].isin(['Pipeline Natural Gas', 'Natural Gas']), 'fuel'] = 1
    plants_xy.loc[plants_xy['fuel_type'] == 'Coal', 'fuel'] = 0
      
    plants_xy = plants_xy.groupby(['id','state','county_fips','lat','lon'])[['treated','fuel','include']].mean()
    plants_xy.reset_index(level=[0,1,2,3,4], inplace=True)
    plants_xy = plants_xy[plants_xy.include > 0.]
    plants_xy.loc[plants_xy.treated > 0, 'treated'] = 1
    
    # merge data
    data = sulfur_data.merge(plants_xy, how='inner', on='id')
    data = data.merge(stack_data, how='inner', on='id')
    data = data[data.year >= 1997]
    data.loc[data.op_time == 0, 'sulfur'] = 0   
                  
    # variable unit conversion:
    
    # SO2 emissions from tonnes per year to grams per second:
    data['sulfur_rate'] = 0
    data.loc[data.op_time > 0, 'sulfur_rate'] = (data.sulfur*1e+6)/(data.op_time*60*60)
    # gas exit temperature from degrees Fahrenheit to degrees Kelvin:
    data['exit_temp'] = (data.exit_temp - 32)*(5/9) + 273
    # stack height from feet to meters:
    data['stack_height'] = data.stack_height*.3048
    # stack exit radius (m) from stack exit area (sqare feet):
    data['radius'] = np.sqrt(data.area*0.09/np.pi)
    # gas exit velocity feet/second to meters/second:
    data['exit_vol'] = data.exit_vol*.3048
       
    # keep only complete years:
    for i in np.unique(data.id):
        for j in np.unique(data.year):
            if len(data[(data.id == i) & (data.year == j)]) < 12:
                data = data.drop(data[(data.id == i) & (data.year == j)].index)   
                   
    plants_xy = data[['lat','lon','id','state']].drop_duplicates()
    
    # Lambert conal projection
    utm = "+proj=utm +zone=17 +datum=NAD83 +units=m"
    # transform to spatial data
    plant_m = gpd.GeoDataFrame(plants_xy,
                               geometry = gpd.points_from_xy(plants_xy.lon, plants_xy.lat),
                               crs = 4269)
    
    plant_m = plant_m.to_crs(utm)
    
    nplants = len(plants_xy)
    
    # extract plant coordinates in meters
    lats = np.floor(np.array([i for i in plant_m.geometry.y]))
    lons = np.floor(np.array([i for i in plant_m.geometry.x]))
      
    # section 1: clean working directory
    del sulfur_data, df, emissions, stack_data
    print('----------------------------------------\nData loaded successfully\n----------------------------------------')


except:
    print('DIRECTORY ENTERED IS NOT VALID.')
    sys.exit(1)

#%% Section 2 -----------------------------------------------------------------
  
# retrieving data from file.

years = [i for i in range(1997,2021)]

wea_dat = []
for i in range(24):
    dat = pd.read_csv("weather_data/weather_data.csv")
    dat['year'] = years[i]
    wea_dat.append(dat)
    
wea_dat = pd.concat(wea_dat)

tmp_dat = pd.read_csv("weather_data/temp_data.txt", sep = "\t").drop(['lat','lon'], axis = 1)

wea_dat = wea_dat.merge(tmp_dat, on = ['id','year','month','day','daytime'], how = 'left')

# replace missing temperature observations with the normals
wea_dat.loc[wea_dat.temp.isnull(), 'temp'] = wea_dat.loc[wea_dat.temp.isnull(), 'temp0']
stations_xy = wea_dat[['id','lon','lat']].drop_duplicates().dropna()

del tmp_dat, dat
  
# section 2: clean working directory

#%% Section 3 -----------------------------------------------------------------

# import Haversine equation for distances between points
from distance import dist_haversine

closest = []
# for each plant 'i'...
for i in range(nplants):
    # compute distance from 'i' to each station 'j'
    dists = np.zeros((461))
    for j in range(461):
        # run haversine program:
        dists[j] = dist_haversine(plants_xy.lon.iloc[i], plants_xy.lat.iloc[i], stations_xy.lon.iloc[j], stations_xy.lat.iloc[j])
    # for each plant 'i' select the closest station
    df = stations_xy.id.iloc[np.argmin(dists)]
    closest.append(df)

del df
print('----------------------------------------\nPlant-station distance matrices computed successfully\n----------------------------------------')

#%% Section 4 -----------------------------------------------------------------

# This section of code generates the input array for GaussMod. GaussMod takes 9 input variables, plus an elevation constant.
# See gauss_input.readme.txt for further documentation on the input file.

# number of days in each month
days = [31,28,31,30,31,30,31,31,30,31,30,31]

# functions to initialize a matrices of inputs for each month

def init_weather(wea_dat, i, j, k):
    
    # initialize matrix of weather inputs:
    wmat = np.zeros((days[k]*2,5)) 
      
    wmat[:,:] = wea_dat.loc[(wea_dat.year == np.unique(wea_dat.year)[j]) & (wea_dat.month == k+1) & (wea_dat.id == closest[i]),\
                            ['dir','speed','cloud','temp','daytime']].to_numpy()           
        
    return wmat

    
def init_sulfur(data, i, j, k):
    
    # initialize matrix of plant inputs:
    pmat = np.zeros(5)
        
    pmat[:] = data.loc[(data.id == np.unique(data.id)[i]) & \
                         (data.year == np.unique(data[data.id == np.unique(data.id)[i]].year)[j]) & \
                         (data.month == k+1), \
                    ['sulfur_rate','stack_height','radius','exit_temp','exit_vol']].to_numpy()
        
    return pmat

print('----------------------------------------\nInput matrices initialized successfully\n----------------------------------------')
      
#%% Section 5 -----------------------------------------------------------------

# load US shapefile
us = gpd.read_file('./shapefiles/cb_2018_us_state_20m.shp')
# filter continental US:
us_cont = us[-us['STUSPS'].isin(['AK','HI','PR'])]
us_cont = us_cont.to_crs(epsg=4269)
# reproject to meters from degrees:
us_cont_m = us_cont.to_crs(utm) 
# generate externality variables showing concentration from stack outside state borders
data['total_sulfur']=0
data['ext_sulfur_sum']=0
data['ext_sulfur_mean']=0
data['ext_sulfur_max']=0

#%% Section 6 ----------------------- MAIN LOOP ------------------------------

# run the GaussMod (Leppert 2022) function: For each power plant, calculate SO2 dispersion in a 50,000 by 50,000 m grid   

from annual_disp import annual_dispersion
import time

start = int(input("ENTER START POSITION: "))

print('\nGaussMod is calculating hourly SO2 dispersion...\nThis will take a while.\nEnsure that your machine has power plugged IN and sleep mode turned OFF')
print('\nWHILE YOU WAIT: Geophysical modelling relies on complete spatial data. When data is missing, spatial interpolation can be used.\nFor more on how to address geographic basis errors in spatial datasets see: Leppert, Dalhaus and Lagerkvist (2021)\n"Accounting for geographic basis risk in heat index insurance: how spatial interpolation can reduce the cost of risk"\nWeather, Climate and Society, 13(2).\n ')

begin = time.perf_counter()
# for each plant i
for i in range(start, nplants):
    # and each year j 
    for j in range(len(np.unique(data[data.id == np.unique(data.id)[i]].year))):
        # and each month k
        for k in range(12):
                     
           pmat = init_sulfur(data, i, j, k)
           wmat = init_weather(wea_dat, i, j, k)
            
           # check that matrices contain all needed data: 
           try:
                # calculate mean annual dispersion
                data = annual_dispersion(lats[i], lons[i], i, j, k, pmat, wmat, plant_m, us_cont_m, data)
            
           # or terminate code with error message:
           except:
                print("DATA ERROR IN ID:", np.unique(data.id)[i], "YEAR:", np.unique(data[data.id == np.unique(data.id)[i]].year)[j], 'MONTH:', k+1)
                sys.exit(1)
                
                
    # print update message to screen
    end = time.perf_counter()
    duration = end-begin
    print(round((i+1)/nplants*100, 1), '% computed. Expected time to completion:', round((duration/(i+1))*(nplants-(i+1))/60/60, 1), 'hours')


# aggregate monthly data to annual
data = data.groupby(['id','year','state','county_fips','lat','lon','treated']).agg(
    {
         # emissions summed over months
         'sulfur' : sum,   
         # heat input summed over months
         'heat_input' : sum,  
         # operating time summed over units
         'op_time' : sum,
         # average annual aggregation of monthly values:
         'sulfur_control' : 'mean',
         'sulfur_rate' : 'mean',
         'fuel' : 'mean',
         'total_sulfur' : 'mean',
         'ext_sulfur_sum' : 'mean',
         'ext_sulfur_mean' : 'mean',
         'ext_sulfur_max' : 'mean'   
    }
)

data = data.reset_index()

print('----------------------------------------\nGaussMod ran successfully\n----------------------------------------')


#%% Section 7 -----------------------------------------------------------------

# load carbon emissions
carbon_df = pd.read_csv('./carbon/carbon_annual.csv', index_col=False)

# load acid rain program permit holdings data
arp_files = os.listdir('./holdings')
arp_list = []
for file in range(len(arp_files)):
    df = pd.read_csv('./holdings/' + arp_files[file])
    df = df[['ORISPL_CODE', 'OP_YEAR', 'ALLOCATED', 'TOTAL_HELD']]
    df.columns = ['id', 'year','allocated', 'arp_permits']
    df = df[df.id.isin(data.id)]
    arp_list.append(df)
    
# Load CAIR and CSAPR data
cair_files = os.listdir('./cair')
cs_g1_files = os.listdir('./csapr1')
cs_g2_files = os.listdir('./csapr2')

cair_list = []
cs_g1_list = []
cs_g2_list = []

# number of permits held per utility and year within the CAIR program
for file in range(len(cair_files)):
    df = pd.read_csv('./cair/' + cair_files[file])
    df = df[['ORISPL_CODE', 'OP_YEAR','TOTAL_HELD']]
    df.columns = ['id', 'year', 'cair_permits']
    df = df[df.id.isin(data.id)]
    cair_list.append(df)
# number of permits held per utility and year within CSAPR program, group 1    
for file in range(len(cs_g1_files)):
    df = pd.read_csv('./csapr1/' + cs_g1_files[file])
    df = df[['ORISPL_CODE', 'OP_YEAR','TOTAL_HELD']]
    df.columns = ['id', 'year', 'cs1_permits']
    df = df[df.id.isin(data.id)]
    cs_g1_list.append(df)
# number of permits held per utility and year within CSAPR program, group 2    
for file in range(len(cs_g2_files)):
    df = pd.read_csv('./csapr2/' + cs_g2_files[file])
    df = df[['ORISPL_CODE', 'OP_YEAR','TOTAL_HELD']]
    df.columns = ['id', 'year', 'cs2_permits']
    df = df[df.id.isin(data.id)]
    cs_g2_list.append(df)

cair_aggregate = pd.concat(cair_list).groupby(['id', 'year'])[['cair_permits']].sum()
cs_g1_aggregate = pd.concat(cs_g1_list).groupby(['id', 'year'])[['cs1_permits']].sum()
cs_g2_aggregate = pd.concat(cs_g2_list).groupby(['id', 'year'])[['cs2_permits']].sum()
arp_aggregate = pd.concat(arp_list).groupby(['id', 'year'])[['allocated','arp_permits']].sum()

arp_aggregate.reset_index(level=[0,1], inplace=True)
cair_aggregate.reset_index(level=[0,1], inplace=True)
cs_g1_aggregate.reset_index(level=[0,1], inplace=True)
cs_g2_aggregate.reset_index(level=[0,1], inplace=True)

del arp_list, cair_list, cs_g1_list, cs_g2_list
#%%
    
# load monthly plant-level generation data
gen_files = os.listdir('./netgen')
gen_list = []
for i in range(len(gen_files)):
    df = pd.read_csv('./netgen/' + gen_files[i])
    
    if i <= 3:
        df = df[['ID','YEAR','GEN01','GEN02','GEN03','GEN04','GEN05','GEN06','GEN07','GEN08','GEN09','GEN10','GEN11','GEN12']]        
    else:
        df = df.iloc[:,[0,13] + [j for j in range(1,13)]]
    
    df.replace('.', np.nan, inplace=True)
    df.columns = ['id', 'year', 'gen01','gen02','gen03','gen04','gen05','gen06','gen07','gen08','gen09','gen10','gen11','gen12']
    df = df.apply(pd.to_numeric)
    df = df.groupby(['id', 'year'])[['gen01','gen02','gen03','gen04','gen05','gen06','gen07','gen08','gen09','gen10','gen11','gen12']].sum()
    df.reset_index(level=[0,1], inplace=True)
    
    cols = ['gen01','gen02','gen03','gen04','gen05','gen06','gen07','gen08','gen09','gen10','gen11','gen12']
    
    df = pd.melt(df, id_vars=['year','id'], value_vars=cols, value_name='generation')

    months = [m for m in range(12)]
    
    for m in range(12):
        df.loc[df.variable == cols[m], 'variable'] = months[m]+1
        
    gen_list.append(df)
    
gen_df = pd.concat(gen_list)
gen_df.columns = ['year','id','month','generation']
gen_df = gen_df.groupby(['id','year'])[['generation']].sum()
gen_df = gen_df.reset_index()

del gen_list, df

#%%

# load annual retail electricity prices:
retail_price = pd.read_csv('./prices/avgprice_annual.csv', header=0)
retail_price = retail_price[retail_price['Industry Sector Category'] == 'Total Electric Industry']
retail_price = retail_price[retail_price['Year'] >= 1997]

retail_price = retail_price[['Year', 'State', 'Total']]
retail_price.columns = ['year', 'state', 'price']

# load SO2 auction spot price data:
spot_price = pd.read_csv('./prices/SO2_spot_prices.csv', header = 0)

# %%

from ppb_to_mgm3 import ppb_to_mgm3

# add variables to final panel dataset
data_final = data.merge(spot_price, on = 'year', how='left')
data_final = data_final.merge(retail_price, on = ['year', 'state'], how='left')
data_final = data_final.merge(carbon_df, on = ['id','year'], how='left')
data_final = data_final.merge(gen_df, on = ['id', 'year'], how='left')
data_final = data_final.merge(arp_aggregate, on = ['id', 'year'], how='left')
data_final = data_final.merge(cair_aggregate, on = ['id', 'year'], how='left')
data_final = data_final.merge(cs_g1_aggregate, on = ['id', 'year'], how='left')
data_final = data_final.merge(cs_g2_aggregate, on = ['id', 'year'], how='left')

data_final['fuel_type'] = 'Majority Oil'
data_final.loc[data_final.fuel < 1, 'fuel_type'] = 'Majority Coal'
data_final.loc[data_final.fuel == 1, 'fuel_type'] = 'Majority Gas'

# ARP permit program ends
data_final.loc[(data_final.year >= 2010) & (data_final.treated == 1), ['allocated', 'arp_permits']] = 0
       
# EPA permit holdings data is complete therefore utilities missing records held no permits in a given year.
data_final[['allocated', 'arp_permits', 'cair_permits', 'cs1_permits', 'cs2_permits']] = data_final[['allocated', 'arp_permits', 'cair_permits', 'cs1_permits', 'cs2_permits']].fillna(0)

data_final['permits'] = 0

data_final['permits'] = data_final[['arp_permits','cair_permits','cs1_permits','cs2_permits']].sum(axis=1)       

# generate indicator variable if CAIR had been announced post-2005:
data_final['CAIR'] = 0
data_final.loc[(data_final.year > 2005), 'CAIR'] = 1
# generate indicator variable if CAIR had been implemented post-2009:
data_final['CSAPR'] = 0
data_final.loc[data_final.year > 2009, 'CSAPR'] = 1
# generate indicator variable for the North Carolina v. EPA ruling of 2008:
data_final['NCvEPA'] = 0
data_final.loc[data_final.year > 2008, 'NCvEPA'] = 1
# generate NAAQS violation indicator:
data_final['naaqs_violate'] = 0
data_final.loc[data_final.ext_sulfur_max*1000 > ppb_to_mgm3(.75, 64.07), 'naaqs_violate'] = 1              
# write final results to csv:
data_final.to_csv('data_final.csv', header=True, sep=',', index=False)
print('----------------------------------------\nProgram terminated successfully\n----------------------------------------')
