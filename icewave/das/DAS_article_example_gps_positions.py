# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 09:11:39 2025

@author: sebas
"""

import numpy as np 
import pickle
import pytz

import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.rw_data as rw
import icewave.geophone.gps_coordinates as geophone_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

#%% Load DAS GPS position 

# change path to GPS DAS positions
file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

file_water_height = file_water_height.replace('.h5','.pkl')
with open(file_water_height,'rb') as pf:
    DAS_water_height = pickle.load(pf)
    
#%% Load DAS results and GPS position

# each active/passive methods have its own dates. 
# in passive method, 'corrected' key corresponds to the correction due to swell orientation

filename = 'U:/Data/Summary/DAS/main_results_active_passive_V2.h5'
results = rw.load_dict_from_h5(filename)

filename = filename.replace('.h5','.pkl')
with open(filename,'rb') as pf:
    results = pickle.load(pf)

#%% Load Stephane and Maddie's GPS measurements

# change paths if needed
# with .h5
filename = 'U:/Data/Summary/DAS/Stephane_gps_measurements_DAS.h5'
gps_results = rw.load_dict_from_h5(filename)

filename = 'U:/Data/Summary/DAS/Maddie_gps_measurements_DAS.h5'
MS_gps = rw.load_dict_from_h5(filename)

# with pickle
filename = 'U:/Data/Summary/DAS/Stephane_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    gps_results = pickle.load(pf)
    
filename = 'U:/Data/Summary/DAS/Maddie_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    MS_gps = pickle.load(pf)
    
#%% Load Geophones GPS positions 

dates = ['0210','0211','0212']

GPS = {'Geophones':{},'thickness':{}}
date = '0210'
year = '2025'

# change paths if needed
path2logfile = f'U:/Data/{date}/Geophones'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
UTC_timearea = pytz.timezone('UTC')

# load geophone GPS coordinates 
new_matrix = geophone_gps.get_GPS_coordinates(path2logfile, geophones_table_path, date, year)
    
#%% Compute GPS_logs for each deployment 

GPS_geoph = np.zeros((len(new_matrix.keys()),len(new_matrix['G01'].keys()),2))

for i,geo_key in enumerate(new_matrix.keys()):
    for j,acqu_key in enumerate(new_matrix[geo_key].keys()):
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
            
            GPS_geoph[i,j,0] = GPS_logs['longitude']
            GPS_geoph[i,j,1] = GPS_logs['latitude']
        else:
            
            GPS_geoph[i,j,0] = GPS_logs['longitude'][0]
            GPS_geoph[i,j,1] = GPS_logs['latitude'][0]
                            
        print(GPS_logs)

