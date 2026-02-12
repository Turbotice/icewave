# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 14:23:59 2026

@author: sebas

Script gathers results obtained from geophone lines 

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pickle 
import scipy

import icewave.geophone.gps_coordinates as geophone_gps
import icewave.geophone.package_geophone as geopack
import icewave.sebastien.set_graphs as set_graphs

#%% Search all results files

# path to geophones table 
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'

year = '2026'
date = '0129' #date format, 'mmdd'
base = 'F:/Rimouski_2026/'
# path2data = os.path.join(base,date,'Geophones/')

path2data = f'{base}{date}/Geophones/'
filelist = glob.glob(f'{path2data}*_results_inversion.pkl')

#%% Load file 

file2load = filelist[0]
with open(file2load,'rb') as pf:
    results = pickle.load(pf)
    print(results.keys())

#%% get acquisition number 
print(file2load)  
acq_numb = file2load.split('acq')[1][:4]

#%% Get gps coordinates of all geophones, for each acquisition

gps_dict = geophone_gps.get_GPS_coordinates(path2data, geophones_table_path, date, year)
GPS_geoph = geophone_gps.build_GPS_array(gps_dict)

#%% Get water height 

#%% Get UTC time of acquisition 
