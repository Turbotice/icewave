# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 17:52:12 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
import csv

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Set fig_folder path 

fig_folder = 'F:/Rimouski_2025/Data/Summary/DAS/Situation_tomo_0210/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load DAS GPS position 

base = 'F:/Rimouski_2025/Data/'

file_water_height = f'{base}0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

#%% Load Stephane and Maddie's GPS measurements

# with pickle
filename = 'F:/Rimouski_2025/Data/Summary/DAS/Stephane_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    gps_results = pickle.load(pf)
    
filename = 'F:/Rimouski_2025/Data/Summary/DAS/Maddie_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    MS_gps = pickle.load(pf)
    
#%% Load geophones gps positions
file2load = 'F:/Rimouski_2025/Data/Summary/DAS/Geophones_positions_0210_0212.pkl'
with open(file2load,'rb') as pf:
    GPS_geoph = pickle.load(pf)
    
#%% Load tomography image

base = 'F:/Rimouski_2025/PIV_images/'
date = '0210'
drone_ID = 'bernache'
exp_ID = '11-situation_map_geophone_003'
movie_ID = 'movie_1'

path2tiff = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/{movie_ID}/'
filelist = glob.glob(f'{path2tiff}*.tiff')
file2read = filelist[0]
print(file2read)

img_tomo = cv.imread(file2read)
img_tomo = cv.cvtColor(img_tomo,cv.COLOR_BGR2RGB)

#%% Load situation picture 

date = '0205'
drone_ID = 'mesange'
exp_ID = '18-doc_010'

path2img = f'F:/Rimouski_2025/Data/{date}/Drones/{drone_ID}/{exp_ID}/'
filelist = glob.glob(f'{path2img}*.jpg')
file2read = filelist[0]
print(file2read)

# fig_folder = f'U:/Data/{date}/Drones/{drone_ID}/Figures/{exp_ID}/'
# if not os.path.isdir(fig_folder):
#     os.mkdir(fig_folder)
    
img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

#%% Set drone parameters 

# Set drone parameters, tomo
param_dict = {'tomo':{}}
param_dict['tomo']['H'] = 91.8 #110.8
param_dict['tomo']['alpha_0'] = 90
param_dict['tomo']['focale'] = 2700 

param_dict['tomo']['latitude'] =  48.34755# 48.3475505902513 
param_dict['tomo']['longitude'] =  -68.81794# -68.8179351905923 
param_dict['tomo']['azimuth'] = 351.5

# set drone parameters situation
param_dict['situ'] = {}
param_dict['situ']['H'] = 110.8 #110.8
param_dict['situ']['alpha_0'] = 22.4
param_dict['situ']['focale'] = 4100 # around 4100 (+- 100) for format 5280 x 2970 # initial guess 4100

# -68.8131039300668
# -68.81312780555555 #- from image data

param_dict['situ']['latitude'] =  48.347697722222215 #- from image data #48.3476967601294 - Column BN True 
param_dict['situ']['longitude'] =  -68.81312780555555 #- from image data #-68.8131401977678 - Column BN True 
param_dict['situ']['azimuth'] = 263

#%% Compute GPS position of tomography image 




