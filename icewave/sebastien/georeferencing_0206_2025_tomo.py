# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 12:56:30 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
from scipy.io import loadmat

import icewave.tools.matlab2python as mat2py
import icewave.drone.drone_projection as dp 

#%%
date = '0206'
drone_ID = 'bernache'
exp_ID = '06-situation_map_geophones_002'

path2data = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/movie_1/'
path2img = f'U:/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/movie_1/'

fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(f'{path2img}*.tiff')
file2read = filelist[0]

#%% Set parameters of drone 

param_dict = {}
param_dict['h_drone'] = 46.4 # drone height
param_dict['alpha_0'] = 37.4*np.pi/180 # rad

param_dict['latitude'] = 48.34689
param_dict['longitude'] = -68.81791

# create UTC_t0
local_time = datetime.strptime('2025-02-06 14:48:05.440' + '000','%Y-%m-%d %H:%M:%S.%f')
tz_quebec = pytz.timezone('America/Montreal')
tz_utc = pytz.utc
local_time = tz_quebec.localize(local_time)
UTC_time = local_time.astimezone(tz_utc)

param_dict['UTC_t0'] = UTC_time

file2save = f'{path2data}Param_drone_{date}_{drone_ID}_{exp_ID}_movie_1.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(param_dict,pf)

#%% 

img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

fig,ax = plt.subplots()
ax.imshow(img)


#%% Set pixel coordinates of geophones 






