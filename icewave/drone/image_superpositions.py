# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:22:25 2025

@author: sebas
"""

import os
import re
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
import scipy
import csv

import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.field.gps as field_gps

import icewave.field.drone as field_drone
import icewave.gps.gps as gps_pack

#%% Get filelist of images 

date = '0210'
drone_ID = 'bernache'
exp_ID = '09-ortho_004'

base = 'U:/Data/'

records = field_drone.get_records('0210',jpg = False)
 # = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'
 
rec = records['drones'][drone_ID]['Drones'][8]
print(rec.keys())

filelist = glob.glob(f'U:/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/*.tiff')
#%% convert .csv flightrecords to .pkl files 
field_drone.convert_flightrecords('0210')

#%% Load flightrecord

file2load = f'{base}{date}/Drones/{drone_ID}/flightrecords/Flightrecord_dict.pkl'
with open(file2load,'rb') as pf:
    records_csv = pickle.load(pf)

# cut flight records using .srt file t_start and t_end 
records_csv = field_drone.cut_flightrecord(rec, records_csv)

#%% Set drone parameters 

param = {}
param['focal'] = 2700
param['h_drone'] = 60 
param['alpha'] = np.pi/2
param['azimuth'] = 263

#%% Convert csv time in time since epoch 

# conversion for csvflightrecord
records_csv['datetime'] = field_drone.get_datetime_from_csvflightrecord(records_csv)
records_csv['t_epoch'] = [d.timestamp() for d in records_csv['datetime']]

# conversion of SRT records
rec['datetime'] = field_drone.get_datetime_from_srtrecord(rec)
rec['t_epoch'] = [d.timestamp() for d in rec['datetime']]

#%% interpolate GPS coordinates and azimuth along time 
I = {}
I['latitude'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.latitude'])
I['longitude'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.longitude'])
I['azimuth'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.yaw [360]'])


#%% Plot several frames 

alpha_0 = param['alpha']
focale = param['focal']


fig, ax = plt.subplots()
N = 1000
step = 100
indices_list = np.arange(0,N,step = step)
extents = []

cmap = ['Reds','Greens','Blues']
for idx in indices_list:
    img_name = filelist[idx]
    img = cv.imread(img_name)
    img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
    
    H = dp.get_height_from_record(rec,indice = idx//100)
    
    Lat_D = I['latitude'](rec['t_epoch'][idx//100])
    Long_D = I['longitude'](rec['t_epoch'][idx//100])
    azimuth = I['azimuth'](rec['t_epoch'][idx//100])
    
    print(idx,H,Lat_D,Long_D,azimuth)
    
    # Lat_D = rec['latitude'][idx//100]
    # Long_D = rec['longitude'][idx//100]
    # print(H,Lat_D,Long_D)
    
    GPS_D = (Lat_D,Long_D,azimuth)
    Lat,Long = dp.georeference_from_param(img,H,alpha_0,focale,GPS_D)
    
    extent = [Long.min(),Long.max(),Lat.min(),Lat.max()]
    extents.append(extent)

    # ax.pcolormesh(Long,Lat,img[:,:,0], cmap = 'gray',rasterized = True)
    ax.pcolormesh(Long,Lat,img[:,:,1],cmap = 'gray',rasterized = True)
    # ax.imshow(np.transpose(img,(1,0,2)),extent = extent,zorder = 2)

extents = np.array(extents)
long_min = np.min(extents[:,0])
long_max = np.max(extents[:,1])
lat_min = np.min(extents[:,2])
lat_max = np.max(extents[:,3])

ax.set_xlim([long_min,long_max])
ax.set_ylim([lat_min,lat_max])

ax.set_aspect(1/np.cos(Lat_D*np.pi/180)) # scaling y/x
# ax.set_aspect('equal')

fig_folder = 'U:/Data/Summary/DAS/'
figname = f'{fig_folder}Images_superposition_DAS_{date}_{drone_ID}_{exp_ID}_N{N}_step{step}_Green'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%%

ny, nx = 200, 300
lon = np.linspace(-0.5, 0.5, nx)
lat = np.linspace(51.0, 51.5, ny)
lon, lat = np.meshgrid(lon, lat)

Z = np.zeros((ny, nx, 3))
Z[..., 0] = np.linspace(0, 1, nx)  # red gradient
Z[..., 1] = np.linspace(1, 0, ny)[:, None]  # green gradient
Z[..., 2] = 0.5  # constant blue

facecolors = Z.reshape(-1, 3)
print(facecolors.shape)

fig, ax = plt.subplots(figsize=(8, 6))

# pcolormesh expects scalar values, so for RGB images, we loop over color channels
# but since it's already RGB, we can use imshow-style plotting using pcolormesh via normalization
# Simplify: flatten and pass as facecolors
facecolors = Z.reshape(-1, 3)
ax.pcolormesh(lon, lat, Z[:,:,0],cmap = 'Reds', shading='auto', alpha = 0.33)
ax.pcolormesh(lon, lat, Z[:,:,1],cmap = 'Greens', shading='auto', alpha = 0.33)
ax.pcolormesh(lon, lat, Z[:,:,2],cmap = 'Blues', shading='auto', alpha = 0.33)


#%%
ax.pcolormesh(Long,Lat,img[:,:,i])