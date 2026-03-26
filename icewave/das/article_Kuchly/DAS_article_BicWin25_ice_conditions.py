# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 14:10:12 2025

@author: sebas

This script gathers aerial pictures to define three different ice conditions encountered during 
BicWin25 field measurement campain. 
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
import icewave.field.drone as field_drone
import icewave.field.gps as field_gps
import icewave.gps.gps as gps_pack

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))


#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/Measurement_map/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
    
#%% path to drone images and records 
date = '0210'
drone_ID = 'bernache'
exp_ID = '09-ortho_004'

base = 'U:/Data/'

records = field_drone.get_records('0210',jpg = False)
 # = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'
 
rec = records['drones'][drone_ID]['Drones'][8]
print(rec.keys())

filelist = glob.glob(f'U:/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/*.tiff')

#%% Load flightrecord

file2load = f'{base}{date}/Drones/{drone_ID}/flightrecords/Flightrecord_dict.pkl'
with open(file2load,'rb') as pf:
    records_csv = pickle.load(pf)

# cut flight records using .srt file t_start and t_end 
records_csv = field_drone.cut_flightrecord(rec, records_csv)

#%% Set drone parameters 

param = {}
param['focal'] = 2700
param['alpha'] = np.pi/2

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

#%% Ice condition : region 1

alpha_0 = param['alpha']
focale = param['focal']

fig, ax = plt.subplots(figsize = (12,6))

N = 3500 #3500
step = 100 #100
indices_list = np.arange(0,N,step = step)
extents = []
for idx in indices_list:
    img_name = filelist[idx]
    img = cv.imread(img_name)
    img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
    gray_img = cv.cvtColor(img, cv.COLOR_RGB2GRAY)

    mean_gray = np.mean(gray_img)
    norm_img = gray_img/mean_gray
    
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

    # ax.pcolormesh(Long,Lat,img[:,:,0], cmap = 'gray')
    # ax.imshow(np.transpose(img,(1,0,2)),extent = extent)
    ax.pcolormesh(Long,Lat,norm_img,cmap = 'gray',vmin = 0,vmax = 1.2,rasterized = True)


# plot DAS positions
# create a norm and a line collection 
points = np.array([DAS_water_height['long'],DAS_water_height['lat']]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
distance = DAS_water_height['s']
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 2)
lc.set_array(distance) # values used for colormap
ax.add_collection(lc)
# ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k')


# superpose deployments points 
ms = 70

# create norm for ice thickness
norm_acq = colors.Normalize(vmin = 1,vmax = 4)
full_cmap = mpl.colormaps['Oranges'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,1,256)))

color_date = {'0210':new_blues(0.2),'0211':new_blues(0.6), '0212': new_blues(0.9)}

# plot geophone lines and tomo
for geo_key in new_matrix.keys():
    for acqu_key in new_matrix[geo_key].keys():
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
        
        current_color = cmap(norm_acq(acqu_key))
        ax.scatter(GPS_logs['longitude'],GPS_logs['latitude'],marker = '^',color = current_color,
                   edgecolors = 'k',zorder = 2)

# plot ice thickness
# for key_date in records.keys():
#     for key in records[key_date].keys():
#         if 'H' in key:
#             h_value = int(key[2:])/100
#             current_color = cmap(norm_thick(h_value))
#             ax.scatter(records[key_date][key]['longitude'],records[key_date][key]['latitude'], marker = '^',
#                        color = current_color)

for key_date in gps_results.keys():
    current_color = color_date[key_date]
    ax.scatter(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],
               marker = 'p',color = current_color,edgecolors = 'k',zorder = 3)


extents = [-68.821900,-68.813627,48.346588,48.348010]
ax.set_xlim([extents[0], extents[1]])
ax.set_ylim([extents[2], extents[3]])

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

Lat0 = DAS_water_height['lat'][0]
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x










