# -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:46:25 2025

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
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather


#%% Load picture

main_path = 'U/Data/'
date = '0211'
drone_ID = 'mesange'
exp_ID = '08-stereo_002'

path2img = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/'
filename = 'DJI_20250211154246_0365_D_exemple.tiff'
file2read = f'{path2img}{filename}'

fig_folder = f'U:/Data/{date}/Drones/{drone_ID}/Figures/{exp_ID}/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

img_name = filename.split('.')[0]

#%% Plot raw image

fig, ax = plt.subplots()
ax.imshow(img)
ax.set_xlabel(r'$x_p $',labelpad = 5)
ax.set_ylabel(r'$y_p $',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 


#%% Set objects position (in pixels)

POS = np.zeros((3,2))
POS[0,:] = [1231.4,1358.5]
POS[1,:] = [1951,352]
POS[2,:] = [3503,780]

# superpose objects on figure 
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.imshow(img)
ax.plot(POS[:,0],POS[:,1],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.set_xlabel(r'$x_p $',labelpad = 5)
ax.set_ylabel(r'$y_p $',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

figname = fig_folder + 'Pixel_coord_' + img_name
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
#%% Set drone parameters

param_dict = {}
param_dict['H'] = 29.9
param_dict['alpha_0'] = 35.4
param_dict['focale'] = 2700

param_dict['latitude'] = 48.3463709995886
param_dict['longitude'] = -68.8267609065429
param_dict['azimuth'] = 14

#%% Georectify picture 

[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['H'],param_dict['alpha_0']*np.pi/180,
                                       param_dict['focale'])

# compute object position
POS_realx,POS_realy = dp.projection_real_space(POS[:,0],POS[:,1],x0,y0,
                                               param_dict['H'],param_dict['alpha_0']*np.pi/180,param_dict['focale'])

#%% Plot georectify image 

fig, ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
ax.plot(POS_realx,POS_realy,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

c.set_rasterized(True)
figname = fig_folder + 'Camera_coord_' + img_name
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Compute GPS position 

# horizontal distance between center of metric coordinates and drone position 
dist2drone = param_dict['H']/np.tan(param_dict['alpha_0']*np.pi/180)
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param_dict['latitude'],param_dict['longitude'],
                                                param_dict['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param_dict['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)

#GPS coordinates of detected objects 
POS_rho,POS_theta = dp.cart2pol(POS_realx,POS_realy)
POS_theta = POS_theta*180/np.pi
POS_azimuth = param_dict['azimuth'] + 90 - POS_theta
POS_Lat,POS_Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                POS_azimuth,POS_rho)

#%% Plot figure with GPS coordinates 

fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x


c.set_rasterized(True)
figname = fig_folder + 'GPS_coord_' + img_name
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')


#%% Save objects position in pixels, meters and GPS


POS_dict = {}
POS_dict['pixels'] = {}
POS_dict['pixels']['xpix'] = POS[:,0]
POS_dict['pixels']['ypix'] = POS[:,1]
POS_dict['meters'] = {}
POS_dict['meters']['X'] = POS_realx
POS_dict['meters']['Y'] = POS_realy
POS_dict['GPS'] = {}
POS_dict['GPS']['latitude'] = POS_Lat
POS_dict['GPS']['longitude'] = POS_Long

file2save = f'{fig_folder}{date}_{drone_ID}_{exp_ID}_{img_name}_phone_positions.h5'
rw.save_dict_to_h5(POS_dict,file2save)


#%% Load objects position 

file2load = f'{fig_folder}{date}_{drone_ID}_{exp_ID}_{img_name}_phone_positions.h5'
test = rw.load_dict_from_h5(file2load)

#%% Check bathymetry

path2bath_interp = 'U:/Data/Bathymetrie/linear_interpolator_bathymetry.pkl'
# load interpolator 
with open(path2bath_interp,'rb') as pf:
    interp_H = pickle.load(pf)
    

#%%

H = weather.get_bathymetry_GPS((POS_dict['GPS']['latitude'][0],POS_dict['GPS']['longitude'][0]), interp_H)

UTC_str = '2025-02-11T15-42-46'
format_date = '%Y-%m-%dT%H-%M-%S'
UTC0 = datetime.strptime(UTC_str,format_date)
tide = weather.tide_from_datetime(UTC0)





