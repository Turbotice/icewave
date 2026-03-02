# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 11:02:03 2025

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

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')


#%% Load picture 

main_path = 'K:/Share_hublot/Repository/'
date = '0306'
drone_ID = 'bernache'
exp_ID = '12-doc_vuedensemble'

path = f'{main_path}{date}/Drones/{drone_ID}/{exp_ID}/'
filename = 'DJI_20240307133501_0198_D_frame_100.tiff'
file2read = f'{path}{filename}'
    
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

POS = np.zeros((4,2))
POS[0,:] = [1237,924]
POS[1,:] = [1507,658]
POS[2,:] = [2466,671]
POS[3,:] = [2659,955]

# superpose objects on figure 
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.imshow(img)
ax.plot(POS[:,0],POS[:,1],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.set_xlabel(r'$x_p $',labelpad = 5)
ax.set_ylabel(r'$y_p $',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

#%% Set drone parameters

param_dict = {}
param_dict['H'] = 37.7
param_dict['alpha_0'] = 19.5
param_dict['focale'] = 2700

param_dict['latitude'] = 48.1564321605613
param_dict['longitude'] = -69.0185743749565
param_dict['azimuth'] = 358.5

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

# You may need to zoom in ! 
fig, ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
ax.plot(POS_realx,POS_realy,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

c.set_rasterized(True)


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

# You may need to zoom in 
fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x




