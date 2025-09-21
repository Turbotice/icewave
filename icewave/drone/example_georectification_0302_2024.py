# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:02:31 2024

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

# one may want to comment these two lines
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Set path to image, object pixel position and drone parameters 

# one may want to modify this path to the directory where georectification images are stored
main_path = 'K:/Share_hublot/Repository/georectification_example/'
date = '0302'
drone_ID = 'bernache'
exp_ID = '06-Calib_oblique_002'

path = f'{main_path}{date}/Drones/{drone_ID}/{exp_ID}/'
filelist_tif = glob.glob(path + '*.tif')
filelist_param = glob.glob(path + 'Param*.mat')
filelist_POS = glob.glob(path + 'POS*.mat')


#%% load image

i0 = 0
frame = filelist_tif[i0]
img_name = frame.split('\\')[-1]
img_name = img_name.split('.tif')[0]

img = cv.imread(frame)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)


#%% Load the .mat file of parameters and convert to a dictionary
 
param_file = filelist_param[i0]

with h5py.File(param_file, 'r') as fmat:
    param_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    param_dict = mat2py.mat_to_dict(fmat['param'],fmat['param'])

# Load the .mat file of objects position (without h5py)

POS_file = filelist_POS[i0]
POS_matdict = loadmat(POS_file)
POS = POS_matdict['POS']

#%% Plot raw image and detected buoys 

fig, ax = plt.subplots()
ax.imshow(img)
ax.plot(POS[:,1],POS[:,2],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.set_xlabel(r'$x_p $',labelpad = 5)
ax.set_ylabel(r'$y_p $',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 


#%% Georectify image and compute object position
    
focale = 2700
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['H'],param_dict['alpha_0']*np.pi/180,focale)

# compute object position
POS_realx,POS_realy = dp.projection_real_space(POS[:,1],POS[:,2],x0,y0,
                                               param_dict['H'],param_dict['alpha_0']*np.pi/180,focale)


#%% Display georectified image 

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

fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x


c.set_rasterized(True)


