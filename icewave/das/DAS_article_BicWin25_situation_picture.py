# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:14:33 2025

@author: sebas
"""


import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
from scipy.io import loadmat

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather

# PARULA COLORMAP 
parula_map = matcmaps.parula()

#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/Situation_picture/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load DAS GPS position 

file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

#%% Situation picture 

main_path = 'U/Data/'
date = '0205'
drone_ID = 'mesange'
exp_ID = '18-doc_010'

path2img = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/'
filelist = glob.glob(f'{path2img}*.jpg')
file2read = filelist[0]

# fig_folder = f'U:/Data/{date}/Drones/{drone_ID}/Figures/{exp_ID}/'
# if not os.path.isdir(fig_folder):
#     os.mkdir(fig_folder)
    
img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

# img_name = filename.split('.')[0]

#%% Plot raw picture 

fig, ax = plt.subplots(figsize = (12,9))
ax.imshow(img)
ax.set_xlabel(r'$x_p$')
ax.set_ylabel(r'$y_p$')

# select pixels coord of beginning and end of fiber 
fiber_tips = np.array([[2395,2969],[2593,566]])
x_fiber = np.linspace(fiber_tips[0,0],fiber_tips[1,0],100)
y_fiber = np.linspace(fiber_tips[0,1],fiber_tips[1,1],100)

ax.plot(x_fiber,y_fiber,'-')

#%% Plot line with a given color


# create a parameter s along the line
s  = np.linspace(0,1,600)
x_fiber = (1-s)*fiber_tips[0,0] + s*fiber_tips[1,0]
y_fiber = (1-s)*fiber_tips[0,1] + s*fiber_tips[1,1]

# create line segments (pais of consecutive points)
points = np.array([x_fiber,y_fiber]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
# distance = np.sqrt((x_fiber - fiber_tips[0,0])**2 + (y_fiber - fiber_tips[0,1])**2)

# apparent distance
length_fiber = 600
app_distance = s*length_fiber

# create a norm and a line collection 
norm = colors.Normalize(vmin = app_distance.min(), vmax = app_distance.max())
lc = LineCollection(segments,cmap = 'viridis', norm = norm, linewidths = 4)
lc.set_array(app_distance) # values used for colormap

set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
ax.imshow(img)
ax.add_collection(lc)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$s \; \mathrm{(m)}$')

# hide labels
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# figname = f'{fig_folder}Situation_picture_0205_18_doc_010'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Try georectify the selected picture 

# Set drone parameters

param_dict = {}
param_dict['H'] = 110.95
param_dict['alpha_0'] = 22.4
param_dict['focale'] = 4100 # around 4100 (+- 100) for format 5280 x 2970 

param_dict['latitude'] = 48.3476967601294
param_dict['longitude'] = -68.8131401977678
param_dict['azimuth'] = 263


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
POS_realx,POS_realy = dp.projection_real_space(x_fiber,y_fiber,x0,y0,
                                               param_dict['H'],param_dict['alpha_0']*np.pi/180,param_dict['focale'])

#%% Plot georectified image 

fig, ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
# ax.plot(POS_realx,POS_realy,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

c.set_rasterized(True)
# figname = fig_folder + 'Camera_coord_situation_picture_0205_18_doc_010' 
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

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

#%% Plot figure using GPS coordinates of fiber 

fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
# ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x


c.set_rasterized(True)

# plot DAS GPS position 
ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'r.')


#%% Position fiber on image using DAS GPS coordinates 

# convert DAS GPS coordinates into (X,Y) coordinates 
X_DAS,Y_DAS = dp.GPS2XY(DAS_water_height['lat'], DAS_water_height['long'], Lat0, Long0, param_dict['azimuth'])

# interpolate new positions along the fiber 
s  = np.linspace(0,1,600)
X_interp = (1-s)*X_DAS[0] + s*X_DAS[-1]
Y_interp = (1-s)*Y_DAS[0] + s*Y_DAS[-1]

# create line segments (pais of consecutive points)
points = np.array([X_interp,Y_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp - X_DAS[0])**2 + (Y_interp - Y_DAS[0])**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = 'viridis', norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)
# ax.plot(X_DAS,Y_DAS,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

ax.add_collection(lc)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$s \; \mathrm{(m)}$')

ax.set_ylim([-240,600])
ax.set_xlim([-400,400])

figname = fig_folder + 'Camera_coord_situation_picture_0205_18_doc_010_DAS_GPS' 
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot initial figure with fiber position in pixel coordinates system 

xpix_interp,ypix_interp = dp.projection_pixel_space(X_interp, Y_interp, x0, y0, param_dict['H'], 
                                               param_dict['alpha_0']*np.pi/180, param_dict['focale'])


# create line segments (pairs of consecutive points)
points = np.array([xpix_interp,ypix_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp - X_DAS[0])**2 + (Y_interp - Y_DAS[0])**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = 'viridis', norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
c = ax.imshow(img)
ax.set_xlabel(r'$x_{pix}$',labelpad = 5)
ax.set_ylabel(r'$y_{pix}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

ax.add_collection(lc)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$s \; \mathrm{(m)}$')

figname = f'{fig_folder}Pixel_coords_situation_picture_0205_18_doc_010_DAS_GPS'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Create figure for presentations

xpix_interp,ypix_interp = dp.projection_pixel_space(X_interp, Y_interp, x0, y0, param_dict['H'], 
                                               param_dict['alpha_0']*np.pi/180, param_dict['focale'])


# create line segments (pairs of consecutive points)
points = np.array([xpix_interp,ypix_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp - X_DAS[0])**2 + (Y_interp - Y_DAS[0])**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
c = ax.imshow(img)
# ax.set_xlabel(r'$x_{pix}$',labelpad = 5)
# ax.set_ylabel(r'$y_{pix}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

ax.add_collection(lc)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$x \; \mathrm{(m)}$')

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

figname = f'{fig_folder}Pixel_coords_situation_picture_0205_18_doc_010_DAS_GPS_powerpoint_version'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)