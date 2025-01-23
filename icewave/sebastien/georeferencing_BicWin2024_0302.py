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

#%%

base = 'Y:/Banquise/Sebastien/Calibration_intrinseque/Calib_oblique_bernache_2/'
filelist_tif = glob.glob(base + '*.tif')
filelist_param = glob.glob(base + 'Param*.mat')
filelist_POS = glob.glob(base  + 'POS*.mat')


fig_folder = 'C:/Users/sebas/OneDrive/Bureau/Article_BicWin2024/' + 'Figures_python_small_size/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%%

# Parameters for plots
font_size_medium = 20
font_size_small = round(0.6*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 100 # dpi to save images 

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%%

i0 = 3
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

figname = fig_folder + 'Pixel_coord_' + img_name
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Georeference image and compute object position
    
focale = 2700
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['H'],param_dict['alpha_0']*np.pi/180,focale)

# compute object position
POS_realx,POS_realy = dp.projection_real_space(POS[:,1],POS[:,2],x0,y0,
                                               param_dict['H'],param_dict['alpha_0']*np.pi/180,focale)


#%% Display georeferenced image 

fig, ax = plt.subplots(figsize = fig_size)
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

fig,ax = plt.subplots(figsize = fig_size)
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x


c.set_rasterized(True)
figname = fig_folder + 'GPS_coord_' + img_name
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')













#%% Plot figure of tiles used for calibration 

base_tile = 'Y:/Banquise/Sebastien/Article_BicWin2024/Figures_calibration_georef/'
path2image = base_tile + 'im_0953.tiff'

img_tile = cv.imread(path2image);
img_tile = cv.cvtColor(img_tile,cv.COLOR_BGR2RGB)

ny,nx,nc = np.shape(img_tile)


#%%
x_pix = np.arange(nx)
x_pix = np.arange(ny)

ymin = 80
ymax = ny
xmin = 1100
xmax = 2590

cut_img = img_tile[ymin:ymax,xmin:xmax,:]
step_ticks = 400
# labels of x and y ticks
xticklabel = np.arange(xmin,xmax,step_ticks)
yticklabel = np.arange(ymin,ymax,step_ticks)

# position of x and y ticks 
xtickposition = np.arange(np.shape(cut_img)[1],step = step_ticks)
ytickposition = np.arange(np.shape(cut_img)[0],step = step_ticks)

#%%
xlabel = ['']*np.size(xticklabel)
for i,label in enumerate(xticklabel) :
    label = '$' + str(label) + '$'
    xlabel[i] = label

ylabel = ['']*np.size(yticklabel)
for i,label in enumerate(yticklabel) :
    label = '$' + str(label) + '$'
    ylabel[i] = label

#%%

fig,ax = plt.subplots(figsize = fig_size)
ax.imshow(cut_img)
ax.set_xticks(xtickposition,labels = xlabel)
ax.set_yticks(ytickposition,labels = ylabel)
ax.set_xlabel(r'$x_p $',labelpad = 15)
ax.set_ylabel(r'$y_p $',labelpad = 15)
ax.set_aspect(1) # set aspect ratio to 1 

figname = fig_folder + 'Tile_image'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')