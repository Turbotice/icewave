# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:37:02 2026

@author: sebas
"""

import numpy as np
import os
import glob
import pandas as pd 
import h5py
import pickle
import cv2 as cv

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.signal as signal
from scipy.interpolate import LinearNDInterpolator, griddata, RegularGridInterpolator
from datetime import datetime,timedelta
import pytz

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.rw_data as rw
import icewave.gps.gps_seb as gps_seb
import icewave.gps.gps as gps
import icewave.geometry.tables as tables

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Import drone data

date = '0226'
base = f'K:/Share_hublot/Data/'

drone_ID = 'mesange'
exp_ID = '12-FRAC_001'
path2drone = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = glob.glob(f'{path2drone}*scaled.h5')[0]

S = rw.load_dict_from_h5(file2load)

#%% Define fig_folder

# fig_folder = f'{base}Summary/SWIIFT_article/'
# if not os.path.isdir(fig_folder):
#     os.mkdir(fig_folder)

fig_folder = f'F:/PhD_Manuscript/ch4/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Import image 

# Load image 
path2img = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/images/*0684*.tiff'
filelist_img = glob.glob(path2img)
file2load = filelist_img[0]

img = cv.imread(file2load)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

#%% Plot image in pixel coordinates 
# Set pixel coordinates of buoys 
POS = {}
POS['pix'] = np.array([[372.8,1179.7],[2173.7,1278.5],[3472.3,1434.6]])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.imshow(img)
ax.plot(POS['pix'][:,0],POS['pix'][:,1],'ro',mec = 'k')

#%% Georectify image and positions 

[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,S['DRONE']['h_drone'],S['DRONE']['alpha_0'],
                                       S['DRONE']['focale'])

# compute object position
POS_realx,POS_realy = dp.projection_real_space(POS['pix'][:,0],POS['pix'][:,1],x0,y0,
                                               S['DRONE']['h_drone'],S['DRONE']['alpha_0'],
                                                                                      S['DRONE']['focale'])

POS['real'] = np.zeros(POS['pix'].shape)
POS['real'][:,0] = POS_realx
POS['real'][:,1] = POS_realy

# flip image and inverse position of buoys
img = np.flip(img,axis = (0,1))
POS['real'][:,0] = - POS['real'][:,0] 

#%% PLot image in metric coordinates 

X0 = Xreal.min()
Y0 = Yreal.max()

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
# c = ax.pcolormesh(Xreal,Yreal,img[:,:,2],shading = 'auto', cmap = 'gray')
#c.set_rasterized(True)
extents = [Xreal[:,0].min() - X0,Xreal[:,-1].max() - X0,Y0 - Yreal[0,:].min(),Y0 - Yreal[-1,:].max()]
ax.imshow(img,extent = extents)
ax.plot(POS['real'][:,0] - X0,Y0 - POS['real'][:,1],'o',
        markerfacecolor = 'r',markeredgecolor = 'k',markersize = 7)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

figname = f'{fig_folder}situation_picture_XY_{drone_ID}_{exp_ID}_origin_lower_left'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute GPS coordinates 

# horizontal distance between center of metric coordinates and drone position 
dist2drone = S['DRONE']['h_drone']/np.tan(S['DRONE']['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(S['GPS']['latitude'],S['GPS']['longitude'],
                                                S['GPS']['azimuth'],dist2drone)

Lat,Long = dp.XY2GPS(Xreal,Yreal,Lat0,Long0,S['GPS']['azimuth'])

POS_Lat,POS_Long = dp.XY2GPS(POS['real'][:,0],POS['real'][:,1],Lat0,Long0,S['GPS']['azimuth'])
POS['GPS'] = np.zeros(POS['pix'].shape)
POS['GPS'][:,0] = POS_Lat
POS['GPS'][:,1] = POS_Long


#%% Plot image in GPS coordinates 

extents = [Long.min(),Long.max(),Lat.min(),Lat.max()]
fig, ax = plt.subplots()
# imsh = ax.imshow(img,extent = extents,origin = 'upper')
imsh = ax.pcolormesh(Long,Lat,img[:,:,2],shading = 'auto',cmap = 'gray')
ax.plot(POS['GPS'][:,1],POS['GPS'][:,0],'o',mfc = 'r',mec = 'k',ms = 7)
ax.set_aspect(1/np.cos(S['GPS']['latitude']*np.pi/180)) # scaling y/x

ax.set_xlim([Long.min(),Long.max()])
ax.set_ylim([Lat.min(),Lat.max()])

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')
imsh.set_rasterized(True)

figname = f'{fig_folder}situation_picture_GPS_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
