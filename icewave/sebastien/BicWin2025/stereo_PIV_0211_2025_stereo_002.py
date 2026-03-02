# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 11:21:04 2026

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
from matplotlib.path import Path

import scipy.signal as signal
from scipy.interpolate import LinearNDInterpolator, griddata
from datetime import datetime,timedelta
import pytz

from shapely.geometry import Polygon, Point
from scipy.spatial import ConvexHull

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.tools.rw_data as rw 

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Function section 

def transpose_PIVmat_fields(m):
    """ Change dimensions of different fields computed using PIVlab """
    key_fields = ['Vx','Vy','Vz','X','Y']

    for key in key_fields:
        m[key] = np.transpose(m[key])
        
    return m

#%% Load data 

date = '0211'
base = 'U:/Data/'

drones = ['bernache','mesange']
data = {}
for i in range(len(drones)):
    path = f'{base}{date}/Drones/{drones[i]}/matData/*stereo_002/*scaled.mat'
    print(path) 
    filelist = glob.glob(path)    
    file2load = filelist[0]
    with h5py.File(file2load, 'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
    
        data[drones[i]] = mat2py.mat_to_dict(fmat['m'],fmat['m'])
        
#%% Transpose fields and supress quadratic noise 

for drone_key in data.keys():
    # data[drone_key] = transpose_PIVmat_fields(data[drone_key]) 

    for velocity_key in ['Vx','Vy']:
        V = data[drone_key][velocity_key]
        V = FT.supress_quadratic_noise(np.transpose(V,(1,0,2)),data[drone_key]['x'],data[drone_key]['y'])
        data[drone_key][velocity_key] = np.transpose(V,(1,0,2))

   
#%% Show velocity fields 

drone_ID = 'mesange'
frame = 80

fig, axs = plt.subplots(ncols = 2, figsize = (12,8), sharey = True)

imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.imshow(data[drone_ID]['Vx'][:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_x$')
        imsh[0].set_clim([-0.5,0.5])
    else:
        imsh.append(ax.imshow(data[drone_ID]['Vy'][:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_y$')
        imsh[1].set_clim([-0.5,0.5])
   
    ax.set_xlabel(r'$x_{pix}$')     

axs[0].set_ylabel(r'$y_{pix}$')
axs[0].set_title(f'{drone_ID} frame = {str(frame)}')
plt.tight_layout()

#%% georectified velocity fields 

drone_ID = 'mesange'
frame = 0
cmap = parula_map

Vx_values = [-0.5,0.5]
Vy_values = [-0.5,0.5]

fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,8))
imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.pcolormesh(data[drone_ID]['X'],data[drone_ID]['Y'],data[drone_ID]['Vx'][:,:,frame],
                                  shading = 'gouraud',cmap = cmap,
                                 vmin = Vx_values[0],vmax = Vx_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_X$')
    else:
        imsh.append(ax.pcolormesh(data[drone_ID]['X'],data[drone_ID]['Y'],data[drone_ID]['Vy'][:,:,frame],
                                  shading = 'gouraud',cmap = cmap,
                                 vmin = Vy_values[0],vmax = Vy_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_Y$')

    ax.set_aspect(1)
    
axs[0].set_ylabel(r'$Y$')
axs[0].set_title(f'{drone_ID} frame = {str(frame)}')
plt.tight_layout()

#%% Compute FK spectrum for a given drone 

# drone_ID = 'mesange'
# i_start = 0
# i_end = 4000

# Efk = FT.space_time_spectrum(data[drone_ID]['Vy'],data[drone_ID]['SCALE']['facq_x'],
#                              S['SCALE']['facq_t'],add_pow2 = [0,0,0])
