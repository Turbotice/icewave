# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:57:05 2025

@author: sebas

Compute water height below the optic fiber. 

GPS position of the optic fiber is computed from script DAS_GPS_bathymetry.py which uses 
waypoints of active sources achieved on Feburary 11th 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob 
import pickle 
import h5py 
import csv
import os 

import gpxpy
import cartopy.crs as ccrs
import re 

from datetime import datetime
import pytz

import scipy

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.gps.gps_seb as gps_seb
import icewave.tools.weather as weather
import icewave.das.DAS_package as DS
import icewave.drone.drone_projection as dp

parula_map = matcmaps.parula()
plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Load DAS GPS position 

path2DAS_GPS = 'U:/Data/0211/DAS/fiber_GPS_bathymetry_0211.pkl'
with open(path2DAS_GPS,'rb') as pf:
    fiber = pickle.load(pf)
    
#%% Function section 

def plot_spatio_temp(spatio,t,s,fiber_length,extents):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - t, numpy array or list, time array 
            - s, numpy array or list, curvilinear coordinate array
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsho, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

    
#%% Load DAS data

date = '0211'
# fig_folder = f'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2025/DAS/{date}/Figures/'

# Load parameters for DAS
path2DAS_param = 'U:/Data/parameters_Febus_2025.pkl'
with open(path2DAS_param,'rb') as pf:
    param = pickle.load(pf)
print('Parameters file loaded')

# Set parameters
# fs = np.shape(strain_rate)[1] # time sampling frequency 
# facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
fs = param[date]['fs']
fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
facq_x = param[date]['facq_x'] 

path2data = f'U:/Data/{date}/DAS/'

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(path2data + '*.h5')
Nb_minutes = 1 # duration of each stack
file_duration = 10 # duration of a file in minute 

water_height = np.zeros((len(filelist)-1,len(fiber['H']),int(file_duration/Nb_minutes)))
for k,file2load in enumerate(filelist[:-1]):
    
    _,_,UTC_stack,_ = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    tide = np.zeros((UTC_stack.shape[0],1))
    for i in range(UTC_stack.shape[0]):
        tide[i] = weather.tide_from_datetime(UTC_stack[i,0])

    # compute water height for each position of fiber
    for j in range(len(tide)):
        water_height[k,:,j] = fiber['H'] + tide[j]

#%% Reshape array of water height 

transpose_hw = np.transpose(water_height,(1,0,2))
test = transpose_hw.reshape((transpose_hw.shape[0],transpose_hw.shape[1]*transpose_hw.shape[2]))

#%%

fig, ax = plt.subplots()
ax.plot(test[0,:],'o')

for k in range(water_height.shape[0]):
    ax.plot(water_height[k,0,:])





#%% Plot spatio 

# set_graphs.set_matplotlib_param('single')
# chunk = 0
# extents = [UTC_stack[chunk,0], UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax , imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:],UTC_stack[chunk,:],s,fiber_length,extents = extents)

# imsh.set(clim = [-1e4, 1e4])

