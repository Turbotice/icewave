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
import os 

from datetime import datetime
import pytz

import scipy

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.gps.gps_seb as gps_seb
import icewave.tools.weather as weather
import icewave.das.DAS_package as DS

parula_map = matcmaps.parula()
plt.rcParams.update({
    "text.usetex": True}) # use latex


#%% Function section 

def plot_spatio_temp(spatio,fiber_length,extents,cmap):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsh, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = cmap,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#%% Load DAS GPS position 

path2DAS_GPS = 'U:/Data/0211/DAS/fiber_GPS_bathymetry_0211.pkl'
with open(path2DAS_GPS,'rb') as pf:
    fiber = pickle.load(pf)
    
    
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
UTC_t = np.empty((len(filelist) - 1,int(file_duration/Nb_minutes)),dtype = object)

for k,file2load in enumerate(filelist[:-1]):
    
    _,_,UTC_stack,_ = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    tide = np.zeros((UTC_stack.shape[0],1))
    for i in range(UTC_stack.shape[0]):
        tide[i] = weather.tide_from_datetime(UTC_stack[i,0])
        UTC_t[k,i] = UTC_stack[i,0]

    # compute water height for each position of fiber
    for j in range(len(tide)):
        water_height[k,:,j] = fiber['H'] + tide[j]
        
#%% Reshape array of water height 

transpose_hw = np.transpose(water_height,(1,0,2))
shaped_hw = transpose_hw.reshape((transpose_hw.shape[0],transpose_hw.shape[1]*transpose_hw.shape[2]))

UTC_new = UTC_t.reshape(UTC_t.shape[0]*UTC_t.shape[1])

# create structure 
DAS_water_height = fiber
DAS_water_height['water_height'] = shaped_hw
DAS_water_height['UTC_t'] = UTC_new
DAS_water_height['units']['water_height'] = 'meter'

# save structure
file2save = f'{path2data}fiber_water_height_GPS_structure_{date}.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(DAS_water_height,pf)

#%% Plot spatio temporal of water height

set_graphs.set_matplotlib_param('single')
extents = [DAS_water_height['UTC_t'][0],DAS_water_height['UTC_t'][-1],
           DAS_water_height['s'][0],DAS_water_height['s'][-1]]
fig, ax, imsh, cbar = plot_spatio_temp(DAS_water_height['water_height'].T,fiber_length, extents,cmap = 'Blues')
cbar.set_label(r'$H \; \mathrm{(m)}$')

imsh.set_clim([2.0,5.5])

figname = f'{fig_folder}DAS_water_height_{date}'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')




#%% Compute water height using weather implemented functions 
   
interp_bathy = weather.get_bathy_interpolator()

test = np.zeros((len(filelist)-1,len(fiber['H']),int(file_duration/Nb_minutes)))
for k,file2load in enumerate(filelist[:-1]):
    
    _,_,UTC_stack,_ = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    for j in range(len(fiber['lat'])):
        for i in range(UTC_stack.shape[0]):
        
            GPS_coords = (fiber['lat'][j],fiber['long'][j])
            UTC_datetime = UTC_stack[i,0]
            test[k,j,i] = weather.get_water_height(GPS_coords, UTC_datetime,interp_bathy = interp_bathy)
           


#%% Check if weather.get_water_height function works

# transpose_hw = np.transpose(test,(1,0,2))
# weather_hw = transpose_hw.reshape((transpose_hw.shape[0],transpose_hw.shape[1]*transpose_hw.shape[2]))

fig, ax = plt.subplots()
ax.plot(UTC_new,shaped_hw[0,:],'o')
# ax.plot(UTC_new,weather_hw[0,:],'^')

