# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 16:48:15 2025

@author: sebas

Synchronization of buoys and drones for 0226/2024

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

import scipy
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

#%% Load data

# base = 'K:/Share_hublot/Data/'
base = 'F:/Rimouski_2024/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '12-FRAC_001'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'

fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

# load data 
file2load = f'{path2data}data_buoys_velocity_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2load,'rb') as pf:
    data = pickle.load(pf)
# data = rw.load_dict_from_h5(file2load)

#%% Create timestamp arrays
N = len(data['drone']['B4'][0,:])
t_drone = np.array([UTC.timestamp() for UTC in data['UTC']['drone'][:N]])
t_buoy = {}
for key_buoy in data['UTC']['buoy'].keys():
    t_buoy[key_buoy] = {}
    for key_acq in data['UTC']['buoy'][key_buoy].keys():
        t_buoy[key_buoy][key_acq] = np.array([UTC.timestamp() for UTC in data['UTC']['buoy'][key_buoy][key_acq]])


#%% Filter signals

# filter buoys signal
fs = 50 # sampling frequency
f_high = 2 # cutoff frequency
f_low = 0.1 
order_filter = 4
b,a = signal.butter(order_filter,[f_low,f_high],'bandpass',fs = fs)

for key_buoy in data['buoy'].keys():
    for key_acq in data['buoy'][key_buoy].keys():
        data['buoy'][key_buoy][key_acq] = signal.filtfilt(b,a,data['buoy'][key_buoy][key_acq],axis = -1)

#%% Plot buoys signals for each acquisition

idx_u = 2 # 0 = u_x, 1 = u_y, 2 = u_z

color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

# outlier = 'B1_1900'
outlier = 'B1acq_1'

fig, axs = plt.subplots(figsize = (12,9),ncols = 1, nrows = 4,constrained_layout = True)
for j,key_acq in enumerate(data['buoy']['B1'].keys()):
    for i,key_buoy in enumerate(data['drone'].keys()):
        current_color = color_dict[key_buoy]

        test = key_buoy + key_acq
        print(test)
        if test != outlier:
            axs[j].plot(data['UTC']['buoy'][key_buoy][key_acq],data['buoy'][key_buoy][key_acq][idx_u,:],
                                '-',color = current_color,label = key_buoy)
        if j == 0:
            axs[j].legend(ncols = 3, bbox_to_anchor = (0.11, 1),
                          loc='lower left')
            
    axs[j].set_ylabel(r'$u_z$')
    axs[j].set_title(f'Acquisition = {key_acq}')

axs[0].legend(ncols = 3, bbox_to_anchor = (0.11, 1),
              loc='lower left')

figname = f'{fig_folder}Comparison_buoys_acquisition'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compare Drone/buoys signals

key_acq = 'acq_1'
color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

N = len(data['drone']['B4'][0,:])
fig, axs = plt.subplots(nrows = 3,ncols = 1)
for i,key_buoy in enumerate(data['drone'].keys()):
    print(key_buoy)
    axs[i].plot(data['UTC']['drone'][:N],data['drone'][key_buoy][1,:],color = color_dict[key_buoy])
    axs[i].plot(data['UTC']['buoy'][key_buoy][key_acq],data['buoy'][key_buoy][key_acq][2,:],'k-')
    axs[i].set_ylim([-0.5,0.5])

axs[0].set_title(key_acq)

#%% Interpolate each acquisition with a sampling frequency of 30 Hz

fs_in = 50
fs_out = 30
down_buoys = {}
down_t = {}
for key_buoy in data['buoy'].keys():
    down_buoys[key_buoy] = {}
    down_t[key_buoy] = {}
    for key_acq in data['buoy'][key_buoy].keys():
        duration = data['buoy'][key_buoy][key_acq].shape[1] / fs_in
        N_out = int(duration * fs_out)
        down_buoys[key_buoy][key_acq] = scipy.signal.resample(data['buoy'][key_buoy][key_acq],N_out,axis = -1)
        down_t[key_buoy][key_acq] = scipy.signal.resample(t_buoy[key_buoy][key_acq],N_out)

#%% Check no phase shift

key_acq = 'acq_0'
color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

N = len(data['drone']['B4'][0,:])
fig, axs = plt.subplots(nrows = 3,ncols = 1)
for i,key_buoy in enumerate(data['drone'].keys()):
    axs[i].plot(down_t[key_buoy][key_acq],down_buoys[key_buoy][key_acq][2,:],color = color_dict[key_buoy])
    axs[i].plot(t_buoy[key_buoy][key_acq],data['buoy'][key_buoy][key_acq][2,:],'k-')
    axs[i].set_ylim([-0.5,0.5])

axs[0].set_title(key_acq)

#%% Compute correlations between drone signal and buoys acquisitions

key_buoy = 'B2'
key_acq = 'acq_1'
C = np.correlate(data['drone'][key_buoy][1,:],down_buoys[key_buoy][key_acq][2,:])

fig, ax = plt.subplots()
ax.plot(C)

# get maximum of correlation
ilag_max = np.argmax(abs(C))
dtlag_max = ilag_max / fs_out
delta_t = timedelta(seconds = dtlag_max)
print(dtlag_max)

#%% Apply time shift 

key_acq = 'acq_1'
color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

N = len(data['drone']['B4'][0,:])
shift_hand = 397.1 - 1.48 # adjusted by hand 
total_shift = -dtlag_max + shift_hand
print(total_shift)
fig, ax = plt.subplots(figsize = (12,6))
key_buoy = 'B2'

ax.plot(t_drone + total_shift,data['drone'][key_buoy][1,:],color = color_dict[key_buoy],label = 'drone')
ax.plot(down_t[key_buoy][key_acq],down_buoys[key_buoy][key_acq][2,:],'k-',lw = 0.5,label = 'buoy')
ax.set_ylim([-1,1])

ax.set_title(f'{key_buoy} - {key_acq}')
ax.set_xlabel(r'$t_{epoch} \; \mathrm{(s)}$')
ax.set_ylabel(r'$u_z \; \mathrm{(m.s^{-1})}$')
ax.legend()

# figname = f'{fig_folder}synchro_buoy_drone_{key_buoy}_{key_acq}_all_length'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Zoom in 

fig, ax = plt.subplots(figsize = (12,6))
key_buoy = 'B2'

ax.plot(t_drone + total_shift,data['drone'][key_buoy][0,:],color = color_dict[key_buoy],label = 'drone')
# ax.plot(down_t[key_buoy][key_acq],down_buoys[key_buoy][key_acq][2,:],'k-',lw = 0.5,label = 'buoy')
ax.plot(t_buoy[key_buoy][key_acq],data['buoy'][key_buoy][key_acq][0,:],'k-',lw = 0.5,label = 'buoy')
ax.set_ylim([-0.25,0.25])
ax.set_xlim([1708975160,1708975480])

ax.set_title(f'{key_buoy} - {key_acq}')
ax.set_xlabel(r'$t_{epoch} \; \mathrm{(s)}$')
ax.set_ylabel(r'$u_z \; \mathrm{(m.s^{-1})}$')
ax.legend()

# figname = f'{fig_folder}synchro_buoy_drone_{key_buoy}_{key_acq}_zoomed'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Apply time shift to all buoys, using 
key_acq = 'acq_1'
key_buoy = 'B1'

fig, ax = plt.subplots()
ax.plot(data['buoy'][key_buoy][key_acq][2,:])





#%%
key_acq = 'acq_0'
fig, axs = plt.subplots(nrows = 3,ncols = 1)
for i,key_buoy in enumerate(data['drone'].keys()):
    axs[i].plot(down_buoys[key_buoy][key_acq][2,:],color = color_dict[key_buoy])
    axs[i].plot(data['buoy'][key_buoy][key_acq][2,:],'k-')
    axs[i].set_ylim([-0.5,0.5])









#%% Concatenate all acquisitions

u_buoy = {}
UTC_buoy = {}
for key_buoy in data['buoy'].keys():
    u_buoy[key_buoy] = np.concatenate([data['buoy'][key_buoy][key_acq] 
                                       for key_acq in data['buoy'][key_buoy].keys()],axis = -1)    

    UTC_buoy[key_buoy] = np.concatenate([data['UTC']['buoy'][key_buoy][key_acq] 
                                       for key_acq in data['UTC']['buoy'][key_buoy].keys()],axis = -1)
    
#%%    

color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

N = len(data['drone']['B4'][0,:])
fig, axs = plt.subplots(nrows = 3,ncols = 1)
for i,key_buoy in enumerate(data['drone'].keys()):
    axs[i].plot(data['UTC']['drone'][:N],data['drone'][key_buoy][1,:],color = color_dict[key_buoy])
    axs[i].plot(UTC_buoy[key_buoy],u_buoy[key_buoy][2,:],'k-')

#%% Correlations  











#%% Load time shift between buoys and drones from '0211/2024'


file2load = 'K:/Share_hublot/Data/0211/Drones/Buoy_drone_synchro_time_delta_i0_520.pkl'

with open(file2load,'rb') as pf:
    synchro_buoy_drone = pickle.load(pf)

UTC_drone = data['UTC']['drone']
UTC_drone = np.array(UTC_drone) + synchro_buoy_drone['shift_timedelta']

#%%

key_acq = 'acq_0'
color_dict = {'B4':'tab:orange','B2':'tab:green','B1':'tab:blue'}

N = len(data['drone']['B4'][0,:])
fig, axs = plt.subplots(nrows = 3,ncols = 1)
for i,key_buoy in enumerate(data['drone'].keys()):
    axs[i].plot(UTC_drone[:N],data['drone'][key_buoy][1,:],color = color_dict[key_buoy])
    axs[i].plot(data['UTC']['buoy'][key_buoy][key_acq],data['buoy'][key_buoy][key_acq][2,:],'k-')









