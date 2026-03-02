# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 17:20:59 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors 
import matplotlib.cm as cm

import h5py 
import glob
import os 
import pickle

from datetime import datetime
import pytz

import scipy
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pywt

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

# plt.rcParams.update({
#     "text.usetex": True}) # use latex

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% Load dispersion relation points 

date = '0211'
main_path = 'U:/'
path2data = f'{main_path}Data/{date}/DAS/'

swell_str = 'swell_corrected'

filepath = f'{main_path}/Data/{date}/DAS/Figures/Wavelet_study*/{date}_hydro_elastic_disp_relation_subpix_{swell_str}_file*.h5'
filelist = glob.glob(filepath, recursive = True)
print(filelist)

main_dict = {}
for file2load in filelist:
    current_disp = rw.load_dict_from_h5(file2load)
    UTC_0 = current_disp['UTC_0']
    main_dict[UTC_0] = current_disp
    
file2save = f'{path2data}disp_curv_from_CWT_{date}.h5'
rw.save_dict_to_h5(main_dict, file2save)

#%% Plot dispersion relation points for a given file

keys = list(main_dict.keys())
key = keys[6]
print(key)

freq_dict = main_dict[key]['freq']
k_dict = main_dict[key]['k']
x_dict = main_dict[key]['xpos']

# create a norm for x position
xvals = np.array(list(x_dict.values()))
norm = plt.Normalize(vmin=xvals.min(), vmax=xvals.max())
cmap = cm.viridis

fig, ax = plt.subplots()

for sub_key in freq_dict.keys():
    
    color = cmap(norm(x_dict[sub_key]))
    
    ax.plot(freq_dict[sub_key],k_dict[sub_key],color = color)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)

cbar.set_label(r'$x \; \mathrm{(m)}$')

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$')

#%% Plot dispersion relation points for a given position and different files


pos_key = list(main_dict[keys[6]]['xpos'].keys())

selected_idx = 0
norm = plt.Normalize(vmin = 0, vmax = len(main_dict.keys()))
cmap = new_blues

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
for i,key in enumerate(main_dict.keys()):
    
    freq = main_dict[key]['freq'][pos_key[selected_idx]]
    k = main_dict[key]['k'][pos_key[selected_idx]]
    
    color = cmap(norm(i))
    
    ax.plot(freq,k,'.',color = color)    
    
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$')  


#%% Average dispersion relation over all files 

# list of position keys
pos_key = list(main_dict[keys[0]]['xpos'].keys())
xvals = np.array(list(main_dict[keys[0]]['xpos'].values()))

results = {}
min_k_value = 0.02

# loop over all positions
for selected_idx in range(len(pos_key)):
    # loop over all files (different time)
    for i,key in enumerate(main_dict.keys()):
      
        freq = main_dict[key]['freq'][pos_key[selected_idx]].tolist()
        if i == 0:
            freq_list = freq
        
        else :
            for f in freq :
                if f not in freq_list:
                    freq_list.append(f)
        
    # sort list in ascending order
    freq_list = sorted(freq_list)
    
    # Build array of wavevectors
    k_array = np.zeros((len(freq_list),len(list(main_dict.keys()))))
    for i,freq in enumerate(freq_list):
        for j,key in enumerate(main_dict.keys()):
            
            current_freq_array = main_dict[key]['freq'][pos_key[selected_idx]]
            current_k_array = main_dict[key]['k'][pos_key[selected_idx]]
            
            mask_freq = np.where(current_freq_array == freq)[0]
            if len(mask_freq) > 0:
                k_array[i,j] = current_k_array[mask_freq][0]

    # Filter outliers
    mask = k_array < min_k_value
    k_array[mask] = np.nan
    
    # compute mean value and standard deviation 
    mean_k = np.nanmean(k_array,axis = 1)
    std_k = np.nanstd(k_array,axis = 1)
    
    freq_array = np.array(freq_list)
    results[pos_key[selected_idx]] = {'freq': freq_array, 
                                      'mean_k': mean_k,
                                      'std_k': std_k,
                                      'xpos': xvals[selected_idx]}

# save results dictionnary
if swell_str == '':
    file2save = f'{path2data}avg_disp_curv_from_CWT_{date}_no_swell_correction.h5'
else : 
    file2save = f'{path2data}avg_disp_curv_from_CWT_{date}_with_swell_correction.h5'
rw.save_dict_to_h5(results, file2save)


#%%
fig, ax = plt.subplots()
ax.plot(freq_list,mean_k, color = 'tab:blue')
ax.fill_between(freq_list,mean_k - std_k, mean_k + std_k, color = 'tab:blue',alpha = 0.2)

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$') 


