# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 11:14:23 2026

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 

import h5py

import scipy 

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs
import icewave.das.DAS_package as DS
import icewave.tools.rw_data as rw

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Set fig_folder path 

fig_folder = 'E:/Rimouski_2025/DAS_article/Figures_answer_reviews/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)    

#%% Load parameters for DAS
main_path = 'E:/Rimouski_2025/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

date = '0211'
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 5 #5 for 0211
file2load = filelist[idx_file]
print(file2load)

Nb_minutes = 10 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
format_date = '%Y-%m-%dT%H-%M-%S'
label_UTC0 = UTC_stack[0,0].strftime(format_date)

# shift curvilinear axis
offset_fiber = 37.5
s = s - offset_fiber

#%% Show spatio-temporal 
 
chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,6))
imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax.set_ylim([0,s[-1]])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'UTC')

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

#%% Compute frequency spectrum 
chunk = 0
TF_spectrum,freq = FT.temporal_FFT_spatio(stack_strain[chunk,:,:].T,fs)

# save spectrum in a dictionnary
my_dict = {}
my_dict[date] = {'TF_spectrum':TF_spectrum,'f':freq}

# =============================================================================
# %% Load data from 0212
# =============================================================================

date = '0212'
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 0 #5 for 0211
file2load = filelist[idx_file]
print(file2load)

Nb_minutes = 10 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
format_date = '%Y-%m-%dT%H-%M-%S'
label_UTC0 = UTC_stack[0,0].strftime(format_date)

# shift curvilinear axis
offset_fiber = 37.5
s = s - offset_fiber

#%% Show spatio-temporal 
chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,6))
imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax.set_ylim([0,s[-1]])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'UTC')

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

#%% Compute frequency spectrum 
chunk = 0
TF_spectrum,freq = FT.temporal_FFT_spatio(stack_strain[chunk,:,:].T,fs,padding_bool = 1)

# save spectrum in a dictionnary
my_dict[date] = {'TF_spectrum':TF_spectrum,'f':freq}

# =============================================================================
# %% Final plot 
# =============================================================================

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for date in my_dict.keys():
    TF_spectrum = my_dict[date]['TF_spectrum']
    freq = my_dict[date]['f']
    
    ax.plot(freq,TF_spectrum,label = date)
    
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylim([1e1,1e3])
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\hat{\dot{\epsilon}} \; \mathrm{(a.u.)}$')
ax.legend()

figname = f'{fig_folder}frequency_spectrum_waves_signal'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

