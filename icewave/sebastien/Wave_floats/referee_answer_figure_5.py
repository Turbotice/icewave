# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 11:13:47 2025

@author: sebas
"""

import numpy as np
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
import pickle
import scipy.signal
from scipy.signal import find_peaks


# import modules 
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
parula_map = matcmaps.parula()

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/Inkscape_drawings/figure5_referees_corrected/'

font_size_medium = 26
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (8,6)

#%% Load data of x-profile

path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/Inkscape_drawings/figure5_referees_corrected/'
file2load = f'{path2data}profile_X_data.mat'

profile_x = {}
with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    profile_x = mat2py.mat_to_dict(fmat['profile_x'],fmat['profile_x'])
    
#%% Load data of y-profile

file2load = f'{path2data}profile_Y_data.mat'

profile_y = {}
with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    profile_y = mat2py.mat_to_dict(fmat['profile_y'],fmat['profile_y'])


#%% Plot both profiles 

fig, ax = plt.subplots()
ax.plot(profile_x['demod']['x'],profile_x['demod']['y'])
ax.plot(profile_x['envelop']['x'],profile_x['envelop']['y'],'.')
ax.set_ylim([-0.15,0.15])
ax.set_xlim([profile_x['demod']['x'].min(), profile_x['demod']['x'].max()])

fig, ax = plt.subplots()
ax.plot(profile_y['demod']['x'],profile_y['demod']['y'])
ax.plot(profile_y['envelop']['x'],profile_y['envelop']['y'],'.')



#%% Change y coordinates for y profile 

# y = profile_y['demod']['x'] - profile_y['demod']['x'].min()

# range of over which profile X is plotted 
range_x = profile_x['demod']['x'].max() - profile_x['demod']['x'].min()
print(range_x)


y = profile_y['demod']['x']

fig, ax = plt.subplots()
ax.plot(y, profile_y['demod']['y'])
ax.set_xlim([-range_x/2,range_x/2])
ax.set_ylim([-0.15,0.15])

#%% Finale figures 
mk_size = 1
# plot profile x 
# set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots(figsize = fig_size)
ax.plot(profile_x['demod']['x'],profile_x['demod']['y'])
ax.plot(profile_x['envelop']['x'],profile_x['envelop']['y'],'r.',ms = mk_size)
ax.set_ylim([-0.15,0.15])
ax.set_xlim([profile_x['demod']['x'].min(), profile_x['demod']['x'].max()])

ax.set_xlabel(r'$x \; \mathrm{(cm)}$')
ax.set_ylabel(r'$\tilde{\eta}^f_{||} \; \mathrm{(cm)}$')

figname = 'new_profile_x'
figname = f'{fig_folder}{figname}'

plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')


# plot profile y 
range_x = profile_x['demod']['x'].max() - profile_x['demod']['x'].min()
print(range_x)

fig, ax = plt.subplots(figsize = fig_size)
ax.plot(profile_y['demod']['x'], profile_y['demod']['y'])
ax.plot(profile_y['envelop']['x'],profile_y['envelop']['y'],'r.',ms = mk_size)
ax.set_xlim([-range_x/2,range_x/2])
ax.set_ylim([-0.15,0.15])
ax.set_xlabel(r'$x \; \mathrm{(cm)}$')
ax.set_ylabel(r'$\tilde{\eta}^f_{\perp} \; \mathrm{(cm)}$')

figname = 'new_profile_y'
figname = f'{fig_folder}{figname}'

plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')







