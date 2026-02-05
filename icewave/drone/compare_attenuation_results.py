# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 15:18:02 2026

@author: sebas

This script aims at comparing results obtained from attenuation scripts : attenuation_ux.py, attenuation_uz.py
We would like to compare the attenuation coefficients obtained from ux, uz fields and compare
detection done on frequencies or wavevectors of the FK spectrum 

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy

import pickle
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.weather as weather
import icewave.drone.drone_projection as dp
import icewave.drone.drone_tools as drone_tools
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter
import icewave.tools.rw_data as rw
import icewave.drone.attenuation_module as att_mod

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')


#%% Load main_results

base = 'F:/Rimouski_2024/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '10-waves_005'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
suffixe = f'{date}_{drone_ID}_{exp_ID}'
filelist = glob.glob(f'{path2data}main_results_{suffixe}.pkl')
print(filelist)

file2load = filelist[0]
with open(file2load,'rb') as pf:
    main_results = pickle.load(pf)
    
#%% Compare dispersion relations 

fig, ax = plt.subplots()
u_keys = ['ux','uz']
for uk in u_keys:
    att_key = f'attenuation_{uk}'
    m = main_results[att_key]
    
    for key in ['time','space']:
        label = f'{uk}_{key}'
        ax.plot(m[key]['k'],m[key]['f'],'.',label = label)
        
ax.legend()
ax.set_xlabel('$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel('$f \; \mathrm{(Hz)}$')

#%% Compare attenuation laws

fig, ax = plt.subplots()
u_keys = ['ux','uz']
for uk in u_keys:
    att_key = f'attenuation_{uk}'
    m = main_results[att_key]
    
    for key in ['time','space']:
        label = f'{uk}_{key}'
        ax.plot(m[key]['f'],m[key]['alpha'],'.',label = label)
        
ax.legend()
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')

