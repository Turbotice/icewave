# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 15:16:07 2025

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
import joblib

from datetime import datetime
import pytz

import scipy
import pywt

import sys
# sys.path.append('C:/Users/sebas/git')
# sys.path.append('/media/turbots/DATA/thiou/homes/skuchly/git/')

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))


global g
g = 9.81

#%% Collect data

main_path = 'U:/Data/'
date_DAS = ['0210', '0211', '0212']

date = '0211'
filepath = f'{main_path}{date}/DAS/Figures/Wavelet_study*/{date}_wavelet_flexural_modulus*.h5'
filelist = glob.glob(filepath, recursive = True)
print(filelist)

# define fig_folder
fig_folder = f'{main_path}{date}/DAS/Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Plot data for a single date

# open .h5 file 
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

# define norm
format_UTC = '%Y-%m-%dT%H-%M-%S'
data = rw.load_dict_from_h5(filelist[0])
t_epoch_min = datetime.strptime(data['UTC_0'], format_UTC).timestamp()
data = rw.load_dict_from_h5(filelist[-1])
t_epoch_max = datetime.strptime(data['UTC_0'], format_UTC).timestamp()

norm = colors.Normalize(vmin = t_epoch_min,vmax = t_epoch_max)
labels_UTC = []

for file2load in filelist:
    data = rw.load_dict_from_h5(file2load)
    UTC_0 = datetime.strptime(data['UTC_0'],format_UTC)
    t_epoch = UTC_0.timestamp()
    current_color = new_blues(norm(t_epoch))
    ax.errorbar(data['x'],data['D'],yerr = data['err_D'],fmt = '-o',color = current_color)
    
    current_label = UTC_0.strftime('%H-%M-%S')
    labels_UTC.append(current_label)
    
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
tick_locs = np.linspace(t_epoch_min,t_epoch_max,len(filelist))
cbar.set_ticks(tick_locs)
cbar.set_ticklabels([d for d in labels_UTC])
cbar.set_label(r'UTC time')
    
ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')

figname = f'{fig_folder}{date}_results_D_VS_x_evolution'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%% Superpose several dates

#%% Superpose D extraction from swell and active sources 















