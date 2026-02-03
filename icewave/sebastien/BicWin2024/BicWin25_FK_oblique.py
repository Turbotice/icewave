# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 14:27:02 2026

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy
import cv2 as cv

import h5py
import pickle
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.datafolders as df
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.drone.drone_projection as dp
import icewave.tools.rw_data as rw
import icewave.sebastien.dispersion_relation_models as disp_rel

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Import data

base = 'D:/Rimouski_2025/Data/'
date = '0214'
drone_ID = 'mesange'
exp_ID = '05-waves_003'

fig_folder = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
filelist = glob.glob(f'{path2data}interpolated*quadratic*.h5')
print(filelist)

idx_file = 0
file2load = filelist[idx_file]
S = rw.load_dict_from_h5(file2load)

#%% Show interpolated field 

frame = 150
fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['grid_x'].T,S['grid_y'].T,S['interp_uz'][:,:,frame].T,shading = 'gouraud',
                     cmap = parula_map,
                    vmin = -0.04, vmax = 0.08)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label('$u_z$')
ax.set_aspect(1)

ax.set_xlabel('$X \; \mathrm{(m)}$')
ax.set_ylabel('$Y \; \mathrm{(m)}$')


#%% Compute FK spectrum 
# i_start = 0
# i_end = -900
Efk = FT.space_time_spectrum(S['interp_uz'],S['artificial_facq_x'],S['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Plot FK spectrum

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
Amin = 2e-6 # change to adjust colormap
Amax = 2e-3 # change to adjust colormap
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log',vmin = Amin, vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
kbounds = [0.05,1.5] # bounds for k axis on Efk plot
fbounds = [0.05,1] # bounds for f axis on Efk plot
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{u}_z| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# show gravity dispersion relation
hw = 4.0
k_th = np.linspace(0,2.0,100)
omega_grav = disp_rel.shallow_water(k_th,hw)
f_grav = omega_grav/2/np.pi
label_grav = '$h_w = ' + f'{hw:.1f}'+ '\; \mathrm{m}$'
ax.plot(k_th,f_grav,'r-',label = label_grav)

# show heavy hydroelastic 
h_ice = 0.3
E = 3e9
nu = 0.3
rho_ice = 917
f_elastic = disp_rel.heavy_hydroelastic(k_th,h_ice,hw,E,nu,rho_ice)/2/np.pi
label_elastic = '$(E,h_{ice}) = (' + f'{E/1e9:.1f},{h_ice:.2f}' + ')$'
ax.plot(k_th,f_elastic,'b-',label = label_elastic)
ax.legend()

figname = f'{fig_folder}FK_spectrum_interp_uz_withfit_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

