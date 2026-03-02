# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 16:32:52 2025

@author: sebas
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import h5py
import glob

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Import data

# one may want to modify the base path
base = 'K:/Share_hublot/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
filelist = glob.glob(f'{path2data}*scaled.mat')
print(filelist)

idx_file = 0
file2load = filelist[idx_file]

# load file 
with h5py.File(file2load, 'r') as fmat:
    S = {}

    print('Top-level keys : ', list(fmat.keys()))

    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    S = mat2py.transpose_PIVmat_fields(S)


#%% Supress quadratic noise

Vx = FT.supress_quadratic_noise(np.transpose(S['Vx'],(1,0,2)),S['x'],S['y'])
Vy = FT.supress_quadratic_noise(np.transpose(S['Vy'],(1,0,2)),S['x'],S['y'])
Vx = np.transpose(Vx,(1,0,2))
Vy = np.transpose(Vy,(1,0,2))

print('Quadratic field supressed')

#%% Plot velocity field 

frame = 0
extents_meter = np.array([S['x'].min(),S['x'].max(),S['y'].min(),S['y'].max()])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(Vx[:,:,frame].T,cmap = parula_map,origin = 'lower',aspect = 'equal',extent = extents_meter)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_x \; \mathrm{(m.s^{-1})}$')
imsh.set_clim([-2,2])

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$y \; \mathrm{(m)}$')

#%% Compute Fourier spectrum of the velocity field 

fps = S['SCALE']['facq_t']
TF_spectrum,freq_time,FFT_t = FT.temporal_FFT(Vx,fps,padding_bool = 1,add_pow2 = 0,output_FFT = True)

#%%
fig, ax = plt.subplots()
ax.loglog(freq_time,abs(TF_spectrum),'-',zorder = 2)
ax.set_ylim(1e-5,1e-1)
ax.grid(zorder = 1)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle_{x,y} \; \mathrm{(u.a.)}$')


#%% Show demodulated field for a given frequency 

f_demod = 0.212 # frequency chosen for demodulation 
idx = np.argmin(abs(freq_time - f_demod))
field = np.real(FFT_t[:,:,idx])

fig, ax = plt.subplots()
c = ax.imshow(field.T, cmap = parula_map, aspect = 'equal', origin = 'lower', 
              interpolation = 'gaussian', extent = extents_meter)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)
cbar.set_label(r'$\hat{V}_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)


#%% Compute space-time spectrum 
    
N = Vx.shape[2]
Efk = FT.space_time_spectrum(Vx,1/S['SCALE']['fx'],S['SCALE']['facq_t'],add_pow2 = [0,0,0])


#%% Plot space time spectrum 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

c = ax.imshow(Efk['E'].T, cmap = parula_map , aspect = 'auto', norm = 'log', 
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['f'].min(),Efk['f'].max(),Efk['k'].min(),Efk['k'].max()))

ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

ax.set_xlim([0.08,1.5])
ax.set_ylim([0.07,4])
ax.set_xscale('log')
ax.set_yscale('log')
c.set_clim([2e-5,6e-3])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$|\hat{V}_x| (k,f) \; \mathrm{(u.a.)}$',labelpad = 5)
