# -*- coding: utf-8 -*-
"""
Created on Mon May  5 16:21:49 2025

@author: sebas
"""

import numpy as np
import os
import glob  
import h5py
import pickle
import scipy
import math
import cv2 as cv
# from skimage.filters import meijering, sato, frangi, hessian

import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
from matplotlib.path import Path
import matplotlib.animation as animation

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter
import icewave.tools.datafolders as df

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% Load data for a given date and experiment
date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

main_path = df.find_path(disk = 'Elements',year = '2024')
path2data = f'{main_path}{date}/Drones/{drone_ID}/matData/{exp_ID}/'

filelist = glob.glob(f'{path2data}*scaled.mat')
print(f'Files available : {filelist}')
file2load = filelist[0]
with h5py.File(file2load,'r') as fmat:
    print('Top-level keys : ', list(fmat.keys()))
    data = mat2py.mat_to_dict(fmat['m'],fmat['m'])

data = mat2py.transpose_PIVmat_fields(data)

#%% Reduce drone noise and compute FFT spectrum 
Vs = FT.supress_quadratic_noise(np.transpose(data['Vx'],(1,0,2)),data['x'],data['y'])
Vs = np.transpose(Vs,(1,0,2))

#%% Load image associated to this movie
img_list = glob.glob(f'{main_path}{date}/Drones/{drone_ID}/{exp_ID}/*exemple.tiff')
file2load = img_list[0]

img = cv.imread(file2load)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

fig, ax = plt.subplots()
ax.imshow(img)

xmin = 0
xmax = img.shape[1] - 1

crop = img[:,xmin:xmax,:]
fig, ax = plt.subplots()
ax.imshow(crop)

# frasil_boxes = np.where(np.logical_and(data['PIXEL']['x_pix'] < xmax,data['PIXEL']['x_pix'] > xmin))[0]
# idx_start = frasil_boxes[0]
# idx_end = frasil_boxes[-1]  
  
#%% Define fig_folder and save folder
Dt = int(data['PIV_param']['Dt'])
fig_folder = f'{path2data}Plots/attenuation_from_space_domain_Dt{Dt}/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

save_folder =f'{path2data}Results/'
if not os.path.isdir(save_folder):
    os.mkdir(save_folder)

#%% Compute histogram
figname = f'{fig_folder}Histogram_{date}_{drone_ID}_{exp_ID}'
FT.histogram_PIV(data['Vx'],data['PIV_param']['w'],figname)

#%% Compute FFT spectrum 
TF_spectrum,freq,FFT_t = FT.temporal_FFT(Vs,data['SCALE']['facq_t'],padding_bool = 1,add_pow2 = 0,output_FFT = True)

# save data in a structure
FFT_spectrum = {}
FFT_spectrum['TF_spectrum'] = TF_spectrum
FFT_spectrum['f'] = freq

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle _{x,y}(f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylim([1e-5,1e0])

figname = f'{fig_folder}TF_spectrum_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute space FFT for a given frequency and extract main peaks 





















#%% Compute space FFT for a given frequency 

select_freq = 0.4
idx_freq = np.argmin(abs(select_freq - freq))
field = FFT_t[:,:,idx_freq]

fig, ax = plt.subplots(figsize = (12,9))
c = ax.imshow(np.real(field).T,cmap = parula_map,aspect = 'equal',norm = 'linear',origin = 'lower',interpolation = 'gaussian',
              extent = (data['x'].min(),data['x'].max(),data['y'].min(),data['y'].max()))

# ax.plot(kx[unravel_coords[0]],ky[unravel_coords[1]],'ro')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\hat{V_x} (x,y) \; \mathrm{(u.a.)}$',labelpad = 1)
ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$k_y \; \mathrm{(rad.m^{-1})}$', labelpad = 5)

# Spatial FFT
facq = [1/data['SCALE']['fx'],1/data['SCALE']['fx']]
shift,kx,ky = FT.fft_2D(field,facq)

fig, ax = plt.subplots(figsize = (12,9))
c = ax.imshow(abs(shift).T,cmap = parula_map,aspect = 'equal',norm = 'log',origin = 'lower')

              # extent = (kx.min(),kx.max(),ky.min(),ky.max()))
c.set_clim([1e-4,2e-2])
# ax.plot(kx[unravel_coords[0]],ky[unravel_coords[1]],'ro')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\hat{V_x} (x,y) \; \mathrm{(u.a.)}$',labelpad = 1)
# ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
# ax.set_ylabel(r'$k_y \; \mathrm{(rad.m^{-1})}$', labelpad = 5)

x,y = np.indices(shift.shape)
x0 = shift.shape[0]//2
y0 = shift.shape[1]//2

# Compute radial distances from the center
r = np.sqrt((x - x0)**2 + (y - y0)**2)

# Flatten the matrix and distances for binning
r_flat = r.flatten()
values_flat = shift.flatten()

# Define bins based on the bin size
bin_size = 1
max_radius = np.max(r_flat)
bin_edges = np.arange(0, max_radius + bin_size, bin_size)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# draw circles 
idx_radius = 24
theta = np.linspace(0,2*np.pi,100)
radius_list = bin_edges[idx_radius: idx_radius + 2]

x_circle = np.array([x0 + radius*np.cos(theta) for radius in radius_list])
y_circle = np.array([y0 + radius*np.sin(theta) for radius in radius_list])

ax.plot(x_circle[0,:],y_circle[0,:],'r')
ax.plot(x_circle[1,:],y_circle[1,:],'r')









