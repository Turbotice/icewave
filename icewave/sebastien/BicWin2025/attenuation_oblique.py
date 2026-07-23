# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:44:09 2026

@author: sebas
"""

#%% 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean

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

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Load data

base = 'U:/Data/'
date = '0211'
drone_ID = 'mesange'
exp_ID = '05-waves_001'
suffixe = f'{date}_{drone_ID}_{exp_ID}'

fig_folder = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load matfile 

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

#%% 

print(S['DRONE']['alpha_0']*180/np.pi)

#%% Supress quadratic noise

Vx = FT.supress_quadratic_noise(np.transpose(S['Vx'],(1,0,2)),S['x'],S['y'])
Vy = FT.supress_quadratic_noise(np.transpose(S['Vy'],(1,0,2)),S['x'],S['y'])
Vx = np.transpose(Vx,(1,0,2))
Vy = np.transpose(Vy,(1,0,2))

#%% Show apparent velocity fields

extents_pix = np.array([S['PIXEL']['x_pix'].min(),S['PIXEL']['x_pix'].max(),
                    S['PIXEL']['y_pix'].min(),S['PIXEL']['y_pix'].max()])

frame = 3600

fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,8))
imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.imshow(Vx[:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_x$')
    else:
        imsh.append(ax.imshow(Vy[:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_y$')
   
    ax.set_xlabel(r'$x_{pix}$')     

axs[0].set_ylabel(r'$y_{pix}$')

plt.tight_layout()

#%% Show georectified apparent velocity fields

cmap = cmocean.cm.balance
frame = 3000

X_bounds = np.array([-50,50])
Y_bounds = np.array([-40,40])
Vx_values = np.array([-1,1])
Vy_values = np.array([-1,2])

fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,8))
imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.pcolormesh(S['X'],S['Y'],Vx[:,:,frame],shading = 'gouraud',cmap = cmap,
                                 vmin = Vx_values[0],vmax = Vx_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_X$')
    else:
        imsh.append(ax.pcolormesh(S['X'],S['Y'],Vy[:,:,frame],shading = 'gouraud',cmap = cmap,
                                 vmin = Vy_values[0],vmax = Vy_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_Y$')

    ax.set_aspect(1)

#%% Check values of coefficients

alpha = S['DRONE']['alpha_0']
h = S['DRONE']['h_drone']
f = S['DRONE']['focale']
Y = S['Y']
X = S['X']

z_star = dp.get_zstar(h,alpha,Y)

coeff = (np.cos(alpha) + (np.sin(alpha)**2)*Y/z_star)*f/z_star

fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],coeff,shading = 'gouraud',cmap = parula_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
ax.set_aspect(1)

#%% Try inversion of uz

frame = 1
cmap = cmocean.cm.balance
uz = Vy/coeff[:,:,None]

uz_values = np.array([-0.08,0.08])
fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],uz[:,:,frame],shading = 'gouraud',cmap = cmap,
                    vmin = uz_values[0],vmax = uz_values[1])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_z$')
ax.set_aspect(1)

#%% Deduce ux
cmap = cmocean.cm.balance
ux = Vx*z_star[:,:,None]/f - uz*np.sin(alpha)*X[:,:,None]/z_star[:,:,None] 

fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],ux[:,:,frame],shading = 'gouraud',cmap = cmap)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_x$')
ax.set_aspect(1)

#%% Define a regular grid

#%% Define grid 

minx = -100
maxx = 100
miny = S['Y'].min()
maxy = S['Y'].max()

artificial_facq_x = 1/0.9 # artificial spatial frequency in box/meter
grid_x, grid_y = np.meshgrid(
    np.linspace(minx, maxx, int((maxx - minx)*artificial_facq_x)),  
    np.linspace(miny, maxy, int((maxy - miny)*artificial_facq_x))
)


#%% Compare spatial sampling 
diff_y = np.diff(S['Y'][0,:])
diff_x = np.mean(np.diff(S['X'],axis = 1),axis = 1)

fig, ax = plt.subplots()
ax.plot(abs(diff_y),'.',label = '$\Delta_y$')
ax.plot(abs(diff_x),'.', label = '$\Delta_x$')
ax.set_ylim([0,2])
ax.legend()

#%% Interpolate uz over the regular grid

points = np.array([S['X'].ravel(),S['Y'].ravel()]).T
interp_uz = np.zeros((grid_x.shape[0],grid_x.shape[1],uz.shape[2]))
for frame in range(uz.shape[2]):
    print(frame)
    interp_uz[:,:,frame] = dp.interpolate_field(points,uz[:,:,frame],grid_x,grid_y) 


#%% Show interpolated uz field 

frame = 3000
fig, ax = plt.subplots()
imsh = ax.pcolormesh(grid_x.T,grid_y.T,interp_uz[:,:,frame].T,shading = 'gouraud',cmap = parula_map,
                    vmin = -0.04, vmax = 0.08)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_y$')
ax.set_aspect(1)


#%% Save interpolated uz field 

data = {}
data['DRONE'] = S['DRONE']
data['GPS'] = S['GPS']
data['uz'] = uz
data['interp_uz'] = interp_uz
data['grid_x'] = grid_x
data['grid_y'] = grid_y
data['SCALE'] = S['SCALE']
data['artificial_facq_x'] = artificial_facq_x
data['t'] = S['t']

file2save = f'{path2data}interpolated_uz_quadratic_correction{date}_{drone_ID}_{exp_ID}.h5'
rw.save_dict_to_h5(data, file2save)