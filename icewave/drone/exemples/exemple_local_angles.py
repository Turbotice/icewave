# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:12:28 2025

@author: sebas
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize

import cv2 as cv
import glob
import os
import pickle 
import h5py
import re

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.drone.drone_projection as dp 

# import seb plotting package 
import icewave.sebastien.set_graphs as set_graphs
parula_map = matcmaps.parula()

#%% Load data from a PIV 

base = 'F:/Rimouski_2024/Data/'

date = '0226'
drone_ID = 'mesange'
exp_ID = '10-waves_005'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = path2data + 'PIV_processed_i00_Dt5_b1_W32_xROI600_width3240_yROI1_height2159_scaled.mat'

with h5py.File(file2load, 'r') as fmat:
    data_wave = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    S = mat2py.transpose_PIVmat_fields(S)

#%% Get a map of angles 

# S['PIXEL'] contains pixel position of center of each PIV box
# create a meshgrid

set_graphs.set_matplotlib_param('single')
xpix,ypix = np.meshgrid(S['PIXEL']['x_pix'],S['PIXEL']['y_pix'],indexing = 'xy')

print(xpix.shape)
fig, ax = plt.subplots()
imsh = ax.imshow(xpix,cmap = parula_map)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

#%%
 
beta_x,beta_y,beta_r = dp.get_angles(xpix,ypix, S['PIXEL']['x0'], S['PIXEL']['y0'], 
                                     S['DRONE']['focale'])

fig, ax = plt.subplots()
imsh = ax.imshow(beta_r,cmap = parula_map)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)








