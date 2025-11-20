# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 16:51:11 2025

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
import imageio as iio
import cv2 as cv
import h5py
import csv

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Load (f,k)

date = '0211'
base = 'U:/Data/'
save_file = f'{base}{date}/DAS/FK_0211_active.pkl'

with open(save_file, 'rb') as fid:
   data = pickle.load(fid)
   
FK = data['FK']
k = data['k']
freq = data['freq']
plt.figure(figsize=(10, 6))
plt.imshow(np.transpose(FK), aspect='auto', extent=[k[0], k[-1], freq[-1], freq[0]], cmap='gray_r')
plt.colorbar(label='Projection Magnitude')
plt.xlabel('Wavenumber (rad/m)')
plt.ylabel('Frequency (Hz)')
plt.title('SVD Projection')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()


#%% Peak extraction 

