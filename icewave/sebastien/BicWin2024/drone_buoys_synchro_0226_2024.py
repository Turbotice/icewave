# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 16:48:15 2025

@author: sebas

Synchronization of buoys and drones for 0226/2024

"""

import numpy as np
import os
import glob
import pandas as pd 
import h5py
import pickle
import cv2 as cv

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.signal as signal
from scipy.interpolate import LinearNDInterpolator, griddata, RegularGridInterpolator
from datetime import datetime,timedelta
import pytz

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.rw_data as rw
import icewave.gps.gps_seb as gps_seb
import icewave.gps.gps as gps
import icewave.geometry.tables as tables

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Load data

base = 'K:/Share_hublot/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '12-FRAC_001'

path2data = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'

fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

# load data 
file2load = f'{path2data}data_buoys_velocity_{date}_{drone_ID}_{exp_ID}.h5'
data = rw.load_dict_from_h5(file2load)


