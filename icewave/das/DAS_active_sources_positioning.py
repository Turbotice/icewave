# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 15:49:58 2025

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

# module for map plotting
import cmocean.cm as cmo
import alphashape
from shapely.geometry import MultiPoint
from shapely.geometry import Point, Polygon
from scipy.interpolate import griddata
import cartopy.io.shapereader as shpreader

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader

# icewave modules
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps
import icewave.das.DAS_package as DS


import icewave.analysis.bathy as bathy
import stephane.display.graphes as graphes
import icewave.gps.gps as gps
import icewave.geometry.tables as tables

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Load sources GPS position and time from Stephane's GPS
 
data_gps = {}
dates = ['0211','0212']
for date in dates:
    data_gps[date] = {}
    
    path2data = f'U:/Data/{date}/GPS/'
    gpx_filelist = glob.glob(path2data + '*.gpx')
    
    gpx_file = gpx_filelist[0]
    # get gpx
    gpx = gps.get_wpts(gpx_file)
    
    # get waypoints position
    # use Map_table to create a dictionnary
    data_gps[date] = tables.dict_from_gpx(gpx,path2data)
    
# keep only waypoints related to sources and borne
sources_gps = {}
for date in data_gps.keys():
    sources_gps[date] = {}
    for key in data_gps[date].keys():
        instrument = key[0]
        
        if instrument == 'S':
            sources_gps[date][key] = data_gps[date][key]
            
        elif key == 'borne_01':
            sources_gps[date][key] = data_gps[date][key]
            
#%% Plot sources each day


#%% Ludovic's sources detection 

sources_xt = {}
sources_xt['0212'] = [
    {'acq_num': 0, 'x0': 58, 'tmin': 164.30},
    {'acq_num': 0, 'x0': 85, 'tmin': 307.29},
    {'acq_num': 1, 'x0': 130, 'tmin': 53.68},
    {'acq_num': 1, 'x0': 152, 'tmin': 135.12},
    {'acq_num': 1, 'x0': 170, 'tmin': 235.50},
    {'acq_num': 1, 'x0': 190, 'tmin': 327.4},
    {'acq_num': 1, 'x0': 211, 'tmin': 404.57},
    {'acq_num': 1, 'x0': 231, 'tmin': 492.65},
    {'acq_num': 2, 'x0': 260, 'tmin': 129.66},
    {'acq_num': 2, 'x0': 290, 'tmin': 308.40},#298.80 #332.70
    {'acq_num': 2, 'x0': 327, 'tmin': 414.66},
    {'acq_num': 2, 'x0': 348, 'tmin': 495.54},
    {'acq_num': 3, 'x0': 387, 'tmin': 71.68},#89.34
    {'acq_num': 3, 'x0': 410, 'tmin': 146.79},
    {'acq_num': 3, 'x0': 426, 'tmin': 230.82},
    {'acq_num': 3, 'x0': 442, 'tmin': 310.76},
    {'acq_num': 4, 'x0': 460, 'tmin': 185.17},
    {'acq_num': 4, 'x0': 490, 'tmin': 321.03},
    {'acq_num': 4, 'x0': 505, 'tmin': 563.85},
    {'acq_num': 5, 'x0': 525, 'tmin': 46.18},
    {'acq_num': 5, 'x0': 545, 'tmin': 122.85},
    {'acq_num': 5, 'x0': 565, 'tmin': 200.61},
    {'acq_num': 5, 'x0': 585, 'tmin': 322.68},
]

sources_xt['0211'] = [
     {'acq_num': 0, 'x0': 40, 'tmin': 125.7},
     {'acq_num': 0, 'x0': 60, 'tmin': 356.1},
     {'acq_num': 0, 'x0': 90, 'tmin': 494.4},
     {'acq_num': 0, 'x0': 120, 'tmin': 560.5},
     {'acq_num': 1, 'x0': 150, 'tmin': 245.1},
     {'acq_num': 1, 'x0': 170, 'tmin': 329.6},
     {'acq_num': 2, 'x0': 190, 'tmin': 9.2},        #  NOT GOOD
     {'acq_num': 2, 'x0': 290, 'tmin': 453.4},      #  NOT GOOD
     {'acq_num': 2, 'x0': 320, 'tmin': 573.7},      #  NOT GOOD
     {'acq_num': 3, 'x0': 330, 'tmin': 16.5},       #  NOT GOOD
     {'acq_num': 3, 'x0': 365, 'tmin': 189.7}, 
     {'acq_num': 3, 'x0': 385, 'tmin': 303.7}, 
     {'acq_num': 3, 'x0': 405, 'tmin': 382.4}, 
     {'acq_num': 3, 'x0': 420, 'tmin': 468.9}, 
     {'acq_num': 3, 'x0': 445, 'tmin': 549}, 
     {'acq_num': 4, 'x0': 485, 'tmin': 106}, 
     {'acq_num': 4, 'x0': 505, 'tmin': 183}, 
     {'acq_num': 4, 'x0': 525, 'tmin': 265.6},]

#%% Compute sources_xt UTC time from UTC file time 

# collect files for a given date
date = '0211'
base = 'U:/Data/'
path2DAS_data = f'{base}{date}/DAS/*_UTC.h5'
filelist = glob.glob(path2DAS_data)


# get UTC initial time for each acquisition 









# Load parameters for DAS
# path2DAS_param = f'{base}/parameters_Febus_2025.pkl'
# date = '0211'
# fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# file2load = filelist[0]
# Nb_minutes = 1
# stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
# format_date = '%Y-%m-%dT%H-%M-%S'
# label_UTC0 = UTC_stack[0,0].strftime(format_date)

