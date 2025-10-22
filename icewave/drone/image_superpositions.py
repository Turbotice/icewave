# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 16:22:25 2025

@author: sebas
"""

import os
import re
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
import scipy
import csv

import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.field.gps as field_gps

import icewave.field.drone as field_drone
import icewave.gps.gps as gps_pack

#%% Get filelist of images 

date = '0210'
drone_ID = 'bernache'
exp_ID = '09-ortho_004'

base = 'U:/Data/'

records = field_drone.get_records('0210',jpg = False)
 # = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'
