# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 16:21:09 2025

@author: sebas

This script aims at getting the orientation of the swell on the 11/02/2025 compared to the deployed optical fiber
To do so we use the orientation of the swell to phones, using a triangle of phone deployed in front of the fiber

Using the GPS coordinates of phones and GPS coordinates of fiber, and swell orientation to the phones, 
we deduce the orientation of swell to the fiber axis 

"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
from scipy.io import loadmat

import icewave.tools.matlab2python as mat2py
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%%








