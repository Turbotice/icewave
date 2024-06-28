# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:17:03 2024

@author: sebas
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate 
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import pickle
import csv 
from scipy import io
import h5py

#%% Load raw matlab Data 

# Load data 
path = 'Y:/Banquise/Aurore/e_5mm/Data/3.5Hz_15mm_e_5mm_sans_laser/'
filename = 'PIV_processed_i00_Dt1_b1_W32_xROI1_width1491_yROI1_height949post-processed.mat'

print('Loading matlab data..')
mat = io.loadmat(path+filename)

f = h5py.File('somefile.mat','r')
data = f.get('data/variable1')
data = np.array(data) # For converting to a NumPy array

#%%

u = mat['u']
v = mat['v']