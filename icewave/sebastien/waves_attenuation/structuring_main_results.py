# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 13:53:03 2025

@author: sebas
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os 
import glob 
import pickle
import h5py
import icewave.tools.matlab2python as mat2py

# import seb plotting package 
import icewave.sebastien.set_graphs as set_graphs

#%% Set plotting parameters 
 
font_size_medium = 30
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 100 # dpi to save images 
marker_size_plot = 12

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%%

dates = ['0223','0226']

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Attenuation_BicWin2024/'


#%% Load results for a given day 
data = {}

date = '0223'
path2data = 'K:/Share_hublot/Data/' + date + '/Drones/'
# search all .pfl file 
filelist = glob.glob(path2data + '**/attenuation*.pkl',recursive = True)

# build a dictionnary of results 
data[date] = {}

#%%
i = 0
current_file = filelist[i]

drone_ID = current_file[-25:-17]
if drone_ID != 'bernache':
    drone_ID = 'mesange'
print(drone_ID)
exp_ID = current_file[-16:-4]
print(exp_ID)
# with open(current_file,'rb') as pf:
#     current_dict = pickle.load(pf)

# print('Top-level keys : ', list(current_dict.keys()))

#%%
data = {}

for date in dates :
    path2data = 'K:/Share_hublot/Data/' + date + '/Drones/'
    # search all .mat file 
    filelist = glob.glob(path2data + '**/attenuation*.pkl',recursive = True)

    # build a dictionnary of results 
    data[date] = {}
    
    for i in range(len(filelist)):
    
        current_file = filelist[i]
        
        with open(current_file, 'rb') as pf:
            current_dict = pickle.load(pf)
            
            print('Top-level keys : ', list(current_dict.keys()))
            
        drone_ID = current_dict['drone_ID']
        if drone_ID == 'Bernache':
            drone_ID = 'bernache'
            current_dict['drone_ID'] = 'bernache'
            
        exp_ID = current_dict['exp_ID']
        
        if drone_ID not in data[date].keys(): # create key if drone_ID does not already exists
            data[date][drone_ID] = {}
            
        data[date][drone_ID][exp_ID] = current_dict # store current_dictionnary































