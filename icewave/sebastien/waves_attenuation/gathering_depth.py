# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 11:39:01 2025

@author: sebas
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from datetime import datetime

import os 
import glob 
import pickle
import csv

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

dates = ['0221','0223','0226']

fig_folder = 'K:/sauvegarde_partielle_hublot_20250206/Share_hublot/Data/Results/Drones/thickness/'

#%% Load results for all selected days

data = {}

for date in dates :
    path2data = 'K:/sauvegarde_partielle_hublot_20250206/Share_hublot/Data/' + date + '/Drones/'
    # search all .mat file 
    filelist = glob.glob(path2data + '**/Drone_PIV_main_results*.mat',recursive = True)

    # build a dictionnary of results 
    data[date] = {}
    
    for i in range(len(filelist)):
    
        current_file = filelist[i]
        
        with h5py.File(current_file, 'r') as fmat:
            current_dict = {}
            
            print('Top-level keys : ', list(fmat.keys()))
        
            current_dict = mat2py.mat_to_dict(fmat['main_results'],fmat['main_results'])
            
        drone_ID = current_dict['drone_ID']
        if drone_ID == 'Bernache':
            drone_ID = 'bernache'
            current_dict['drone_ID'] = 'bernache'
            
        exp_ID = current_dict['exp_ID']
        
        if drone_ID not in data[date].keys(): # create key if drone_ID does not already exists
            data[date][drone_ID] = {}
            
        data[date][drone_ID][exp_ID] = current_dict # store current_dictionnary
            
#%% Load records for drones, extracted from .SRT files

records = {}
for date in dates:
    path2records = f'K:/sauvegarde_partielle_hublot_20250206/Share_hublot/Data/{date}/Summary/'
    
    file2load = f'{path2records}records_{date}.pkl'
    
    with open(file2load,'rb') as pf:
        current_records = pickle.load(pf)
    
    record_drones = {}
    for drone_ID in data[date].keys():
        record_drones[drone_ID] = {}
        for exp_ID in data[date][drone_ID]:
            if drone_ID == 'bernache':
                new_drone_ID = 'Bernache'
            else :
                new_drone_ID = 'mesange'
            
            record_drones[drone_ID][exp_ID] = current_records['drones'][new_drone_ID][exp_ID][0]
            
    records[date] = record_drones

#%% Create dictionnary with UTC t0, date, water depth 

summary = {}

for date in data.keys():
    summary[date] = {}
    for drone_ID in data[date].keys():
        summary[date][drone_ID] = {}
        for exp_ID in data[date][drone_ID].keys():
            
            S_disp = data[date][drone_ID][exp_ID]['disp_relation']
            hw = S_disp['param']['h_w']
            txt =f'{date}_{drone_ID}_{exp_ID}, hw = {hw} m'
            print(txt)

            time = records[date][drone_ID][exp_ID]['time'][0]
            date_drone = records[date][drone_ID][exp_ID]['date'][0]
            latitude = records[date][drone_ID][exp_ID]['latitude'][0]
            longitude = records[date][drone_ID][exp_ID]['longitude'][0]
            
            UTC_t0 = f'{date_drone} {time}'
            print(UTC_t0)
            
            summary[date][drone_ID][exp_ID] = {'UTC_t0':UTC_t0,'latitude':latitude,'longitude':longitude,'hw':hw}

file2save = f'{fig_folder}summary_thicknesses.pkl'

with open(file2save,'wb') as pf:
    pickle.dump(summary,pf)
    
#%% Write summary in a .csv file 

thickness_csv = f'{fig_folder}thicknes_summary.csv'

headers = ['key_exp','UTC_t0','latitude','longitude','hw']


with open(thickness_csv,'w',newline = '') as file:
    writer = csv.DictWriter(file, fieldnames=headers)
    writer.writeheader()
    for date in summary.keys():  
        for drone_ID in summary[date].keys():
            for exp_ID in summary[date][drone_ID].keys():
                # create a small dictionnary to be written
                
                small_dict = {}
                key_exp = f'{date}_{drone_ID}_{exp_ID}'
                small_dict['key_exp'] = key_exp
                for header in headers[1:]:
                    small_dict[header] = summary[date][drone_ID][exp_ID][header]
                    
                writer.writerow(small_dict)
                    
                 
                

                
#                 UTC_t0 = 
                
#                 file.write(f'{key_exp} : {UTC_t0')
#     for key, value in my_dict.items():
#         file.write(f"{key}: {value}\n")
            

        

# with open('output.csv', 'w', newline='') as file:
#     writer = csv.DictWriter(file, fieldnames=my_dict.keys())
    
#     writer.writeheader()       # Write the column headers (keys)
#     writer.writerow(my_dict)   # Write the values as a single row
# records['drones']['Bernache']['12-waves_010'][0]['longitude'][0]







