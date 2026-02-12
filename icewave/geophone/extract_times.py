# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 17:34:21 2025

@author: sebas

Extract initial times for geophone line deployment 

"""

#%% Import modules 
import numpy as np
from datetime import datetime
import mplcursors
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy import read
import os
import pickle
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import warnings
from matplotlib.patches import PathPatch
from scipy.interpolate import make_interp_spline
import warnings
from matplotlib.patches import PathPatch
import time
import csv 

import icewave.geophone.package_geophone as geopack


#%% Set parameters 
year = '2026'
date = '0129' #date format, 'mmdd'
acqu_numb = '0003' #acquisition number 

path2data = os.path.join('F:/Rimouski_2026/',date,'Geophones/')

# set path to geophone correspondence table
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
# channel = 0  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files

geophones_spacing = 6 # space between geophones, in meter 
signal_length = 1 # duration in seconds 
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}

#%% Set plotting parameters 

# Parameters for plots
font_size_medium = 24
font_size_small = 18
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 1000 # dpi to save images 

plt.rcParams.update({
    "text.usetex": True}) # use latex

fig_folder = path2data + acqu_numb + '/' + 'Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

###################################################################
#%% -------------------- Loading geophones Data -------------------
###################################################################
seismic_data_streams,datetime_values,fs = geopack.build_data_streams(path2data + acqu_numb +'/',geophones_table_path)

##############################################################################################
#%%----------------------- PLOTTING SELECTED CHANNEL FOR ALL GEOPHONES ------------------ 
##############################################################################################
 
channel = 1 #0 for E, 1 for N, 2 for Z. 

selected_indices = range(channel, len(seismic_data_streams),3)
fig, ax = plt.subplots(figsize = fig_size)
for k, idx in enumerate(selected_indices):
    current_stream = seismic_data_streams[idx]
    print(current_stream[0])
    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) ) + geophones_spacing*k + geophones_spacing, 
            label=f"Stream {k}")
    
    # ax.plot(datetime_values, 
    #         (current_stream[0].data-np.mean(current_stream[0].data))  + geophones_spacing*k + geophones_spacing, 
    #         label=f"Stream {k}")

#################################################################################
#%%----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
#################################################################################
# Select signal segments manually for each track and each source 
# Store selected times in a dictionnary and save it as a .pkl file 


# open .pkl file and load dictionnary 
filename = 't1_to_time_' + date + '_' + year + '.pkl'
base = path2data
file2save = base + filename
if os.path.isfile(file2save):
    print('Time dictionnary already exists')
    with open(file2save,'rb') as f:
        time_dict = pickle.load(f)
else:
    print('No time dictionnary saved yet')
    time_dict = {}

composante = 'N' #Z , E or N -> direction de la source

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:21:48.00")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:23:05.68")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:24:23.20")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:45:24.85")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:18:59.75")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2026-01-29T17:20:05.80")

# S107 , S108 , S109
# key = 'd' + date + 'a' + acqu_numb + 'tS' + '107' + composante 
# time_dict[key] = UTCDateTime("2025-02-03T16:01:15.70")
# key = 'd' + date + 'a' + acqu_numb + 'tS' + '108' + composante 
# time_dict[key] = UTCDateTime("2025-02-03T16:02:42.58")
# key = 'd' + date + 'a' + acqu_numb + 'tS' + '109' + composante 
# time_dict[key] = UTCDateTime("2025-02-03T16:03:54.90")


# Save t0 dictionnary in pickle file 

with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)

print('Time dictionnary saved')


#%% Load time dictionnary 

signal_length = 1 # duration in seconds

# load data of intial times 
composante = 'Z'
channel = 2
direction = 1 # 1 ou 2

flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave
 
# assign a string to S values depending on the direction
if direction == 1 :
    S1 = '101' 
    S2 = '102' 
    S3 = '103'
if direction == 2: 
    S1 = '104' 
    S2 = '105' 
    S3 = '106'

# time dictionnary to be loaded

pkl_path = f'{path2data}t1_to_time_{date}_{year}.pkl'
with open(pkl_path, 'rb') as f:
    loaded_data = pickle.load(f)


# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)
t1_values = [t1,t2,t3]

# Create a matrix to store the seismic data
selected_indices = np.arange(channel, len(seismic_data_streams),3)
ta = t1 + signal_length 
num_samples = int(signal_length * fs) # Number of points generated in time_vector 
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

num_traces = len(selected_indices)

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((num_traces, third_dim_len, num_samples)) # matrix dim : geo_indices x nb_sources x time

# Extract the relevant data for each trace and populate the 3D matrix
for i, stream_index in enumerate(selected_indices):
    stream = seismic_data_streams[stream_index]
    print(stream[0])
    for t1_index, t1 in enumerate(t1_values):
        for j, trace in enumerate(stream):
            start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
            # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index, :] = trace.data[start_sample:end_sample]

print('Matrix computed for channel : ' + composante)

# Plot the data
fig, ax = plt.subplots()
#print(seismic_matrix.shape[0]) # Seismic_matrix 1000 values [2], 3 sources [1], 16 geophones [0]
for k in range(seismic_matrix.shape[0]): # size of the array(nb of rows) = nb of geophones
    for t1_index in range(seismic_matrix.shape[1]): #size of the array(nb of collumns) = nb of sources
        ax.plot(time_vector, seismic_matrix[k, t1_index, :] / max(np.abs(seismic_matrix[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}") # : is each values (1000) for one geophones and 1 source

# Visualization of the plot
# Show time values for t1, t2, t3 ???
ax.set_title(f"Seismic Data - {date}-{year}-{acqu_numb}\n" 
             f"Channel: {channel}, Source: {composante}, Direction: {direction}", fontsize = font_size_medium)
ax.set_ylabel('Normalized Seismic Data', fontsize = font_size_small)
ax.tick_params(axis='x',labelsize = font_size_small)
ax.tick_params(axis='y',labelsize = font_size_small)
plt.show()














