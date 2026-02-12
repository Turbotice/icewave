# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:12:01 2025

@author: sebas

Script to study geophones data from 0914 during Refuge Arctic Expedition on 14/09/2024
"""

#%% Import modules 
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy import read
import os
import pickle
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
from scipy.interpolate import make_interp_spline
from scipy.signal import ShortTimeFFT, windows, iirnotch, filtfilt
import scipy.signal 
import warnings
from matplotlib.patches import PathPatch
import time
import csv 

import icewave.tools.time_package as time_package
import icewave.sebastien.geophones.package_geophone as geophone
import icewave.tools.Fourier_tools as FT
import icewave.tools.matlab_colormaps as matcmaps

parula_map = matcmaps.parula()

#%% Set parameters 

year = '2024'
date = '0914' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 

# path2data = os.path.join('U:/Data/',date,'Geophones/')
path2data = f'F:/Amundsen_RA_2024/Data/{year}/{date}/Geophones/'

# set path to geophone correspondence table
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
# channel = 0  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files

geophones_spacing = 8 # space between geophones, in meter 
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
    
#%% FUNCTION SECTION 


###################################################################
#%% -------------------- Loading geophones Data -------------------
###################################################################

path2miniseeds = f'{path2data}{acqu_numb}/'
seismic_data_streams, datetime_values, fs = geophone.build_data_streams(path2miniseeds, geophones_table_path)


##############################################################################################
#%%----------------------- PLOTTING SELECTED CHANNEL FOR ALL GEOPHONES ------------------ 
##############################################################################################
 
channel = 0 #0 for E, 1 for N, 2 for Z. 

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
    

################################################################################
#%% ----------------- LOAD TIME DICTIONNARY ------------------------------------
################################################################################

path2times = f'{path2data}t1_to_time_{date}_{year}_first_hit.pkl'
with open(path2times, 'rb') as f:
    time_dict = pickle.load(f)



#%% Extract signal associated to a single source

signal_length = 10 #duration in seconds

# choose source
source_idx = '008'
channel = 1 # channel on which we look at recorded signals
composante = channel_dic[channel] # type of excitation performed 

key_source = f'd{date}a{acqu_numb}tS{source_idx}{composante}'
t0 = time_dict.get(key_source)

# select channel for all geophones 
selected_indices = np.arange(channel,len(seismic_data_streams),3)
tend = t0 + signal_length # final time of time segment

# initialize a 2D matrix 
num_samples = int(signal_length * fs) # Number of points generated in time_vector 
time_vector = np.linspace(t0.timestamp, tend.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

# Dynamically determine the length of the third dimension
seismic_matrix = np.zeros((len(selected_indices), num_samples)) # matrix dim : geo_indices x time

# Extract the relevant data for each trace and populate a 3D matrix
for i, stream_index in enumerate(selected_indices):
    stream = seismic_data_streams[stream_index]
    print(stream[0])
    for j, trace in enumerate(stream):
        start_sample = int((t0 - trace.stats.starttime) * fs)
        # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
        end_sample = start_sample + num_samples
        seismic_matrix[i, :] = trace.data[start_sample:end_sample]


# Plot seismic 2D matrix 

# fig, ax = plt.subplots(figsize = fig_size)
# for k in range(seismic_matrix.shape[0]): 
#     ax.plot(time_vector, seismic_matrix[k, :] / max(np.abs(seismic_matrix[k, :])) + geophones_spacing*k)
    
# title = f'Source {source_idx} channel {composante} - Raw data'
# ax.set_title(title, fontsize = font_size_small)
# ax.tick_params(axis='x',labelsize = font_size_small)
# ax.tick_params(axis='y',labelsize = font_size_small)
# plt.show()

#############################################################################################
#%% ---------------- Filter out frequencies related to Amundsen  ----------------------------
#############################################################################################

# Load frequencies to cut from the Amundsen 
path2frequencies = f'{path2data}Frequencies_Amundsen_{year}_{date}_channel_{channel_dic[channel]}.pkl'
print(path2frequencies)
with open(path2frequencies,'rb') as pf:
    fc = pickle.load(pf)


fs = 1000
Q = 30.0  # Quality factor
filtered_matrix = seismic_matrix
for f0 in fc:

    # Design notch filter
    b, a = iirnotch(f0, Q, fs)
    filtered_matrix = filtfilt(b,a,filtered_matrix)

fig, ax = plt.subplots(figsize = fig_size)
for k in range(seismic_matrix.shape[0]): 
    ax.plot(time_vector, filtered_matrix[k, :] / max(np.abs(filtered_matrix[k, :])) + geophones_spacing*k)
    
title = f'Source {source_idx} channel {composante} - Filtered data'
ax.set_title(title, fontsize = font_size_small)
ax.tick_params(axis='x',labelsize = font_size_small)
ax.tick_params(axis='y',labelsize = font_size_small)
plt.show()

##############################################################################################
#%% ----------------------- Extract Young modulus and Poisson coefficient --------------------
##############################################################################################

# select geophone indices (#1  to #16)
geoph_selection = [13,14,15,16]
sub_matrix = np.zeros((len(geoph_selection),filtered_matrix.shape[1]))

for i,idx_geoph in enumerate(geoph_selection):
    sub_matrix[i,:] = filtered_matrix[idx_geoph - 1,:]
    
fig, ax = plt.subplots(figsize = fig_size)
for k in range(sub_matrix.shape[0]): 
    ax.plot(sub_matrix[k, :] / max(np.abs(sub_matrix[k, :])) + geophones_spacing*k)
    
title = f'Source {source_idx} channel {composante} - Selected data'
ax.set_title(title, fontsize = font_size_small)
ax.tick_params(axis='x',labelsize = font_size_small)
ax.tick_params(axis='y',labelsize = font_size_small)
plt.show()

#%% Compute correlations between first signal of sub_matrix 
# select indices 
idx_start = 50
idx_end = 200

# compute correlation 
correlation = np.zeros((sub_matrix.shape[0],int(2*(idx_end - idx_start) - 1)))
for i in range(4):
    correlation[i,:] = scipy.signal.correlate(sub_matrix[0,idx_start:idx_end],sub_matrix[i,idx_start:idx_end])

lag = scipy.signal.correlation_lags(len(sub_matrix[0,idx_start:idx_end]),len(sub_matrix[0,idx_start:idx_end]))
tlag = lag/fs
fig, ax = plt.subplots()
for i in range(correlation.shape[0]):
    ax.plot(tlag,correlation[i,:])

# search for maximum of each correlation 
max_corr = np.argmax(abs(correlation),axis = 1)
delta_t = abs(tlag[max_corr])

x = np.arange(0,4)*geophones_spacing

p = np.polyfit(delta_t,x,1)
print(p)
dt_th = np.linspace(delta_t.min(),delta_t.max())

fig, ax = plt.subplots()
ax.plot(delta_t,x,'o')
ax.plot(dt_th,np.polyval(p,dt_th),'r-')


#%%

correlation = scipy.signal.correlate(sub_matrix[0,idx_start:idx_end], sub_matrix[0,idx_start:idx_end])
fig, ax = plt.subplots()
ax.plot(correlation)

#%% Compute FFT of the sub_matrix

facq = [1/geophones_spacing , fs]
shift, k, omega = FT.fft_2D(sub_matrix[:,5500:5800],facq,add_pow2 = [2,0])

fig, ax = plt.subplots()
ax.imshow(abs(shift), cmap = parula_map)






###############################################################################################
#%% ------------------------------ DRAFT ------------------------------------------------------
###############################################################################################





