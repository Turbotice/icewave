# -*- coding: utf-8 -*-
"""

This script can be decomposed into 2 main Parts : 
    Part 1: extraction of initial times for different sources performed on field 
    Part 2: 
"""
#%% Import modules 
import numpy as np
import math
from datetime import datetime
import mplcursors
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy import read
import shutil
import os
import re
import tkinter as tk
from tkinter import messagebox
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
from scipy.interpolate import interp1d
#%%
%matplotlib qt
#plt.rcParams['text.usetex'] = False

#%% Set parameters 
year = '2025'
date = '0220' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 

ordi = 'dell_vasco'

if ordi == 'babasse':
    path2data = os.path.join('E:/Data/',date,'Geophones/')
elif ordi == 'dell_vasco':
    path2data = f'S:/Data/{date}/Geophones/'

###################
path2data = f'D:/BicWin2025/Data/{date}/Geophones/'
###################

fig_folder = path2data + 'Figures/'  # folder where figures are saved 
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

if ordi=='babasse':    
    geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
elif ordi=='dell_vasco':
    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared

channel = 0  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files

geophones_spacing = 3 # space between geophones, in meters 
signal_length = 1 # duration in seconds 
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}
ch = channel_dic[channel]

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
    "text.usetex": False})

fig_folder = path2data + acqu_numb + '/' + 'Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%%  import sismo functions and load geophones data
# import sismo functions #################
curdir = os.getcwd()
os.chdir('Z:/vasco/Geophones/BicWin2025/python_functions/')
from sismo_analysis import *
os.chdir(curdir)
"""
##########################################
# load geophones data :
geophones_dict, seismic_data_streams, first_stream, first_trace, tstart, tend, fs, selected_indices, channel_dic, ch, time_array, time_vector, datetime_values, data_vector, start_time_utc = load_geophone_data(date,acqu_numb,path2data,geophones_table_path,channel)
#########################################
"""

#%% -------------------- Loading geophones Data -------------------
###################################################################


# Load geophones table into a dictionary
geophones_dict = {}
with open(geophones_table_path, 'r') as table_file:
    for line in table_file:
        if line.startswith("Num"):
            continue  # Skip header line
        num, geo_sn = line.strip().split()
        # Assuming the last five digits are relevant for comparison
        last_five_digits = geo_sn[-5:]
        geophones_dict[last_five_digits] = num


# Read MiniSEED file directly
seismic_data_streams, miniseed_files = read_data(path2data +'/' +acqu_numb)

# Iterate over streams and rename traces
for i, stream in enumerate(seismic_data_streams):
    seismic_data_streams[i] = rename_traces(stream, geophones_dict)


# Sort the seismic_data_streams based on the custom sorting function
seismic_data_streams = sorted(seismic_data_streams, key=sort_key)

# Find the largest start time (tstart) and the smallest end time (tend) among all streams
tstart = max([stream[0].stats.starttime for stream in seismic_data_streams])
tend = min([stream[-1].stats.endtime for stream in seismic_data_streams])

# Truncate all streams between tstart and tend
for i, stream in enumerate(seismic_data_streams):
    for trace in stream:
        if trace.stats.starttime < tstart:
            trace.trim(tstart, tend, pad=True, fill_value=0)
        elif trace.stats.endtime > tend:
            trace.trim(tstart, tend, pad=True, fill_value=0)


# Extract time_vector and data_vector for the first stream outside the loop
first_stream = seismic_data_streams[0]
first_trace = first_stream[0]
time_array = first_trace.times() # experimental time array 
time_vector = convert_to_utc_times(first_trace, first_trace.times())
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
data_vector = first_trace.data
start_time_utc = UTCDateTime(first_trace.stats.starttime) 
fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 
#%%----------------------- PLOTTING SELECTED CHANNEL FOR WHOLE RECORDING ------------------ 
##############################################################################################
 
# c'est avec ce plot qu'on va choisir les temps pour les sources

#channel = 1 #0 for E, 1 for N, 2 for Z. 

# selected indices of the corresponding channel
# seismic_data_G16 = seismic_data_streams[45:]
# current_stream  = seismic_data_G16[channel] 
# print(current_stream[0])
# ax.plot(datetime_values, 
#         current_stream[0].data / max(np.abs(current_stream[0].data) ), 
#         label=f"Stream {k}")
# selected_indices = range(channel, len(seismic_data_streams) - 3,3)

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
    

#fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
# Create an instance of ZoomHandler
zoom_handler = ZoomHandler(ax, time_vector, data_vector)
fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder


#%%----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
#################################################################################
# Select signal segments manually for each track and each source 
# Store selected times in a dictionnary and save it as a .pkl file 


# open .pkl file and load dictionnary 
filename = 't1_to_time_flighttimeanalysis_' + date + '_' + year + '.pkl'
base = path2data
file2save = base + filename
if os.path.isfile(file2save):
    print('Time dictionnary already exists')
    with open(file2save,'rb') as f:
        time_dict = pickle.load(f)
else:
    print('No time dictionnary saved yet')
    time_dict = {}

"""
composante = 'Z'
# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:09:33.30")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:10:21.80")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:10:56.90")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:12:57.70")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:14:15.20")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:15:17.80")
"""


composante = 'E' #Z , E or N -> direction de la source

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-20T20:58:19.00")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-20T20:58:28.70")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-20T20:59:26.30")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-20T21:00:51.50")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-20T21:01:02.40")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-10T21:01:10.30")

"""
composante = 'N' #Z , E or N -> direction de la source

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-10T21:10:01.80")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-10T21:10:44.60")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-10T21:11:33.30")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:13:32.10")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:14:49.40")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-04T21:15:47.90")
"""

# Save t0 dictionnary in pickle file 

with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)

print('Time dictionnary saved')


#%%

flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave
direction = 1  
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
if ordi=='babasse': 
    base = f'E:/Data/{date}/Geophones/'
elif ordi=='dell_vasco':
    base = f'S:/Data/{date}/Geophones/'

###############################
base = f'D:/BicWin2025/Data/{date}/Geophones/'

pkl_path = base + 't1_to_time_flighttimeanalysis_' + date + '_' + year  + '.pkl'

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
num_samples = int(signal_length * first_trace.stats.sampling_rate) # Number of points generated in time_vector 
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

num_traces = len(selected_indices)

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((num_traces, third_dim_len, num_samples)) # matrix dim : geo_indices x nb_sources x time

# add G16 in first place 
# idx_G16 = selected_indices[-1]
# stream = seismic_data_streams[idx_G16]
# print(stream[0])
# for t1_index, t1 in enumerate(t1_values):
#     for j, trace in enumerate(stream):
#         start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
#         # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
#         end_sample = start_sample + num_samples
#         seismic_matrix[0, t1_index, :] = trace.data[start_sample:end_sample]


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
plt.grid(True)

# %% partie où on correle les signaux

# seismic matrix has shape : (num_traces , number of t0 times , signal length)



geophones_spacing = 1 # en metres

t1_index = 0


wave_speed_values = []

for i in range(1,num_traces):

        geophone_pair = (0,i)


        time_shift = detect_time_shift(seismic_matrix[geophone_pair[0],t1_index,:],seismic_matrix[geophone_pair[1],t1_index,:],freq_acq=fs,plot=True)
        # for non dispersive modes :
        wave_speed = geophones_spacing * (geophone_pair[1] - geophone_pair[0])/time_shift
        print('wave velocity (assuming non dispersive mode) :', wave_speed)
        wave_speed_values.append(wave_speed)

wave_speed_values = np.array(wave_speed_values)
# %%
