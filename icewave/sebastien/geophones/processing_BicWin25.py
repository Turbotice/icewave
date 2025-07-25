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

#%% Set parameters 
year = '2025'
date = '0204' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 

path2data = os.path.join('U:/Data/',date,'Geophones/')

# set path to geophone correspondence table
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
# channel = 0  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files

geophones_spacing = 3 # space between geophones, in meter 
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

#----------------------------------------------------------------------------------------------
def read_data(path2data):
    """ Read path to folder where data are stored, 
    Extract each stream and create a list of streams 
    Inputs :
        - path2data : path to data files
    Outputs : 
        - streams : list of stream of all geophones 
        - minseed_files : list of files .miniseed """
    miniseed_files = []

    # Iterate over files in the specified directory
    for filename in os.listdir(path2data):
        if filename.endswith(".miniseed"):
            file_path = os.path.join(path2data, filename)
            miniseed_files.append(file_path)

    # Read MiniSEED files
    streams = []
    for file_path in miniseed_files:
        stream = read(file_path)
        streams.append(stream)

    return streams, miniseed_files
#----------------------------------------------------------------------------------------------
def convert_to_utc_times(trace, time_vector):
    """ Create a UTC-time time vector for a given trace 
    Inputs :
        - trace : trace object (obspy) 
        - time_vector : array of time elapsed since beginning of the acquisition 
    Output : 
        - array of UTC time (datetime type)"""
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]
#----------------------------------------------------------------------------------------------
def rename_traces(stream, geophones_dict):
    """ Rename geophone traces using geophone_table 
    Inputs : 
        - stream : stream containing traces of each geophone 
        - geophones_dict : dictionnary of geophone correspondance 
    Output : 
        - sorted_stream : stream sorted by name using geophone table of correspondance """
    for trace in stream:
        # Extract the last 5 digits after "SS."
        last_five_digits = trace.stats.station[-5:]

        # Check if the last five digits exist in the geophones dictionary
        if last_five_digits in geophones_dict:
            # Get the corresponding 2-digit number from the geophones dictionary
            new_station_code = geophones_dict[last_five_digits]

            # Replace the last 5 digits in the existing station code with the new ones
            trace.stats.station = f"{trace.stats.station[:-5]}{new_station_code}.{trace.stats.channel}"
        else:
            print(
                f"Warning: No entry found for {last_five_digits} in the geophones table.")

    # Sort traces based on the new station codes
    sorted_stream = sorted(stream, key=lambda trace: trace.stats.station)

    return sorted_stream

#---------------------------------------------------------------------------------------------
def sort_key(trace):
    """ Create a sorting key based on stream stations """
    return trace[0].stats.station

#----------------------------------------------------------------------------------------------

class ZoomHandler:
    """ Create a class to perform zoom on matplotlib graphs """
    def __init__(self, ax, time_vector, data_vector):
        self.ax = ax
        self.time_vector = time_vector
        self.data_vector = data_vector
        self.original_xlim = ax.get_xlim()
        self.original_ylim = ax.get_ylim()

        # Initialize rectangle selector
        self.rs = RectangleSelector(ax, self.on_rectangle_select, drawtype='box', useblit=True, button=[1],
                                    minspanx=5, minspany=5, spancoords='pixels', interactive=True)

    def on_click(self, event):
        if event.dblclick:
            self.reset_zoom()

    def on_rectangle_select(self, eclick, erelease):
        # Extract rectangle coordinates
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Apply zoom to the selected area
        self.ax.set_xlim(min(x1, x2), max(x1, x2))
        self.ax.set_ylim(min(y1, y2), max(y1, y2))
        plt.draw()

    def reset_zoom(self):
        self.ax.set_xlim(self.original_xlim)
        self.ax.set_ylim(self.original_ylim)
        plt.draw()

def disconnect_toolbar_events(fig):
    toolbar = fig.canvas.toolbar
    toolbar.toolmanager.remove_tool('zoom')
    toolbar.toolmanager.remove_tool('pan')

# Function to reconnect the toolbar events
def reconnect_toolbar_events(fig):
    toolbar = fig.canvas.toolbar
    toolbar.toolmanager.add_tool('zoom', plt._tools.Zoom)
    toolbar.toolmanager.add_tool('pan', plt._tools.Pan)

###################################################################
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

# unsorted_seismic_data = seismic_data_streams.copy()
# streams_G11 = unsorted_seismic_data[33:36]
# streams_G12 = unsorted_seismic_data[30:33]
# unsorted_seismic_data[30:33] = streams_G11
# unsorted_seismic_data[33:36] = streams_G12

# seismic_data_streams = unsorted_seismic_data

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

##############################################################################################
#%%----------------------- PLOTTING SELECTED CHANNEL FOR ALL GEOPHONES ------------------ 
##############################################################################################
 
channel = 2 #0 for E, 1 for N, 2 for Z. 

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
time_dict[key] = UTCDateTime("2025-02-12T21:08:48.20")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-12T21:09:53.50")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-12T21:10:55.80")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-12T21:13:12.50")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-12T21:14:08.40")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-12T21:15:08.85")

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















