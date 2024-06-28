#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:38:04 2024

@author: moreaul
"""
%matplotlib qt

import os
import shutil
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime, date, time 
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate 

##############################################################
#%% --------- SETTING PARAMETERS AND PATHS -------------------
##############################################################

path2data = 'E:/Rimouski_2024/Data/2024/0211/Geophones/'
geophones_table_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Geophones/geophones_table'

#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0001' 
channel = 2 #0 for E, 1 for N, 2 for Z. 
channel_correspondance = ['E','N','Z']
gain = 12 # gain in dB 
scale_factor = 10**(-gain/10)
geophones_spacing = 3 # space between geophones, in meter 

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
img_quality = 1200 # dpi to save images 


fig_folder = path2data + 'Couder_conference/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
######################################################
#%% ----------------------- FUNCTIONS ------------------ 
######################################################

#warnings.filterwarnings("ignore", category=np.ComplexWarning)
def fn_svd(signals, fs, xs, name, issaving, *varargin):
    threshold_user = None
    rang = None

    if varargin:
        if varargin[0] == 'threshold':
            threshold_user = varargin[1]
        elif varargin[0] == 'rang':
            rang = varargin[1]
        else:
            print('varargin(1) unknown')
            return
    else:
        print('varargin empty')

    Nreceiv, Nt, Nemit = signals.shape

    # time domain fft
    Nf = 1024
    f = (np.arange(Nf//2 + 1) / Nf) * (fs if fs else 1)
    f_axename = 'f/fs' if not fs else 'f'
    SIGNALS = fft(signals, Nf, axis=1)
    SIGNALS = SIGNALS[:, :Nf//2 + 1, :]

    # svd
    U = np.zeros((Nreceiv, Nf//2 + 1, Nemit))
    S = np.zeros((Nemit, Nf//2 + 1, Nemit))
    V = np.zeros((Nemit, Nf//2 + 1, Nemit))
    D = np.zeros((Nemit, Nf//2 + 1))

    for ii in range(Nf//2 + 1):
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False)
        D[:, ii] = np.diag(S[:, ii, :])

    for ne in range(Nemit):
        titi = 20 * np.log10(D[ne, :] / np.max(D[0, :]))
        plt.plot(f, titi, label=f'Slice {ne}')

    if threshold_user is None:
        plt.xlabel('frequency (Hz)')
        plt.ylabel('Singular values (in dB of peak value)')
        [fcut, sigmacutsup] = plt.ginput(5)
        fcut = [f[0]] + fcut.tolist() + [f[-1]]
        sigmacutsup = [sigmacutsup[0]] + sigmacutsup.tolist() + [sigmacutsup[-1]]
        sigmacutsup = np.interp(f, fcut, sigmacutsup)
    else:
        sigmacutsup = np.full(int(Nf/2+1), threshold_user)
        sigmacutsup = threshold_user
        #print(sigmacutsup)              
    #plt.plot(f, sigmacutsup)
    
    
    # uncomment #U[:, idx, ne] = 0 for applying the threshold, otherwise
    for ne in range(Nemit):
        titi = 20 * np.log10(D[ne, :] / np.max(D[0, :]))
        idx = np.where(titi <= sigmacutsup)[0]
        U[:, idx, ne] = 0

    if issaving:
        plt.savefig(name + '_sv')

    # projection onto each singular vector
    Nk = 2048
    k = (np.arange(Nk//2 + 1) / Nk) * (2 * np.pi / xs) if xs else np.arange(Nk//2 + 1)
    k_axename = 'k/ks' if not xs else 'k'

    projections = ifft(U, Nk, axis=0)
    projections = projections[:Nk//2 + 1, :, :]

    projections_sum = np.zeros((Nk//2 + 1, Nf//2 + 1, Nemit))
    #for kemit in range(Nemit):
    #    for ii in range(Nf//2 + 1):
    #        projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]) ** 2
            
    #kemit_values = [1,2,3]
    kemit_values = [0,1,2]
    for kemit in kemit_values:
        for ii in range(Nf//2 + 1):
            max_value = 1#np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2
            

    projections_sum = np.abs(np.mean(projections_sum, axis=2))

    if issaving:
        plt.savefig(name + '_svd')

    return f, k, projections_sum

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


def read_data(path2data):
    """ Read path to folder where data are stored, 
    Extracte each stream and create a list of streams """
    
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

def convert_to_utc_times(trace, time_vector):
    """ Create a UTC-time time vector for a given trace """
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]


def rename_traces(stream, geophones_dict):
    """ Rename geophone traces using geophone_table """
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
            print(f"Warning: No entry found for {last_five_digits} in the geophones table.")

    # Sort traces based on the new station codes
    sorted_stream = sorted(stream, key=lambda trace: trace.stats.station)

    return sorted_stream

def sort_key(trace):
    """ Create a sorting key based on stream stations """
    return trace[0].stats.station

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

####################################################
#%% Plot a single channel for a given geophone
####################################################
channel = 0
idx_geophone = 3

current_stream = seismic_data_streams[(idx_geophone - 1)*3 + channel]
print(current_stream[0])
normalized_data = current_stream[0].data * scale_factor

fig, ax  = plt.subplots(figsize = fig_size)
ax.plot(first_trace.times(),normalized_data)
ax.set_ylabel('$V \; \mathrm{(m.s^{-1})}$')
ax.set_xlabel('$t \; \mathrm{(s)}$')


#%%
figname = fig_folder + 'Scaled_signal_channel_' + channel_correspondance[channel] + 'idx_geophone_' + str(idx_geophone) +'_' + acqu_numb + '_gain_' + str(gain)
plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality, bbox_inches = 'tight')

##############################################################################################
#%% ---------------------- Plot all channels for a single geophone ---------------------------
##############################################################################################

idx_geophone = 2 # index of the geophone to be plotted

selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]
plt.rcParams.update({
    "text.usetex": True})
fig, ax = plt.subplots(3,figsize = fig_size)
for k, chan in enumerate([2,0,1]):
    current_stream = selected_streams[chan]
    print(current_stream[0])
    ax[k].plot(first_trace.times(),current_stream[0].data / max(np.abs(current_stream[0].data) ),label = channel_correspondance[chan])
    # ax[k].plot(first_trace.times(),current_stream[0].data*scale_factor,label = channel_correspondance[chan])
    ax[k].legend()
    
    ax[k].set_ylabel(r'$V/V_{max}$')
    ax[k].set_ylim([-1.1, 1.1])
    # fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
    # # # Create an instance of ZoomHandler
    # # zoom_handler = ZoomHandler(ax, time_vector, data_vector)
    # # fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
    # # # Enable interactive labels using mplcursors
    # # mplcursors.cursor(hover=True)
    # # Adjust the backend to make it work better in Spyder
    # plt.ion()
    # plt.show(block=True)  # Use block=True to make it work better in Spyder

ax[2].set_xlabel('$t \; \mathrm{(s)}$')
fig.tight_layout()

#%% Set xlim for all subplots

Mydate = date(2024,2,11)
time_start = datetime.combine(Mydate,time(18,47,7))
time_end = datetime.combine(Mydate,time(18,47,35))
xlimits = [time_start, time_end]

for k in range(0,3):
    ax[k].set_xlim([time_start, time_end])
    
#%% Save current figure
fig.tight_layout()

figname = 'Stream_ZEN_geophone_' + str(idx_geophone) + '_acqu_1'
figname = fig_folder + figname

plt.savefig(figname + '.pdf',dpi = img_quality,bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality,bbox_inches = 'tight')

###################################################
#%% OPTIONAL - Show the integrated signal on a given channel 
###################################################

channel = 2 
idx_geophone = 2
facq_t = 1000 
current_stream = seismic_data_streams[(idx_geophone - 1)*3 + channel]
displacement = integrate.cumulative_trapezoid(current_stream[0].data*scale_factor)/facq_t

fig, ax = plt.subplots(figsize = (16,10))
ax.plot(datetime_values[:-1],displacement)
ax.set_xlabel(r'UTC time')
ax.set_ylabel(r'$\delta \; \mathrm{(m)}$')

fig, ax = plt.subplots(figsize = (16,10))
ax.plot(datetime_values,current_stream[0].data * norm_factor)

#%%
Mydate = date(2024,2,11)
time_start = datetime.combine(Mydate,time(20,30,0))
time_end = datetime.combine(Mydate,time(20,40,0))
ax.set_xlim([time_start, time_end])

figname = fig_folder + 'Vertical_velocity_0211_G2'
plt.savefig(figname + '.pdf')
plt.savefig(figname + '.png')

###############################################################
#%% OPTIONAL -Separated subplots for one single channel and all geophones 
###############################################################

geophone_location = [i*3 for i in range(0,16,1)]

t0 = 68.6
t1 = 70.1

# mask = np.where(np.logical_and(first_trace.times() >= t0, first_trace.times() <= t1))
# selected_time = first_trace.times()[mask]
# selected_time = selected_time - selected_time[0]

time_array = first_trace.times() - t0
color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

fig, ax = plt.subplots(len(selected_indices),1, figsize = fig_size)
for idx_geo,k in enumerate(selected_indices): # k : index of the stream 
    current_stream = seismic_data_streams[k]
    print(current_stream[0])
    idx_subplot = 15 - idx_geo
    if idx_geo <= 9:
        color_idx = idx_geo
    else: 
        color_idx = idx_geo % 10 # take remainer of euclidian division 
    ax[idx_subplot].plot(time_array,current_stream[0].data / max(np.abs(current_stream[0].data) ) + idx_geo*3, color = color_cycle[color_idx])
    ax[idx_subplot].set_xlim([0.1,t1-t0])
    ax[idx_subplot].tick_params(bottom = False)
    ax[idx_subplot].set_yticks([geophone_location[idx_geo]])
    ax[idx_subplot].spines[['top', 'right', 'bottom' , 'left']].set_visible(False)

ax[15].spines['bottom'].set_visible(True)
ax[15].tick_params(bottom = True)
# ax2 = ax[15].twinx()
# ax2.set_yticks(np.array([-1,1]) * max(np.abs(seismic_data_streams[selected_indices[0]][0].data)) * scale_factor)
# # ax[15].set_xticks([0,0.2,0.4,0.6,0.8,1.0,1.2,1.4])
# ax2.spines[['top','left']].set_visible(False)

plt.tight_layout()

# set labels 
ax[15].set_xlabel(r'$t \; \mathrm{(s)}$')
ax[8].set_ylabel(r'Distance to G1 (m)')
# ax2.set_ylabel('$V_z \; \mathrm{(m.s^{-1})}$')

#%%
############################################################################################
figname = fig_folder + 'Superposition_geophones_Z3_S3_fancy'

plt.savefig(figname + '.pdf',dpi = img_quality,bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = img_quality,bbox_inches = 'tight')

#%%----------------------- PLOTTING SELECTED CHANNEL FOR WHOLE RECORDING ------------------ 
##############################################################################################

fig, ax = plt.subplots(figsize = fig_size)
for k in selected_indices:
    current_stream = seismic_data_streams[k]
    print(current_stream[0])
    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) ) + k, 
            label=f"Stream {k}")
    

#fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
# Create an instance of ZoomHandler
zoom_handler = ZoomHandler(ax, time_vector, data_vector)
fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder


#%% Select xlim of the plot 
Mydate = date(2024,2,11)
time_start = datetime.combine(Mydate,time(18,47,7))
time_end = datetime.combine(Mydate,time(18,47,35))
ax.set_xlim([time_start, time_end])

#%% Once correclty zoomed, we can save the figure 
figname = 'Streams_all_geophones_' + channel_correspondance[channel] + acqu_numb
figname = fig_folder + figname

plt.savefig(figname + '.pdf')
plt.savefig(figname + '.png')

#################################################################################
#%%----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
#################################################################################
# Select signal segments manually for each track and each source 

# initial time of each segment, dimensions : 
    #1 = direction 
    #2 = channel (0:E,1:N,2:Z)
    #3 = source 
    
t0_segments = np.empty((2,3,3))

# direction 1
 
# source Z 
t0_segments[0,2,:] = [UTCDateTime("2024-02-11T18:44:29.20"), # S101 Z3
      UTCDateTime("2024-02-11T18:44:58.45"), # S102 Z3
      UTCDateTime("2024-02-11T18:45:31.05")] # S103 Z2

# source E
t0_segments[0,0,:] = [UTCDateTime("2024-02-11T18:44:37.32"), # S101 E2
      UTCDateTime("2024-02-11T18:45:08.30"), # S102 E3
      UTCDateTime("2024-02-11T18:45:38.42")] # S103 E1

# source N
t0_segments[0,1,:] = [UTCDateTime("2024-02-11T18:44:48.14"), # S101 N2
      UTCDateTime("2024-02-11T18:45:16.48"), # S102 N3
      UTCDateTime("2024-02-11T18:45:48.15")] # S103 N2

# direction 2

# source Z 
t0_segments[1,2,:] = [UTCDateTime("2024-02-11T18:47:08.80"), # S101 Z1
      UTCDateTime("2024-02-11T18:47:40.80"), # S102 Z3
      UTCDateTime("2024-02-11T18:48:09.50")] # S103 Z2

# source E
t0_segments[1,0,:] = [UTCDateTime("2024-02-11T18:47:16.90"), # S101 E2
      UTCDateTime("2024-02-11T18:47:51.00"), # S102 E3
      UTCDateTime("2024-02-11T18:48:16.90")] # S103 E1

# source N
t0_segments[1,1,:] = [UTCDateTime("2024-02-11T18:47:30.30"), # S101 N2
      UTCDateTime("2024-02-11T18:47:57.90"), # S102 N3
      UTCDateTime("2024-02-11T18:48:20.60")] # S103 N2

########################
#%% Observation of the selected signal 
########################

signal_length = 1.5 # duration in seconds

### Select direction on which we extract data #0 = dir1, #1 = dir2
direction = 0
### Select the channel on which we extract data
channel = 0 # #0 = E, 1 = N, #2 = Z 
### Select source index 0, 1 or 2 (S1, S2 or S3)
source_index = 0

# selected indices of the corresponding channel 
selected_indices = range(channel, len(seismic_data_streams),3)

t1 = t0_segments[direction,channel,source_index]
# if channel == 0:
#     t1 = tE[source_index]
# elif channel == 1:
#     t1 = tN[source_index]
# else :
#     t1 = tZ[source_index]
    
t2 = t1 + signal_length
num_samples = int((t2 - t1) * first_trace.stats.sampling_rate)
# time_vector = np.linspace(t1.timestamp, t2.timestamp, num_samples)
time_vector = np.linspace(t1, t2, num_samples)
UTC_time_vector = [datetime.utcfromtimestamp(t) for t in time_vector]
# Create a matrix to store the seismic data
seismic_matrix = np.zeros((len(selected_indices), num_samples))

# Extract the relevant data for each trace and populate the matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)):
    stream = seismic_data_streams[stream_index]
    for j, trace in enumerate(stream):
        start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
        end_sample = start_sample + num_samples
        seismic_matrix[i, :] = trace.data[start_sample:end_sample]

fig, ax = plt.subplots(figsize = (16,10))
for k in range(seismic_matrix.shape[0]):
    ax.plot(UTC_time_vector, seismic_matrix[k, :] / max(np.abs(seismic_matrix[k, :])) + geophones_spacing*k, label=f"Stream {k}")
ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Distance to G1 $(\mathrm{m^{-1}})$')
# ax.set_title('Seismic Data for Selected Streams, channel ' + channel_correspondance[channel])
#ax.legend()
plt.show()

#%% Save the current figure 
figname = fig_folder + 'Signal_evolution_' + channel_correspondance[channel] + str(source_index)
plt.savefig(figname + '.pdf',dpi = 1200)
plt.savefig(figname + '.png',dpi = 1200)

########################################################################
#%% OPTIONAL - FREQUENCY SPECTRUM VS WAVENUMBER for the selected signal
########################################################################

Nfft_t= 1024
Nfft_k = 1024
 
fs = 1000
ks = 2*np.pi/geophones_spacing
space_vector = np.arange(0,seismic_matrix.shape[0],geophones_spacing)

# Define the spatial and temporal frequencies
wavenum = np.fft.fftfreq(Nfft_k, d=1/ks)
freq = np.fft.fftfreq(Nfft_t, d=1/fs)

# Create a 2D Fourier transform of the seismic matrix
seismic_fft = np.fft.fft2(seismic_matrix, (Nfft_k, Nfft_t))

# Shift zero frequency component to the center
seismic_fft_shifted = np.fft.fftshift(seismic_fft)

fig, ax = plt.subplots()
# Plot the 2D Fourier transform
plt.imshow(np.transpose(np.abs(seismic_fft_shifted)), extent=(wavenum.min(), wavenum.max(),freq.min(), freq.max()), aspect='auto', cmap='viridis')
#plt.colorbar(label='Amplitude')
plt.ylabel(r'$f\; \mathrm{(Hz)}$')
plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
plt.title('2D Fourier Transform of Seismic')
plt.show()
ax.set_ylim([0 , fs/4])


##################################################################################
#%%----------------------- CONSTRUCTING MATRIX FOR FK WITH SVD ------------------
##################################################################################


### Select direction on which we extract data #0 = dir1, #1 = dir2
direction = 0
# choose channel from which we extract signals : #0 = E, #1 = N, #2 = Z
channel = 0

t1_values = t0_segments[direction,channel,:]
# if channel == 0 :
#     t1_values = tE
# elif channel == 1:
#     t1_values = tN
# else:
#     t1_values = tZ
    

# Create a 3D matrix to store the seismic data
num_traces = len(seismic_data_streams)
num_samples = int(signal_length * first_trace.stats.sampling_rate)  

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((len(selected_indices), third_dim_len, num_samples)) # matrix dim : geo_indices x signal_number x time

# Extract the relevant data for each trace and populate the 3D matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)):
    stream = seismic_data_streams[stream_index]
    
    for t1_index, t1 in enumerate(t1_values):
        for j, trace in enumerate(stream):
            # start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
            start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index, :] = trace.data[start_sample:end_sample]

print('Matrix computed for channel : ' + channel_correspondance[channel])

#%% Plot the data (superposition of three sources)

fig, ax = plt.subplots(figsize = fig_size)
for k in range(seismic_matrix.shape[0]):
    for t1_index in range(seismic_matrix.shape[1]):
        ax.plot(UTC_time_vector, seismic_matrix[k, t1_index, :] / max(np.abs(seismic_matrix[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams at Different t1 Values')
plt.show()

#######################################################################
#%% OPTIONAL ---------------------- INTERPOLATING 3D MATRIX ------------------
#######################################################################

# Interpolate along the first dimension



xr = np.arange(0, 46, 3)
xs = np.arange(5, 12, 3)
xr_new = np.arange(0, 46, 1)
t_new = time_vector


xx, tt = np.meshgrid(xr, t_new, sparse=True)


seismic_matrices_interp = []

for k in range(seismic_matrix.shape[1]):  # Loop over the second dimension of seismic_matrix
    data_to_interpolate = np.transpose(seismic_matrix[:, k, :])
    f = interpolate.interp2d(xr, t_new, data_to_interpolate, kind='cubic')
    seismic_matrix_interp = f(xr_new, t_new)
    
    # Append the interpolated matrix to the list
    seismic_matrices_interp.append(seismic_matrix_interp)

# Convert the list to a 3D numpy array
seismic_matrices_interp = np.stack(seismic_matrices_interp, axis=1)
seismic_matrices_interp = np.transpose(seismic_matrices_interp, (2, 1,0))


fig, ax = plt.subplots()
for k in range(seismic_matrices_interp.shape[0]):
    for t1_index in range(seismic_matrices_interp.shape[1]):
        ax.plot(time_vector, seismic_matrices_interp[k, t1_index, :] / max(np.abs(seismic_matrices_interp[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams at Different t1 Values')
plt.show()

#%% Select points to determine shear velocity 

points_shear = plt.ginput(10, timeout=-1)
k_mode, f_mode = zip(*points_shear)

#%%
# fit by a 1D polynome
deg = 1
p = np.polyfit(k_mode,f_mode,deg)

C_shear = p[0]*2*np.pi
print(C_shear)


#%% Select points to determine longitudinal velocity 

points_longi = plt.ginput(10, timeout=-1)
f_mode, k_mode = zip(*points_longi)

#%%
# fit by a 1D polynome
deg = 1
p = np.polyfit(f_mode,k_mode,deg)

C_longi = 2*np.pi*p[0]
print(C_longi)

#%%
# C_shear = np.mean([2*np.pi*200/0.85])
C_shear = np.mean([2*np.pi*220/0.5])
# C_longi = np.mean([2*np.pi*100/0.59])
C_longi = np.mean([2*np.pi*220/1.0])

#%%
rho_ice = 917
nu = 1-2*(C_shear/C_longi)**2
E = rho_ice*C_longi**2*(1-nu**2)

print(E*1e-9)
print(nu)









