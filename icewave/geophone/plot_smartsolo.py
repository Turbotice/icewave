#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:38:04 2024

@author: moreaul
"""
#%matplotlib qt
import os
import shutil
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime
import numpy as np



path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0210/Geophones'
acqu_numb = '0001'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

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

def read_data(path2data):
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

class ZoomHandler:
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

def convert_to_utc_times(trace, time_vector):
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]







# Read geophones table into a dictionary
geophones_dict = {}
with open(geophones_table_path, 'r') as table_file:
    next(table_file)  # Skip the header line
    for line in table_file:
        num, geo_sn = line.strip().split('\t')
        geophones_dict[geo_sn] = num

def read_data(path2data):
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
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]

# Load geophones table into a dictionary
geophones_dict = {}
with open('/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table', 'r') as table_file:
    for line in table_file:
        if line.startswith("Num"):
            continue  # Skip header line
        num, geo_sn = line.strip().split()
        # Assuming the last five digits are relevant for comparison
        last_five_digits = geo_sn[-5:]
        geophones_dict[last_five_digits] = num


def rename_traces(stream, geophones_dict):
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





# Read MiniSEED file directly
seismic_data_streams, miniseed_files = read_data(path2data +'/' +acqu_numb)

# Iterate over streams and rename traces
for i, stream in enumerate(seismic_data_streams):
    seismic_data_streams[i] = rename_traces(stream, geophones_dict)

def sort_key(trace):
    return trace[0].stats.station

# Sort the seismic_data_streams based on the custom sorting function
seismic_data_streams = sorted(seismic_data_streams, key=sort_key)

# Find the largest start time (t0) and the smallest end time (t1) among all streams
t0 = max([stream[0].stats.starttime for stream in seismic_data_streams])
t1 = min([stream[-1].stats.endtime for stream in seismic_data_streams])

# Truncate all streams between t0 and t1
for i, stream in enumerate(seismic_data_streams):
    for trace in stream:
        if trace.stats.starttime < t0:
            trace.trim(t0, t1, pad=True, fill_value=0)
        elif trace.stats.endtime > t1:
            trace.trim(t0, t1, pad=True, fill_value=0)



# Extract time_vector and data_vector for the first stream outside the loop
first_stream = seismic_data_streams[1]
first_trace = first_stream[0]
time_vector = convert_to_utc_times(first_trace, first_trace.times())
data_vector = first_trace.data
start_time_utc = UTCDateTime(first_trace.stats.starttime)
fs = first_trace.stats.sampling_rate

num_traces_to_plot = 48
# Calculate datetime values only once outside the loop
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector]


channel_nb = 0 #0 for E, 1 for N, 2 for Z.  
fig, ax = plt.subplots()
for k in range(0, num_traces_to_plot,3):
    stream_index =  k+channel_nb  # Adjust index to match Python's 0-based indexing
    current_stream = seismic_data_streams[stream_index]
    print(current_stream[0])

    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) )+ k, 
            label=f"Stream {stream_index}")

fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)


# Create an instance of ZoomHandler
zoom_handler = ZoomHandler(ax, time_vector, data_vector)
fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder




## EXTRACTING AND PLOTTING SIGNALS IN TIME OF INTEREST 
# Define the desired time range
t1 = UTCDateTime("2024-02-10T12:50:52") 
t2 = t1 + 1.5

# Create a matrix to store the seismic data
num_traces = len(seismic_data_streams)
num_samples = int((t2 - t1) * first_trace.stats.sampling_rate)

seismic_matrix = np.zeros((num_traces, num_samples))

# Extract the relevant data for each trace and populate the matrix
for i, stream in enumerate(seismic_data_streams):
    for j, trace in enumerate(stream):
        start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
        end_sample = start_sample + num_samples
        seismic_matrix[i, :] = trace.data[start_sample:end_sample]



time_vector = np.linspace(t1.timestamp, t2.timestamp, num_samples)


selected_indices = range(2, num_traces_to_plot,3)

fig, ax = plt.subplots()

for k in selected_indices:
    ax.plot(time_vector, seismic_matrix[k, :] / max(np.abs(seismic_matrix[k, :])) + k, label=f"Stream {k}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams')
#ax.legend()
plt.show()



## PLOTTING WAVENUMBER VS FREQUENCY SPECTRUM
Nfft_t= 1024
Nfft_k = 1024
geohpones_spacing = 3
fs = 1000
ks = 2*np.pi/geohpones_spacing
space_vector = np.arange(0,48,geohpones_spacing)

# Define the spatial and temporal frequencies
freq_spatial = np.linspace(0,ks/2,Nfft_k)
freq_temporal = np.linspace(0,fs/2,Nfft_t)

freq_spatial = np.fft.fftfreq(Nfft_k, d=1/ks)
freq_temporal = np.fft.fftfreq(Nfft_t, d=1/fs)


# Create a 2D Fourier transform of the seismic matrix
seismic_fft = np.fft.fft2(seismic_matrix, (Nfft_k, Nfft_t))

# Shift zero frequency component to the center
seismic_fft_shifted = np.fft.fftshift(seismic_fft)

# Plot the 2D Fourier transform
plt.imshow(np.abs(seismic_fft_shifted), extent=(freq_temporal.min(), freq_temporal.max(), freq_spatial.min(), freq_spatial.max()), aspect='auto', cmap='viridis')
#plt.colorbar(label='Amplitude')
plt.xlabel('Temporal Frequency (Hz)')
plt.ylabel('Spatial Frequency')
plt.title('2D Fourier Transform of Seismic')
plt.show()


# Define three different values of t1
t1_values = [UTCDateTime("2024-02-10T12:50:52"), UTCDateTime("2024-02-10T12:50:54.75"), UTCDateTime("2024-02-10T12:50:58")]

# Calculate corresponding t2 values
t2_values = [t1 + 1.5 for t1 in t1_values]

# Create a 3D matrix to store the seismic data
num_traces = len(seismic_data_streams)

# Determine the maximum number of samples among all time windows
num_samples_max = max(int((t2 - t1) * seismic_data_streams[0][0].stats.sampling_rate) for t1, t2 in zip(t1_values, t2_values))

seismic_3D_matrix = np.zeros((len(t1_values), num_traces, num_samples_max))

# Extract the relevant data for each trace, for each t1 value, and populate the 3D matrix
for k, (t1, t2) in enumerate(zip(t1_values, t2_values)):
    num_samples = int((t2 - t1) * seismic_data_streams[0][0].stats.sampling_rate)
    
    for i, stream in enumerate(seismic_data_streams):
        for j, trace in enumerate(stream):
            start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_3D_matrix[k, i, :num_samples] = trace.data[start_sample:end_sample]

fig, axs = plt.subplots(len(t1_values), 1, figsize=(8, 6), sharex=True, sharey=True)

# Plot seismic data for each time window
for k, (t1, t2) in enumerate(zip(t1_values, t2_values)):
    seismic_matrix_slice = seismic_3D_matrix[k, :, :]
    
    # Extract time vector
    num_samples = len(seismic_matrix_slice[0])
    time_vector = np.linspace(t1.timestamp, t2.timestamp, num_samples)
    
    # Plot seismic data for selected traces
    selected_indices = range(2, num_traces, 3)
    for trace_idx in selected_indices:
        axs[k].plot(time_vector, seismic_matrix_slice[trace_idx, :] / max(np.abs(seismic_matrix_slice[trace_idx, :])) + trace_idx, label=f"Stream {trace_idx} - Time Window {k+1}")

    axs[k].set_title(f'Seismic Data - Time Window {k+1}')

axs[-1].set_xlabel('Time (UTC)')
axs[len(t1_values)//2].set_ylabel('Normalized Seismic Data')
plt.tight_layout()
plt.show()


# Accumulate Fourier transform values for each time window
average_fft = np.zeros((Nfft_k, Nfft_t))

for k, t1 in enumerate(t1_values):
    seismic_matrix_slice = seismic_3D_matrix[k, :, :]
    
    seismic_fft = np.fft.fft2(seismic_matrix_slice, (Nfft_k, Nfft_t))
    seismic_fft_shifted = np.fft.fftshift(seismic_fft)

    # Accumulate Fourier transform values
    average_fft += np.abs(seismic_fft_shifted)

# Calculate the average
average_fft /= len(t1_values)

# Plot the averaged spectrum
fig, ax_avg = plt.subplots(figsize=(8, 6))
ax_avg.imshow(average_fft, extent=(freq_temporal.min(), freq_temporal.max(), freq_spatial.min(), freq_spatial.max()), aspect='auto', cmap='viridis')
ax_avg.set_title('Averaged 2D Fourier Transform')

ax_avg.set_xlabel('Temporal Frequency (Hz)')
ax_avg.set_ylabel('Spatial Frequency')

plt.tight_layout()
plt.show()
