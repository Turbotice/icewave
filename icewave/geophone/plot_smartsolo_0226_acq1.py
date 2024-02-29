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
from datetime import datetime
import numpy as np
from scipy import interpolate




path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0226/Geophones'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0002'
channel = 0 #0 for E, 1 for N, 2 for Z. 


#----------------------- FUNCTIONS ------------------ 

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
first_stream = seismic_data_streams[1]
first_trace = first_stream[0]
time_vector = convert_to_utc_times(first_trace, first_trace.times())
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector]
data_vector = first_trace.data
start_time_utc = UTCDateTime(first_trace.stats.starttime)
fs = first_trace.stats.sampling_rate

 
selected_indices = range(channel, len(seismic_data_streams),3)

#----------------------- PLOTTING SELECTED CHANNEL FOR WHOLE RECORDING ------------------ 
fig, ax = plt.subplots()
for k in selected_indices:
    current_stream = seismic_data_streams[k]
    print(current_stream[0])
    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) )+ k, 
            label=f"Stream {k}")
fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
# Create an instance of ZoomHandler
zoom_handler = ZoomHandler(ax, time_vector, data_vector)
fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder



#----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
# Define time range for "zooming"




t1 = UTCDateTime("2024-02-26T18:18:52.7") # S104 Z3
t1 = UTCDateTime("2024-02-26T18:19:54.8") # S105 Z3
t1 = UTCDateTime("2024-02-26T18:20:37.75") # S106 Z3

t1 = UTCDateTime("2024-02-26T18:19:09.24") # S104 E?
t1 = UTCDateTime("2024-02-26T18:20:42.59") # S104 E?
t1 = UTCDateTime("2024-02-26T18:20:53.44") # S104 E?

t1 = UTCDateTime("2024-02-26T18:19:17.84") # S104 N2
t1 = UTCDateTime("2024-02-26T18:20:15.42") # S104 N
t1 = UTCDateTime("2024-02-26T18:21:01.86") # S104 N


t1 = UTCDateTime("2024-02-26T18:16:43.62") # S101 N2
t1 = UTCDateTime("2024-02-26T18:17:16.70") # S102 N2
t1 = UTCDateTime("2024-02-26T18:17:58.1") # S103 N3


t1 = UTCDateTime("2024-02-26T18:16:41.75") # S103 N3
t1 = UTCDateTime("2024-02-26T18:16:36.15") # S103 N3

signal_length = 1.5 # duration in seconds
t2 = t1 + signal_length
num_samples = int((t2 - t1) * first_trace.stats.sampling_rate)
time_vector = np.linspace(t1.timestamp, t2.timestamp, num_samples)

# Create a matrix to store the seismic data
num_traces = len(seismic_data_streams)



seismic_matrix = np.zeros((len(selected_indices), num_samples))
# Extract the relevant data for each trace and populate the matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)):
    stream = seismic_data_streams[stream_index]
    for j, trace in enumerate(stream):
        start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
        end_sample = start_sample + num_samples
        seismic_matrix[i, :] = trace.data[start_sample:end_sample]

fig, ax = plt.subplots()
for k in range(seismic_matrix.shape[0]):
    ax.plot(time_vector, seismic_matrix[k, :] / max(np.abs(seismic_matrix[k, :])) + 3*k, label=f"Stream {k}")
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
space_vector = np.arange(0,seismic_matrix.shape[0],geohpones_spacing)

# Define the spatial and temporal frequencies
wavenum = np.fft.fftfreq(Nfft_k, d=1/ks)
freq = np.fft.fftfreq(Nfft_t, d=1/fs)

# Create a 2D Fourier transform of the seismic matrix
seismic_fft = np.fft.fft2(seismic_matrix, (Nfft_k, Nfft_t))

# Shift zero frequency component to the center
seismic_fft_shifted = np.fft.fftshift(seismic_fft)

# Plot the 2D Fourier transform
plt.imshow((np.abs(seismic_fft_shifted)), extent=(freq.min(), freq.max(), wavenum.min(), wavenum.max()), aspect='auto', cmap='viridis')
#plt.colorbar(label='Amplitude')
plt.xlabel('Temporal Frequency (Hz)')
plt.ylabel('Spatial Frequency')
plt.title('2D Fourier Transform of Seismic')
plt.show()









#----------------------- CONSTRUCTING MATRIX FOR FK WITH SVD ------------------




t1_values = [UTCDateTime("2024-02-26T18:18:52.7"), 
              UTCDateTime("2024-02-26T18:19:54.8"), 
              UTCDateTime("2024-02-26T18:20:37.75")]  #sources Z dir2

t1_values = [UTCDateTime("2024-02-26T18:19:09.24"), 
              UTCDateTime("2024-02-26T18:20:42.59"), 
              UTCDateTime("2024-02-26T18:20:53.44")]  #sources E dir2


t1_values = [UTCDateTime("2024-02-26T18:19:17.84"), 
              UTCDateTime("2024-02-26T18:20:15.42"), 
              UTCDateTime("2024-02-26T18:21:01.86")]  #sources N dir2


t1_values = [UTCDateTime("2024-02-26T18:16:43.62"), 
              UTCDateTime("2024-02-26T18:17:16.70"), 
              UTCDateTime("2024-02-26T18:17:58.1")]  #sources N dir1





# Create a 3D matrix to store the seismic data
num_traces = len(seismic_data_streams)
num_samples = int(signal_length * first_trace.stats.sampling_rate)  

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((len(selected_indices), third_dim_len, num_samples))

# Extract the relevant data for each trace and populate the 3D matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)):
    stream = seismic_data_streams[stream_index]
    
    for t1_index, t1 in enumerate(t1_values):
        for j, trace in enumerate(stream):
            start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index, :] = trace.data[start_sample:end_sample]

# Plot the data
fig, ax = plt.subplots()
for k in range(seismic_matrix.shape[0]):
    for t1_index in range(seismic_matrix.shape[1]):
        ax.plot(time_vector, seismic_matrix[k, t1_index, :] / max(np.abs(seismic_matrix[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams at Different t1 Values')
plt.show()



#----------------------- INTERPOLATING 3D MATRIX ------------------

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



C_shear = np.mean([2*np.pi*48/0.417])
C_longi = np.mean([2*np.pi*100/0.59])


rho_ice = 917
nu = 1-2*(C_shear/C_longi)**2
E = rho_ice*C_longi**2*(1-nu**2)










