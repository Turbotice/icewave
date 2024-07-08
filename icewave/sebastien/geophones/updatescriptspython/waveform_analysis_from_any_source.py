#!/usr/bin/env python
# coding: utf-8

# In[1]:


#get_ipython().system(' pip install mplcursors')
#get_ipython().system(' pip install obspy')
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
# get_ipython().run_line_magic('matplotlib', 'qt')
plt.close('all')

# ## Input

# In[50]:


year = '2024'
date = '0211' #date format, 'mmdd'
acqu_numb = '0003' #acquisition number 

direction = 2 # 1 ou 2 
composante = 'Z' #Z , E or N -> direction de la source
channel = 2  # 0 for E, 1 for N, 2 for Z. 
flexure_wave = 1 # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = 0
#files need to be organised as: data/0210/Geophones/0001/minised files


path2data = 'F:/Rimouski_2024/Data/2024'
geophones_table_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Geophones/geophones_table'

# ## Fonctions

# In[51]:


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


#----------------------------------------------------------------------------------------------
def fn_svd ( signals , fs , xs , rang , name , issaving , *varargin ) : #varagin ??
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

    Nreceiv, Nt, Nemit = signals.shape # (16, 1000, 3) (geophones, nb de valeurs, sources) # Nt not used

    # Time domain fft
    Nf = 2048
    f = (np.arange(Nf) / Nf) * (fs if fs else 1)
    f_axename = 'f/fs' if not fs else 'f'
    SIGNALS = fft(signals, Nf, axis=1)
    SIGNALS = SIGNALS[:, :Nf + 1, :]

    # svd
    # Creqtion matrice U S V,  D???
    U = np.zeros((Nreceiv, Nf, Nemit), dtype=complex) # 16, 2048, 3
    S = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    V = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    D = np.zeros((Nemit, Nf), dtype=complex)

    for ii in range(Nf):
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False)
        D[:, ii] = np.diag(S[:, ii, :])
    for ne in range(Nemit):
        titi = 20 * np.log10(np.abs(D[ne, :]) / np.max(np.abs(D[0, :])))
        #plt.plot(f, titi, label=f'Slice {ne}')

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

    for ne in range(Nemit):
        titi = 20 * np.log10(D[ne, :] / np.max(D[0, :]))
        idx = np.where(titi <= sigmacutsup)[0]
        U[:, idx, ne] = 0

    if issaving:
        plt.savefig(name + '_sv')

    # projection onto each singular vector
    Nk = 2048
    k = (np.arange(Nk) / Nk) * (2 * np.pi / xs)  if xs else np.arange(Nk + 1)
    k_axename = 'k/ks' if not xs else 'k'

    projections = ifft(U, Nk, axis=0)#np.fft.fftshift(ifft(U, Nk, axis=0), axes=0)
    projections_sum = np.zeros((Nf, Nk, Nemit))

    for kemit in rang:
        for ii in range(Nf):
            max_value = 1  # np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2

    projections_sum = np.abs(np.mean(projections_sum, axis=2))

    return f, k, projections_sum
#----------------------------------------------------------------------------------------------


# ## Plot of the 3 wave signals 

# In[52]:


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
seismic_data_streams, miniseed_files = read_data (path2data +'/'+ date + '/Geophones/' + acqu_numb)

# Iterate over streams and rename traces
for i, stream in enumerate(seismic_data_streams):
    seismic_data_streams[i] = rename_traces(stream, geophones_dict)
    

def sort_key(trace):
    return trace[0].stats.station
# Sort the seismic_data_streams based on the custom sorting function
seismic_data_streams = sorted(seismic_data_streams, key=sort_key)

first_stream = seismic_data_streams[1]  #first stream is a obspy.core.trace object
first_trace = first_stream[0] #first stream is SS.01.SSN..SSN

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


fs = first_trace.stats.sampling_rate # Needed, fn_svd
selected_indices = range(channel, len(seismic_data_streams), 3) # Needed for seismic_matrix

#%%

# If "Unknown format for file 0005\._01.0005.2024.02.11.20.06.55.000.E.miniseed" appears
# Create a NEW 000n file only with the miniseed files. 
# assign a string to the channel values 
signal_length = 1 # duration in seconds 
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}
ch = channel_dic[channel]

# assign a string to S values depending on the direction
if direction == 1 :
    S1 = '101' 
    S2 = '102' 
    S3 = '103'
if direction == 2: 
    S1 = '104' 
    S2 = '105' 
    S3 = '106'

# get the pck file that assigns date + year and the t1_to_time... file
base = 'C:/Users/sebas/git/icewave/sebastien/geophones/updatescriptspython/'
var_path = base + 'variable_to_dict.pkl'
with open(var_path, 'rb') as m:
    variable_to_dict = pickle.load(m)
if variable_to_dict.get(date+'_'+year) is None:
    messagebox.showerror("Error", "Date not found")
    raise RuntimeError("Date not found")
    
# get the pkl file assigned to the date and year value
pkl_path = base + variable_to_dict.get(date+'_'+year) + ".pkl"
with open(pkl_path, 'rb') as f:
    loaded_data = pickle.load(f)
# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)

# Get the value of "time_vector" for the x abscisse in plot, 1000 values between t1 and ta

ta = t1 + signal_length 
num_samples = int(signal_length * first_trace.stats.sampling_rate) # Number of points gereted in tim_vector (~1000)
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

t1_values = [t1,t2,t3]
# Create a 3D matrix to store the seismic data
# Dynamically determine the length of the third dimension, number of sources
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((len(selected_indices), third_dim_len, num_samples)) # create the size of seismic_matrix

# Extract the relevant data for each trace and populate the 3D matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)): # stream index ??
    stream = seismic_data_streams[stream_index]
    for t1_index, t1 in enumerate(t1_values): #t1_index = indice de t1_values (0,1,2) t1 = 3 different t0
        for j, trace in enumerate(stream): #stream cell above ? in loop, stream = seismic_data_streams[stream_index]
            start_sample = int((t1 - trace.stats.starttime)
                               * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index,
                           :] = trace.data[start_sample:end_sample]

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
             f"Channel: {ch}, Source: {composante}, Direction: {direction}", fontsize = 25)
ax.set_ylabel('Normalized Seismic Data', fontsize = 20)
ax.tick_params(axis='x',labelsize = 20)
ax.tick_params(axis='y',labelsize = 20)
plt.show()
#plt.grid(True)





# Calling fn_SVD
rang = [0,1,2]
geophones_spacing = 3 # in meters
signals = np.transpose(seismic_matrix, (0, 2, 1))

f, k, FK = fn_svd(signals, fs, geophones_spacing , rang ,'ExampleName', 0, 'threshold',-90) #plot valeurs  singuliere/freq
F, K = np.meshgrid(f, k)

# Normalize FK based on direction
if direction == 1:
    FK_normalized = FK / np.max(FK)
elif direction == 2:
    FK_normalized = np.flipud(FK / np.max(FK))

# Unwrapping Spectrum
FK_uwp = np.vstack((FK_normalized, FK_normalized))
k_uwp = np.linspace(0, 2 * max(k), FK_uwp.shape[0])
F, K_uwp = np.meshgrid(f, k_uwp)

#%%
#-----------------------------Plotting--------------------------------------------------
vmin = 0     # Minimum value for colormap
vmax = 1     # Maximum value for colormap (adjust as needed)
xmin = 0  # Minimum frequency
xmax = 500  # Maximum frequency
ymin = 0  # Minimum wavenumber
ymax = 4  # Maximum wavenumber

fig, axes = plt.subplots(1, 1, figsize=(15, 10))
ax1 = axes  # No need for indexing when there's only one subplot
c1 = ax1.contourf(F, K_uwp, FK_uwp, cmap='gnuplot2', vmin=vmin, vmax=vmax)

ax1.set_xlabel('Frequency (Hz)', fontsize=30)
ax1.set_ylabel('Wavenumber (rad/m)', fontsize=30)
ax1.tick_params(axis='x', labelsize=25, length=10)
ax1.tick_params(axis='y', labelsize=25, length=10)
#ax1.set_xlim(xmin, xmax)
#ax1.set_ylim(ymin, ymax)

ax1.set_title(f"Spectrum with SVD filter - {date}-{year}-{acqu_numb}\n"
              f"Channel: {ch}, Source: {composante}, Direction: {direction}", fontsize=30)
colorbar = plt.colorbar(c1, ax=ax1)
colorbar.set_label('Spectrum Magnitude', fontsize=25)
colorbar.ax.tick_params(labelsize=20)

#-----------------------------------------------
# Set parameters 
threshold = 0.2
precision_k = [0.09,0.025] # precision over k values (when searching for a local maximum)
prec = precision_k[flexure_wave]
# 0: horizontal_wave 1: flexure_wave
semirange_freq = 5 # semi-range of frequencies for horizontal waves
#-----------------------------------------------

plt.show()
if flexure_wave:
   
    points = plt.ginput(10, timeout=-1)
    # Extract x and y coordinates of the selected points
    f_mode, k_mode = zip(*points)
    f_mode = np.array(f_mode)
    k_mode = np.array(k_mode)
    # Create a spline interpolation of the points
    spline = make_interp_spline(f_mode, k_mode)
    # Generate 50 equally spaced points between the minimum and maximum of the original x values
    f_mode_new = np.linspace(f_mode.min(), f_mode.max(), 50)
    k_mode_new = spline(f_mode_new)
    k_p = k_mode_new + prec
    k_m = k_mode_new - prec
    # Plot the original points


    
    # Extract kmax values
    kmax_values = []
    for i, f in enumerate(f_mode_new): # loop over interpolated frequencies 
        f_idx = np.argmin(np.abs(F - f))  # Find the closest frequency index in F
        k_range = K_uwp[(K_uwp >= k_m[i]) & (K_uwp <= k_p[i])]  # k values within bounds
        k_indices = np.where((K_uwp >= k_m[i]) & (K_uwp <= k_p[i]))[0]  # indices of these k values
        FK_values = FK_uwp[k_indices, f_idx]  # corresponding FK values within the bounds
        # Filter out values below the threshold
        valid_indices = np.where(FK_values >= threshold)[0]
        if valid_indices.size > 0:  # Check if there are any valid indices left
            k_range_valid = k_range[valid_indices]
            FK_values_valid = FK_values[valid_indices]
            kmax_idx = np.argmax(FK_values_valid)  # index of maximum FK value within valid bounds
            kmax = k_range_valid[kmax_idx]  # corresponding k value
            kmax_values.append(kmax)
        else:
            kmax_values.append(np.nan)  # Append NaN if no values are above the threshold

    kmax_values = np.array(kmax_values)
    [f_mode, k_mode] = [f_mode_new, kmax_values]

    plt.plot(f_mode_new, kmax_values, linestyle='--', color='k', label='kmax values')
    plt.xlabel('f_mode_new')
    plt.ylabel('kmax')
    plt.legend()
    plt.show()
    plt.scatter(f_mode, k_mode, color='b', label='Original points')
    plt.plot(f_mode_new, k_p, linestyle='--', color='g', label='sup bound line')
    plt.plot(f_mode_new, k_m, linestyle='--', color='g', label='inf bound line')
    plt.xlabel('f_mode')
    plt.ylabel('k_mode')
    plt.legend()
    plt.show()

    
elif horizontal_wave:
    # ------------------------------------------------------------------------------------------
    points = [(0, 0), (0, )]
    points[1] = plt.ginput(1, timeout=-1)[0]
    # Extract x and y coordinates of the selected points
    f_points, k_points = zip(*points)

    # Assign start and end frequencies
    f_points_sorted = sorted(f_points)
    f_start = int(f_points_sorted[0])
    f_end = int(f_points_sorted[1])

    # Assign start and end k_values
    k_points_sorted = sorted(k_points)
    k_start1 = k_points_sorted[0]
    k_end1 = k_points_sorted[1]

    angle = math.atan((k_end1 - k_start1) / (f_end - f_start))

    freq_values_liste = []
    k_values_liste = []

    # Finding the maximum spectrum magnitude within the range
    for frequency in range(f_start, f_end, 1):
        # Find the index corresponding to the current frequency
        frequency_index = np.argmin(np.abs(f - frequency))

        # -----------redefining the range---------------------

        k_start = (k_start1 - prec) + (frequency - f_start) * math.tan(angle)
        k_end = (k_start1 + prec) + (frequency - f_start) * math.tan(angle)
        # --------------------------------------------

        # Find the indices corresponding to the selected range
        k_indices_range = np.where((k >= k_start) & (k <= k_end))[0]

        frequency_index_range = [frequency_index - semirange_freq, frequency_index + semirange_freq]

        # Find the maximum spectrum magnitude values within the range
        max_spectrum_magnitudes = np.sort(FK_normalized[k_indices_range, frequency_index])[::-1]
        max_spectrum_indices = np.argsort(FK_normalized[k_indices_range, frequency_index])[::-1]

        # Plot the points with the highest spectrum magnitude for the current frequency
        max_spectrum_magnitude = max_spectrum_magnitudes[0]
        max_spectrum_index = max_spectrum_indices[0]
        k_value = k_uwp[k_indices_range[max_spectrum_index]]
        ax1.scatter(frequency, k_start, color='green', marker='1')
        ax1.scatter(frequency, k_end, color='green', marker='1')

        # Plot the point if it falls within the specified limits
        if max_spectrum_magnitude > threshold:  # doesn't take in account the lower spectrum magnitude within the range
            ax1.scatter(frequency, k_value, color='red')
            freq_values_liste.append(frequency)
            k_values_liste.append(k_value)
    plt.show()

    # Plot the linear regression line
    coefficients = np.polyfit(freq_values_liste, k_values_liste, 1)
    poly_function = np.poly1d(coefficients)

    plt.plot(freq_values_liste, poly_function(freq_values_liste), color='blue', linestyle='--', label='Linear Regression')

    # Redraw the graph with the new limits, added scatter points, and linear regression line
    plt.show()

    wave_speed = 2 * np.pi * (1 / coefficients[0])
    print("Speed of the wave:", wave_speed, "m/s")
    

save = 0
if save:
    file2save = path2data  + '/' + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir.pkl'
    with open(file2save, 'wb') as f:
        pickle.dump(wave_speed, f)
    
    file2save = path2data  + '/' + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir.pkl'
    with open(file2save, 'wb') as f:
        pickle.dump(wave_speed, f)

    file2save = path2data  + '/' + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
    with open(file2save, 'wb') as f:
        pickle.dump([f_mode, k_mode], f)




