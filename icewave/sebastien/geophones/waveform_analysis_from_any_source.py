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

get_ipython().run_line_magic('matplotlib', 'qt')
plt.close('all')

# ## Input

# In[50]:


year = '2024'
date = '0211' #date format, 'mmdd'
acqu_numb = '0003' #acquisition number 

direction = 2 # 1 ou 2 
composante = 'N' #Z , E or N -> direction de la source
channel = 1  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files


path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'



# ## Fonctions

# In[51]:


#----------------------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------------
def convert_to_utc_times(trace, time_vector):
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]
#----------------------------------------------------------------------------------------------
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
            print(
                f"Warning: No entry found for {last_five_digits} in the geophones table.")

    # Sort traces based on the new station codes
    sorted_stream = sorted(stream, key=lambda trace: trace.stats.station)

    return sorted_stream


#----------------------------------------------------------------------------------------------
def onclick(event):
    if len(points) < 2:
        points.append((event.xdata, event.ydata))
        if len(points) == 2:
            # Finalize the line after the second point is clicked
            line.set_data([points[0][0], points[1][0]], [points[0][1], points[1][1]])
            #fig.canvas.draw()
            fig.canvas.mpl_disconnect(cid_click)
            fig.canvas.mpl_disconnect(cid_move)
            # After selecting the points, check if there are enough points before accessing them
            if len(points) == 2:
                f_points, k_points = zip(*points)

#----------------------------------------------------------------------------------------------
def onmove(event):
    if len(points) == 1:
        line.set_data([points[0][0], event.xdata], [points[0][1], event.ydata])
        #line1.set_data([points[0][0] , event.xdata + pres], [points[0][1], event.ydata + pres ])
        #line2.set_data([points[0][0] , event.xdata + pres], [points[0][1], event.ydata + pres ])
        fig.canvas.draw()


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


# If "Unknown format for file 0005\._01.0005.2024.02.11.20.06.55.000.E.miniseed" appears
# Create a NEW 000n file only with the miniseed files. 
# assign a string to the channel values 
signal_length = 1
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
var_path = "variable_to_dict.pkl"
with open(var_path, 'rb') as m:
    variable_to_dict = pickle.load(m)
if variable_to_dict.get(date+'_'+year) is None:
    messagebox.showerror("Error", "Date not found")
    raise RuntimeError("Date not found")
    
# get the pkl file assigned to the date and year value
pkl_path =variable_to_dict.get(date+'_'+year) + ".pkl"
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
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)

ax1.set_title(f"Spectrum with SVD filter - {date}-{year}-{acqu_numb}\n"
              f"Channel: {ch}, Source: {composante}, Direction: {direction}", fontsize=30)
colorbar = plt.colorbar(c1, ax=ax1)
colorbar.set_label('Spectrum Magnitude', fontsize=25)
colorbar.ax.tick_params(labelsize=20)

plt.show()
#------------------------------------------------------------------------------------------
points = [(0,0),(0,)]
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

angle = math.atan((k_end1-k_start1)/(f_end-f_start))

freq_values_liste = []
k_values_liste = []

# Finding the maximum spectru magnitude within the range
for frequency in range(f_start, f_end, 1):
    # Find the index corresponding to the current frequency
    frequency_index = np.argmin(np.abs(f - frequency))
    
    #-----------redefining the range---------------------
    pres = 0.09 # permet d'ajuster la largeur de la zone delimité
    k_start = (k_start1 - pres) + (frequency - f_start)* math.tan(angle)
    k_end = (k_start1 + pres) + (frequency - f_start)* math.tan(angle)
    #--------------------------------------------
    
    # Find the indices corresponding to the selected range
    k_indices_range = np.where((k >= k_start) & (k <= k_end))[0]
    
    frequency_index_range = [frequency_index - 5 , frequency_index +5 ]
    
    # Find the maximum spectrum magnitude values within the range
    max_spectrum_magnitudes = np.sort(FK_normalized[k_indices_range, frequency_index])[::-1]
    max_spectrum_indices = np.argsort(FK_normalized[k_indices_range, frequency_index])[::-1]
    
    # Plot the points with the highest spectrum magnitude for the current frequency
    max_spectrum_magnitude = max_spectrum_magnitudes[0]
    max_spectrum_index = max_spectrum_indices[0]
    k_value = k_uwp[k_indices_range[max_spectrum_index]]
    ax1.scatter(frequency, k_start, color='green',marker = '1')
    ax1.scatter(frequency, k_end, color='green', marker = '1')
    
    # Plot the point if it falls within the specified limits
    if max_spectrum_magnitude > 0.6:  # doesn't take in acount the lower spectrum magnitude within the range
        ax1.scatter(frequency, k_value, color='red')
                    #, label=f'Highest spectrum magnitude ({frequency} Hz)')
        freq_values_liste.append(frequency)
        k_values_liste.append(k_value)
plt.show()

# Plot the linear regression line
coefficients = np.polyfit(freq_values_liste, k_values_liste, 1)
poly_function = np.poly1d(coefficients)

plt.plot(freq_values_liste, poly_function(freq_values_liste), color='blue', linestyle='--', label='Linear Regression')

# Redraw the graph with the new limits, added scatter points, and linear regression line
plt.show()

wave_speed = 6.28 * (1 / coefficients[0])
print("Speed of the wave:", wave_speed, "rad/s")



