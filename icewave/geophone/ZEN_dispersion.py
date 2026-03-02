# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 09:52:37 2026

@author: sebas

Compute seismic waves dispersion relations from geophones seismic measurements. 
Arrays of time sources must have been already generated using the script "extract_times.py"

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
import math

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

fig_folder = path2data + acqu_numb + '/Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% FUNCTION SECTION 

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]
    
###################################################################
#%% -------------------- Loading geophones Data -------------------
###################################################################

seismic_data_streams,datetime_values,fs = geopack.build_data_streams(path2data + acqu_numb +'/',geophones_table_path)

##############################################################################################
#%%----------------------- PLOTTING SELECTED CHANNEL FOR ALL GEOPHONES ------------------ 
##############################################################################################
 
channel = 0 #0 for E, 1 for N, 2 for Z. 

selected_indices = range(channel, len(seismic_data_streams),3) # change step of arange if needed!! 
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
    
#%% Load time dictionnary 

signal_length = 1 # duration in seconds

# load data of intial times 
composante = 'N'
channel = 1
direction = 2 # 1 ou 2

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
# t1_values = [t2,t3]

# Create a matrix to store the seismic data
selected_indices = np.arange(channel, len(seismic_data_streams),3) # change step of arange !! 
ta = t2 + signal_length 
num_samples = int(signal_length * fs) # Number of points generated in time_vector 
time_vector = np.linspace(t2.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

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

#%%  Calling fn_SVD

#####################################
##### IMPORTANT PARAMETERS ##########
#####################################

rang = [0,1,2]
signals = np.transpose(seismic_matrix, (0, 2, 1))

f, k, FK = geopack.fn_svd(signals, fs, geophones_spacing , rang ,'ExampleName', 0, 'threshold',-60) #plot valeurs  singuliere/freq
F, K = np.meshgrid(f, k)

# Normalize FK based on direction
if direction == 1:
    FK_normalized = FK / np.max(FK)
elif direction == 2:
    FK_normalized = np.flipud(FK / np.max(FK))

# Unwrapping Spectrum
nb_stacking = 2 # number of times we want to stack the FK plot horizontally
idx_stacking = 1
FK_uwp = np.vstack((FK_normalized, FK_normalized))
while idx_stacking < nb_stacking :
    FK_uwp = np.vstack((FK_uwp, FK_normalized))
    idx_stacking += 1
    
k_uwp= np.linspace(0,(nb_stacking + 1)*max(k),FK_uwp.shape[0])
F, K_uwp = np.meshgrid(f, k_uwp)

################################################
#%%---------- Plotting FK ---------------------
################################################

vmin = 0     # Minimum value for colormap
vmax = 1     # Maximum value for colormap (adjust as needed)
ymin = 0  # Minimum frequency
ymax = 400  # Maximum frequency
xmin = 0  # Minimum wavenumber
xmax = 4  # Maximum wavenumber


#-----------------------------------------------
# Set parameters 
threshold = 0.2 # minimum relative amplitude to detect a local maximum in the dispersion relation FK plot 
precision_k = [0.03,0.025] # precision over k values (when searching for a local maximum)
prec = precision_k[flexure_wave]
# 0: horizontal_wave 1: flexure_wave
semirange_freq = 5 # semi-range of frequencies for horizontal waves
nb_points = 12 # points to select on graph FK flexural
#-----------------------------------------------

fig, ax1 = plt.subplots(1, 1, figsize=fig_size)
plt.show()
if flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
    #ax1.set_ylim([-1.5, 1.5])
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')

    ax1.set_ylim([0, 200])
    print(f"Select {nb_points} points on the graph")
    points = plt.ginput(nb_points, timeout=-1)

    # Extract x and y coordinates of the selected points
    k_mode, f_mode = zip(*points)
    f_mode = np.array(f_mode)
    k_mode = np.array(k_mode)
    # Create a spline interpolation of the points
    spline = make_interp_spline(k_mode, f_mode)
    # Generate 50 equally spaced points between the minimum and maximum of the original x values
    k_mode_new = np.linspace(k_mode.min(), k_mode.max(), 50)
    f_mode_new = spline(k_mode_new)
    k_p = k_mode_new + prec
    k_m = k_mode_new - prec
    # Plot the original points

    # Extract values of k for which FK is maximum at a given frequency 
    kmax_values = []
    for i, freq in enumerate(f_mode_new): # loop over interpolated frequencies 
        f_idx = np.argmin(np.abs(F - freq))  # Find the closest frequency index in F
        k_range = K_uwp[(K_uwp >= k_m[i]) & (K_uwp <= k_p[i])]  # k values within bounds
        k_indices = np.where((K_uwp >= k_m[i]) & (K_uwp <= k_p[i]))[0]  # indices of these k values
        FK_values = FK_uwp[k_indices, f_idx]  # corresponding FK values within the bounds
        # Filter out values below the amplitude threshold
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

    # plt.plot(kmax_values, f_mode_new, linestyle='--', color='k', label='kmax values')
    # plt.ylabel('f_mode_new')
    # plt.xlabel('kmax')
    # plt.legend()
    # plt.show()
    plt.scatter(k_mode, f_mode, color='w', label='Selected points')
    plt.plot( k_p, f_mode_new, linestyle='--', color='g', label='sup bound line')
    plt.plot( k_m, f_mode_new, linestyle='--', color='g', label='inf bound line')
    plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
    plt.ylabel(r'$f \; \mathrm{(Hz)}$')
    plt.legend()
    plt.show()

#--------------------------------------------------------------------------------
elif horizontal_wave:

    # Computes extents for imshow
    x = extents(k)
    y = extents(f)

    c1 = ax1.imshow(FK_normalized.T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

    #ax1.set_ylim([-1.5, 1.5])
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    # ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}|_{max}}(f,k)$')
    ax1.tick_params(axis='both', which='major', pad=7)
    ax1.set_ylim([ymin, ymax])
    # ------------------------------------------------------------------------------------------
    print('Select a single point on the graph')
    points = [(0, 0), (0, )]
    points[1] = plt.ginput(1, timeout=-1)[0]
    # Extract x and y coordinates of the selected points
    k_points, f_points = zip(*points)

    # Assign start and end frequencies
    f_points_sorted = sorted(f_points)
    f_start = int(f_points_sorted[0])
    f_end = int(f_points_sorted[1])

    # # Assign start and end k_values
    k_points_sorted = sorted(k_points)
    k_start1 = k_points_sorted[0]
    k_end1 = k_points_sorted[1]

    # angle = math.atan((k_end1 - k_start1) / (f_end - f_start))
    angle = math.atan((f_end - f_start)/(k_end1 - k_start1))

    freq_values_list = []
    k_values_list = []

    # Finding the maximum spectrum magnitude within the range
    for frequency in range(f_start, f_end, 1):
        # Find the index corresponding to the current frequency
        frequency_index = np.argmin(np.abs(f - frequency))

        # -----------redefining the range---------------------

        k_start = (k_start1 - prec) + (frequency - f_start) / math.tan(angle)
        k_end = (k_start1 + prec) + (frequency - f_start) / math.tan(angle)
        # --------------------------------------------

        # Find the indices corresponding to the selected range
        k_indices_range = np.where((k >= k_start) & (k <= k_end))[0]

        frequency_index_range = [frequency_index - semirange_freq, frequency_index + semirange_freq]

        # Find the maximum spectrum magnitude values within the range
        max_spectrum_magnitudes = np.sort(FK_normalized[k_indices_range, frequency_index])[::-1]
        max_spectrum_indices = np.argsort(FK_normalized[k_indices_range, frequency_index])[::-1]

        # Plot the points that corresponds to window limits
        max_spectrum_magnitude = max_spectrum_magnitudes[0]
        max_spectrum_index = max_spectrum_indices[0]
        k_value = k[k_indices_range[max_spectrum_index]]
        ax1.scatter( k_start, frequency, color='green', marker='1')
        ax1.scatter( k_end, frequency, color='green', marker='1')

        # Plot the point if it falls within the specified limits
        if max_spectrum_magnitude > threshold:  # doesn't take in account the lower spectrum magnitude within the range
            ax1.scatter(k_value, frequency, color='red')
            freq_values_list.append(frequency)
            k_values_list.append(k_value)
    plt.show()

    # Plot the linear regression line
    
    # fit by a 1D polynome
    deg = 1
    p,V = np.polyfit(k_values_list,freq_values_list,deg,cov = True )
    poly_function = np.poly1d(p)
    C = np.zeros((2,1))
    C[0] = p[0]*2*np.pi
    C[1] = np.sqrt(2*np.pi*np.diag(V)[0]) # standard deviation
    print("Phase velocity :" ,C[0], " m/s")

    plt.plot(k_values_list, poly_function(k_values_list), color='blue', linestyle='--', label='Linear Regression')

    # Redraw the graph with the new limits, added scatter points, and linear regression line
    plt.show()
    
#%% Plot and save a clean graph

fig, ax1 = plt.subplots(1, 1, figsize=fig_size)

if horizontal_wave:
    # Computes extents for imshow
    x = extents(k)
    y = extents(f)

    c1 = ax1.imshow(FK_normalized.T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)
    plt.plot(k_values_list, poly_function(k_values_list), color='white', linestyle='--', linewidth = 3, label='Linear Regression')
    #ax1.set_ylim([-1.5, 1.5])
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    # ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label= r'$|\hat{s}|(f,k)$')
    ax1.tick_params(axis='both', which='major', pad=7)
    ax1.set_ylim([ymin, ymax])
    
    if channel == 0:
        wave_type = 'SH0'
    elif channel == 1:
        wave_type = 'QS0'
        
    figname = fig_folder + 'Horizontal_FKplot_' + wave_type + '_acqu_' + acqu_numb + '_dir'  + str(direction) + '_sig_length_' + str(signal_length).replace('.','p')

    plt.savefig(figname + '.pdf', dpi = 300, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = 300, bbox_inches = 'tight')
    
    
elif flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
    plt.plot(k_mode, f_mode, color='white', linestyle='-', linewidth = 3, label='Selected points')
    #ax1.set_ylim([-1.5, 1.5])
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    # ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label= r'$|\hat{s}|(f,k)$')
    ax1.set_ylim([ymin, ymax])
    
    wave_type = 'QS'
    figname = fig_folder + 'Flexural_FKplot_' + wave_type + '_acqu_' + acqu_numb + '_dir'  + str(direction)

    plt.savefig(figname + '.pdf', dpi = 300, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = 300, bbox_inches = 'tight')
    
#%% Save phase velocity data in a dictionnary  

if horizontal_wave :
    pkl_file = path2data + 'Phase_velocity_dictionnary_acqu_' + acqu_numb + '_sig_length_' + str(signal_length).replace('.','p') + '.pkl'
    if os.path.isfile(pkl_file):
        print('Phase velocity dictionnary already exists')
        with open(pkl_file,'rb') as pfile:
            s = pickle.load(pfile)
    else:
        print('No phase velocity dictionnary saved yet')
        s = {}
        s['dir1'] = {}
        s['dir2'] = {}
        
    keys_phase_velocity = {'N':'C_longi','E':'C_shear'}
    
    key_dir = 'dir' + str(direction)
    key = keys_phase_velocity[composante]
    key_uncertainty = 'u' + key
    s[key_dir][key] = C[0] # value 
    s[key_dir][key_uncertainty] = C[1] # standard deviation 
        
    with open (pkl_file,'wb') as pfile:
        pickle.dump(s,pfile)
    
#%% Save (f,k) points associated to QS mode 

if flexure_wave :
    file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' + str(direction) +'.pkl'
    with open(file2save, 'wb') as file:
         pickle.dump([f_mode, k_mode], file)
         
#%% Load phase velocity data
pkl_file = path2data + 'Phase_velocity_dictionnary_acqu_' + acqu_numb + '_sig_length_' + str(signal_length).replace('.','p') + '.pkl'
with open (pkl_file,'rb') as pfile:
    s = pickle.load(pfile)    


#%% Compute horizontal modes phase velocity in both directions

wave_speed = 0.5*(s['dir1']['C_longi'] + s['dir2']['C_longi'])
file2save = path2data  +  year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir.pkl'
with open(file2save, 'wb') as file:
    pickle.dump(wave_speed, file)
    
wave_speed = 0.5*(s['dir1']['C_shear'] + s['dir2']['C_shear'])
file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir.pkl'
with open(file2save, 'wb') as file:
    pickle.dump(wave_speed, file)
