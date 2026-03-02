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
from datetime import datetime
import numpy as np


date = '0221'
year = '2024'


path2data = f'/YOUR_DATA_PATH_HERE/{year}_BICWIN/{date}/Geophones'
acqu_numb = '0001'
geophones_table_path = '/icewave/geophone/geophones_table'



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

channel = 2 #0 for E, 1 for N, 2 for Z.  
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



