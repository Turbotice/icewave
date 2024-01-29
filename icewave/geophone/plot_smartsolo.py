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


path2data = '/Users/moreaul/Desktop/Geophones/Grenoble_20240118_dehors'
acqu_numb = '0001'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

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


num_traces_to_plot = 16
fig, ax = plt.subplots()
# Calculate datetime values only once outside the loop
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector]

for k in range(1, num_traces_to_plot + 1):
    stream_index = 1 + k  # Adjust index to match Python's 0-based indexing
    current_stream = seismic_data_streams[stream_index]

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


# Define the desired time range



t1 = UTCDateTime("2024-01-18T10:35:14") 
t2 = t1 + 1


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

# Plot for k = 2, 5, 8, ..., 48
selected_indices = range(2, 49, 3)

fig, ax = plt.subplots()

for k in selected_indices:
    ax.plot(time_vector, seismic_matrix[k, :] / max(np.abs(seismic_matrix[k, :])) + k, label=f"Stream {k}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams')
#ax.legend()
plt.show()







