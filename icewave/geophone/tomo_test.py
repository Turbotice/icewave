#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 14:32:00 2024

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
import pickle
import h5py





#path2coordinates = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0211/GPS'
path2coordinates = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0211/Geometry'


path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0211/Geophones'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'



pickle_file_path = os.path.join(path2coordinates, 'geometry_2024_0211Y_m_vs_X_m_general.pkl')
with open(pickle_file_path, 'rb') as f:
    saved_variables = pickle.load(f)


for i in range(86):  # assuming you have 86 pairs based on the provided data
    x_data = saved_variables[f'xdata_{i}']
    y_data = saved_variables[f'ydata_{i}']
    plt.scatter(x_data, y_data, label=f'Data {i}')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Scatter Plot of Data')
plt.legend()
plt.show()



#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0005'
channel = 2 #0 for E, 1 for N, 2 for Z. 


#----------------------- FUNCTIONS ------------------ 


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
        self.t_begin = None
        self.t_end = None

        # Initialize rectangle selector
        self.rs = RectangleSelector(ax, self.on_rectangle_select, drawtype='box', useblit=True, button=[1],
                            minspanx=5, minspany=5, spancoords='pixels', interactive=True)


        # Connect the right-click event to the on_right_click method
        self.right_click_cid = ax.figure.canvas.mpl_connect('button_press_event', self.on_right_click)

    def on_right_click(self, event):
        if event.button == 3:  # Check for right-click (button 3)
            self.disconnect()

    def on_rectangle_select(self, eclick, erelease):
        # Extract rectangle coordinates
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Apply zoom to the selected area
        self.ax.set_xlim(min(x1, x2), max(x1, x2))
        self.ax.set_ylim(min(y1, y2), max(y1, y2))
        plt.draw()

        # Check if shift key is held down during the click event
        if eclick.key == 'shift':
            # Convert x-coordinates to time values
            t_begin = np.interp(min(x1, x2), self.ax.get_xlim(), self.time_vector)
            t_end = np.interp(max(x1, x2), self.ax.get_xlim(), self.time_vector)

            # Update t_begin and t_end
            self.t_begin = t_begin
            self.t_end = t_end

    def reset_zoom(self):
        self.ax.set_xlim(self.original_xlim)
        self.ax.set_ylim(self.original_ylim)
        plt.draw()

    def disconnect(self):
        self.ax.figure.canvas.mpl_disconnect(self.right_click_cid)
        plt.show()





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

 
#selected_indices = range(channel, len(seismic_data_streams),3)
selected_indices = range(2, 6,3)

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
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder

print("t_begin:", zoom_handler.t_begin)
print("t_end:", zoom_handler.t_end)
