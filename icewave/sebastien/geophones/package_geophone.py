# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 10:55:55 2025

@author: sebas

Gather all useful functions for geophones data 

"""

import numpy as np
from datetime import datetime
from obspy.core import UTCDateTime
from obspy import read
import os
import pickle

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

#---------------------------------------------------------------------------------------------

def build_data_streams(path2miniseeds,geophones_table_path):
    """ Build seismic data streams matrix, datetime values as well as sampling frequency. 
    The data matrix is ordered as an increasing order of geophone index, with each of the three 
    channels following for a single geophone : [E1,N1,Z1,E2,N2,Z2,...]
    Inputs : - path2miniseeds, path to miniseeds files
             - geophones_table_path, path to geophone correspondance table
    Outputs : - seismic_data_streams, matrix of all geophones streams
              - datetime_values, array of date time 
              - fs, float, sampling frequency
        """
    
    
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
    seismic_data_streams, miniseed_files = read_data(path2miniseeds)

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
    time_vector = convert_to_utc_times(first_trace, first_trace.times())
    datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
    fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 

    
    return seismic_data_streams,datetime_values,fs
    