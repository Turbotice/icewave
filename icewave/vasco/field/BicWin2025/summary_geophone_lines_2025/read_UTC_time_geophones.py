#%% Import modules 
import numpy as np
from datetime import datetime
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy import read
import os
import pickle
import time
import csv 
import sys
sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/field/BicWin2025/')
from BicWin2025.python_functions.sismo_analysis import *


#%% define function to read start and end time

def read_tstart_tend(year='2025', date = '0210', acqu_numb = '0001', ordi = 'dell_vasco'):

    if ordi == 'babasse':
        path2data = os.path.join('E:/Data/',date,'Geophones/')
    elif ordi == 'dell_vasco':
        path2data = f'B:/Data/{date}/Geophones/'
    #    path2data = f'D:/BicWin2025/Data/{date}/Geophones/'
    elif ordi == 'adour':
        path2data = f'/media/turbots/Backup25/Data/{date}/Geophones/'


    fig_folder = path2data + 'Figures/'  # folder where figures are saved 
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)

    if ordi=='babasse':    
        geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
    elif ordi=='dell_vasco':
        geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
    elif ordi=='adour':
        geophones_table_path = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'

    channel = 0  # 0 for E, 1 for N, 2 for Z. 

    #files need to be organised as: data/0210/Geophones/0001/minised files

    geophones_spacing = 3 # space between geophones, in meters 
    signal_length = 1 # duration in seconds 
    channel_dic = {
        1: "N",
        2: "Z",
        0: "E",}
    ch = channel_dic[channel]


    ######
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

    return tstart, tend
# %%

#read_tstart_tend()