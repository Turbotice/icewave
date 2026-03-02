

import numpy as np
import math
from datetime import datetime
from datetime import timedelta
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
import csv 
from scipy.interpolate import interp1d

import scipy
from scipy.signal import find_peaks

#########################################
# fonctions geophones (Ludovic)

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





def load_geophone_data(date,acqu_numb,path2data,geophones_table_path,channel):
    #table into a dictionary
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
    #signal_length = 1
#    channel_dic = {
#        1: "N",
#        2: "Z",
#        0: "E",}
#    ch = channel_dic[channel]

    ####################### 
    first_stream = seismic_data_streams[0]
    first_trace = first_stream[0]
    time_array = first_trace.times() # experimental time array 
    time_vector = convert_to_utc_times(first_trace, first_trace.times())
    datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
    data_vector = first_trace.data
    start_time_utc = UTCDateTime(first_trace.stats.starttime) 
    fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 

    return geophones_dict, seismic_data_streams, first_stream, first_trace, tstart, tend, fs, selected_indices, channel_dic, ch, time_array, time_vector, datetime_values, data_vector, start_time_utc#, signal_length, 







##########################################
# fonctions Vasco

def find_9shocks(signal,typical_signal,params={'width':2 , 'distance':100 , 'threshold':1e-2 , 'height':100}):# adjust height for find_peaks
    
    """
    function to find the 3*3 = 9 shocks in the desired signal
    inputs : 
    - signal (not normalized)
    - typical signal with which we want to correlate the signal to find the times of the 9 shocks (not normalized)
    Remark : returns indices (not times in UTC) that index the array
    """
    typical_signal = typical_signal/np.max(np.abs(typical_signal)) # typical signal is normalized
    typical_signal_withzeros = np.zeros(len(signal))
    typical_signal_withzeros[:len(typical_signal)] = typical_signal

    ycorr = scipy.signal.correlate(signal,typical_signal_withzeros)
    lags = scipy.signal.correlation_lags(len(signal),len(typical_signal_withzeros))

    pp = find_peaks(ycorr)
    #ppp = find_peaks(ycorr[pp[0]],width=2,distance=100,threshold=1e-2,height=80)
    ppp = find_peaks(ycorr[pp[0]],width=params['width'],distance=params['distance'],threshold=params['threshold'],height=params['height'])
    peak_heights = ppp[1]['peak_heights']
    argsort = np.argsort(peak_heights)
    nine_highest_peaks_idx = argsort[:9]
    plt.figure()
    plt.plot(lags,ycorr)
    plt.show()
    plt.plot(lags[pp[0]],ycorr[pp[0]])
    plt.vlines(lags[pp[0]][ppp[0]][nine_highest_peaks_idx],ymin=0,ymax=np.max(ycorr),colors='red',linestyles='--')
    print(lags[pp[0]][ppp[0]][nine_highest_peaks_idx])
    plt.show()
    if len(peak_heights)<9:
        print('didnt find all the shocks')
    return lags[pp[0]][ppp[0]]



def LongLat2deltaEdeltaN(coord1,coord2,latitude_approx,R=6371e3):
    coord1_rad = np.radians(np.array(coord1))
    coord2_rad = np.radians(np.array(coord2))
    latitude_approx_rad = np.radians(latitude_approx)
    
    deltaE = R * np.sin(np.pi/2-latitude_approx_rad) * (coord2_rad[0]-coord1_rad[0])
    deltaN = R * (coord2_rad[1]-coord1_rad[1])
    
    return deltaE,deltaN


def find_LongLat_from_idx_geophone(idx_geophone,path2coordfile = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/Geophones/GPS_coordinates_acq4.txt',path2geophonetable = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'):
    geophones_table = np.loadtxt(path2geophonetable,skiprows=1)
    geo_SN = int(geophones_table[:,1][np.where(idx_geophone==geophones_table)[0][0]])
    geo_SN_str = str(geo_SN)
    geo_SN_last4str = geo_SN_str[-4:]
    geo_coord_table = np.loadtxt(path2coordfile,skiprows=1,dtype=str)
    indrow = np.where(geo_coord_table=='IGU_'+geo_SN_last4str)[0][0]
    latitude = float(geo_coord_table[indrow,1])
    longitude =  float(geo_coord_table[indrow,2])
    return longitude,latitude

def read_data_choosing_geophone_channel(seismic_data_streams, idx_geophone,chan,datetime_indices_to_plot):
    """
    fonction permettant de specifier un geophone et de lire son contenu : en input, donner 
    seismic_data_streams qui est donné par la fonction read_data, l'indice du geophone,
    le channel auquel on s'intéresse et (0 pour E, 1 pour N, 2 pour Z) et la plage d'inces de temps
    datetime_indices_to_plot. Output : signal (array 1d). 
    """
    selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]

    current_stream = selected_streams[chan]
    print(current_stream[0])
    
    #ax[k].plot(np.array(datetime_values)[datetime_indices_to_plot],current_stream[0].data[datetime_indices_to_plot] / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) ),label = channel_dic[chan])        
    # ici il ne faut pas normaliser le signal (plus haut c'est fait seulement pour plotter tous les signaux)
    signal = current_stream[0].data[datetime_indices_to_plot]# / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) )
    return signal


def compute_mean_angle(list_deltaE,list_deltaN):
    '''
    Fonction qui calcule l'angle moyen que l'on considère pour la tomographie "sur une ligne".
    Attention : il faut bien que tous les géophones ainsi que la source soient alignés !!!
    '''
    angles = np.zeros(len(list_deltaE))
    for i in range(len(list_deltaE)):
        angles[i] = np.arctan2(list_deltaN[i],list_deltaE[i])
    mean_angle = np.mean(angles) # on peut le calculer plus précisément...
    return mean_angle

"""
def extract_longi_transv_signal(idx_geophone,direction):
    '''
    prend en argument idx_geophone (le numéro du geophone considéré),
    et direction (qui peut être 'L' pour longi et 'T' pour transverse)
    '''
    chan = 0 # for E
    signal_E = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
    chan = 1 # for N
    signal_N = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
    if direction == 'L':
        signal_longi = np.cos(mean_angle) * signal_E + np.sin(mean_angle) * signal_N
        return signal_longi
    elif direction == 'T':
        signal_transv = np.sin(mean_angle) * signal_E - np.cos(mean_angle) * signal_N
        return signal_transv
"""


def detect_time_shift(signal1_cut,signal2_cut,freq_acq=1000,plot=False):
    cross_corr = scipy.signal.correlate(signal2_cut,signal1_cut)
    lags = scipy.signal.correlation_lags(len(signal2_cut),len(signal1_cut))
    if plot:
        plt.figure()
        plt.plot(signal1_cut)
        plt.plot(signal2_cut)
        plt.plot(lags,cross_corr)
        plt.show()

    idx_max = np.argmax(cross_corr)

    # detection subpixellaire du max :
    x1 = lags[idx_max-1]
    y1 = cross_corr[idx_max-1]
    x2 = lags[idx_max]
    y2 = cross_corr[idx_max]
    x3 = lags[idx_max+1]
    y3 = cross_corr[idx_max+1]

    A = np.array([[x1**2,x1,1],[x2**2,x2,1],[x3**2,x3,1]])
    B = np.array([y1,y2,y3])
    Coefficients = np.linalg.solve(A,B)

    x_max_subpix = -Coefficients[1]/(2*Coefficients[0])
    delta_t_subpix = x_max_subpix/freq_acq
    return delta_t_subpix

