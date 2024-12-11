# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:20:24 2024

@author: sebas

This script can be decomposed into 2 main Parts : 
    Part 1: extraction of initial times for different sources performed on field 
    Part 2: SVD decomposition of each geophone channel (E,N,Z) in order to compute phase velocity
    of modes QS0 and SH0, and infer Young modulus E and Poisson coefficient nu. 
    
    Extraction of coordinates (f,k) from FK plot of channel Z (flexural waves) is also performed in order to infer 
    ice thickness and ice density in an other Python script : 'Inv_disp_curves_bidir.py'

"""

#%% Import modules 
import numpy as np
import math
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
from scipy.signal import ShortTimeFFT, windows, iirnotch, filtfilt
from scipy.interpolate import make_interp_spline

import time

#%% Set parameters 
year = '2024'
date = '0914' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 

path2data = os.path.join('F:/Amundsen_RA_2024/Data/',year,date,'Geophones/')

fig_folder = path2data + 'Figures/'  # folder where figures are saved 
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
channel = 2  # 0 for E, 1 for N, 2 for Z. 

#files need to be organised as: data/0210/Geophones/0001/minised files
geophones_spacing = 8 # space between geophones, in meter 
signal_length = 1.0 # duration in seconds 
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}
ch = channel_dic[channel]

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
img_quality = 200 # dpi to save images 

plt.rcParams.update({
    "text.usetex": True})

fig_folder = path2data + acqu_numb + '/' + 'Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    


#%% Function definition 

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


#----------------------------------------------------------------------------------------------
def fn_svd ( signals , fs , xs , rang , name , issaving , *varargin ) : #varagin ??
    """ Computes spatio-temporal Fourier Transform using SVD decomposition.
    Arguments : 
        - signals : matrix #dim 1 : nb_recepters, #dim 2 : time, #dim 3 : nb_emitters (sources)
        - fs : sampling frequency in time 
        - xs : distance between geophones (inverse of spatial sampling)
        - rang : row used to perform decomposition using SVD method (number of singulavr values used to perform decomposition)
        - name : title of figure showing singular values 
        - issaving : boolean to choose to save figue showing singular values or not 
        - varargin : optionnal argument
        
    Returns : 
        
        - f : array of frequencies, scaled
        - k : array of wavenumber, scaled 
        - projections_sum : amplitude of FFT of the signal in time and space, using SVD method 
        #dim 1 : k
        #dim 2 : f """
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

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

###################################################################
#%% -------------------- Loading geophones Data -------------------
###################################################################


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
time_array = first_trace.times() # experimental time array 
time_vector = convert_to_utc_times(first_trace, first_trace.times())
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
data_vector = first_trace.data
start_time_utc = UTCDateTime(first_trace.stats.starttime) 
fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 

#%% ---------------- Plot a single channel for a given geophone ----------------------

idx_geophone = 2
chan = 0
selected_stream = seismic_data_streams[(idx_geophone - 1)*3 + chan]
print(selected_stream[0])


fig, ax = plt.subplots(figsize = fig_size)
ax.plot(datetime_values,selected_stream[0].data / max(np.abs(selected_stream[0].data)))


#%% Check FFT of several geophones 

channel = 2
N = len(first_trace.times())
selected_indices = range(channel,len(seismic_data_streams),3)

# build matrix of raw signals 
raw_matrix = np.zeros((len(selected_indices),N))

for k,idx in enumerate(selected_indices):
    current_stream = seismic_data_streams[idx]
    print(current_stream[0])
    raw_matrix[k,:] = current_stream[0].data
    
# Compute FFT spectrum of all signals 
FFT_matrix = np.fft.fft(raw_matrix, norm = 'forward')
FFT_matrix = FFT_matrix[:,:N//2]
freq = np.fft.fftfreq(N,1/fs)
freq = freq[:N//2]
    
#%% plot FFT spectrum of all signals

fig, ax = plt.subplots()
for k in range(1,5):
    ax.loglog(freq,abs(FFT_matrix[k,:]))
ax.set_xlabel(r'$f \: \rm (Hz)$')
ax.set_ylabel(r'$\hat{V}(f)$')
ax.grid()


#%% Filter Amundsen noise for all geophones
                      
# fc = [8.333 , 16.666, 25 , 2*16.666 , 41.66, 3*16.666, 58.33, 4*16.666, 75, 5*16.66, 
#       91.66, 6*16.66, 7*16.66, 8*16.66] 

# Load frequencies to cut 
file_Amundsen_frequencies = path2data + 'Frequencies_Amundsen_' + year + '_' + date + '_channel_'+ channel_dic[channel]  + '.pkl'
if os.path.isfile(file_Amundsen_frequencies):
    print('File of Amundsen frequencies already exists, loading..')
    with open(file_Amundsen_frequencies,'rb') as pfile:
        fc = pickle.load(pfile)
else :
    print('File of Amundsen frequencies does not exist, saving..')
    with open(file_Amundsen_frequencies,'wb') as pfile:
        pickle.dump(fc,pfile)

# use band-cut filter 
Q = 30.0  # Quality factor
filtered_matrix = raw_matrix
for f0 in fc:
    # Design notch filter
    b, a = iirnotch(f0, Q, fs)
    # Filter frequencies
    filtered_matrix = filtfilt(b,a,filtered_matrix)

#%% Plot filtered and raw signals

fig, ax = plt.subplots(figsize = fig_size)
for k in range(np.shape(filtered_matrix)[0]):
    current_stream = filtered_matrix[k,:]
    print(k)
    # norm_data = (current_stream - np.min(current_stream)) / (np.max(current_stream) - np.min(current_stream))
    norm_data = current_stream / np.max(np.abs(current_stream))
    ax.plot(datetime_values,norm_data + k * geophones_spacing)

ax.set_xlabel(r'$t \: \rm (s)$')                      
ax.set_ylabel(r'$x \: \rm (m)$')                          
ax.set_title(r'Filtered Data channel ' + channel_dic[channel])


fig, ax = plt.subplots(figsize = fig_size)
for k in range(np.shape(raw_matrix)[0]):
    current_stream = raw_matrix[k,:]
    print(k)
    # norm_data = (current_stream - np.min(current_stream)) / (np.max(current_stream) - np.min(current_stream))
    norm_data = current_stream / np.max(np.abs(current_stream))
    ax.plot(datetime_values,norm_data + k * geophones_spacing)

ax.set_xlabel(r'$t \: \rm (s)$')                      
ax.set_ylabel(r'$x \: \rm (m)$')                          
ax.set_title(r'Raw Data channel ' + channel_dic[channel])


#################################################################################
#%%----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
#################################################################################
# Select signal segments manually for each track and each source 
# Store selected times in a dictionnary and save it as a .pkl file 


# open .pkl file and load dictionnary 
filename = 't1_to_time_' + date + '_' + year + '_first_hit.pkl'
base = path2data
file2save = base + filename
if os.path.isfile(file2save):
    print('Time dictionnary already exists')
    with open(file2save,'rb') as f:
        time_dict = pickle.load(f)
else:
    print('No time dictionnary saved yet')
    time_dict = {}

composante = 'N' #Z , E or N -> direction de la source


key = 'd' + date + 'a' + acqu_numb + 'tS' + '001' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:04:28.94")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '002' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:05:39.57")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '003' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:06:57.76")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '004' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:08:06.03")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '005' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:09:41.43")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '006' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:11:57.72")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '007' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:13:07.85")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '008' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:14:11.69")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '009' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:15:22.65")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '010' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:16:32.45")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '011' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:18:50.00")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '012' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:20:08.92")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '013' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:21:32.80")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '014' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:24:17.26")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '015' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:25:50.77")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '016' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:29:27.40")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '017' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:30:28.51")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '018' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:31:33.68")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '019' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:32:37.65")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '020' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:33:50.10")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '021' + composante 
time_dict[key] = UTCDateTime("2024-09-18T17:35:10.64")




#%% Save t0 dictionnary in pickle file 

with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)

print('Time dictionnary saved')








