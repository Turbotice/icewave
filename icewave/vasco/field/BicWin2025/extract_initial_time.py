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
import shutil
import os
import re
#import tkinter as tk
#from tkinter import messagebox
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
#%%
%matplotlib qt
#plt.rcParams['text.usetex'] = False

#%% Set parameters 
year = '2025'
date = '0227' #date format, 'mmdd'
acqu_numb = '0002' #acquisition number 

ordi = 'dell_vasco'

if ordi == 'babasse':
    path2data = os.path.join('E:/Data/',date,'Geophones/')
elif ordi == 'dell_vasco':
#    path2data = f'B:/Data/{date}/Geophones/'
    path2data = f'D:/copie_BicWin25_geophones/Data/{date}/Geophones/'
elif ordi == 'adour':
    path2data = f'/media/turbots/Shack25/Data/{date}/Geophones/'


fig_folder = path2data + 'Figures/'  # folder where figures are saved 
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

if ordi=='babasse':    
    geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
elif ordi=='dell_vasco':
    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
elif ordi=='adour':
    geophones_table_path = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'

channel = 2  # 0 for E, 1 for N, 2 for Z. 
composante = 'Z'

#files need to be organised as: data/0210/Geophones/0001/minised files

geophones_spacing = 3 # space between geophones, in meters 
signal_length = 1 # duration in seconds 
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
img_quality = 1000 # dpi to save images 

plt.rcParams.update({
    "text.usetex": False})

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
    # Creation matrice U S V,  D???
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

#----------------------------------------------------------------------------------------------

def wavenumbers_stein( rho_ice, h, E, nu,freq,c_w,rho_w):
    """ This function computes the wave vectors associated to a given array of frequencies
    It takes as arguments : 
        - rho_ice : ice density 
        - h : a given thickness of ice 
        - E : Young modulus of ice
        - nu : Poisson coefficient of ice
        - freq : an array of frequencies, to which will correspond wave vectors 
        - c_w : waves phase velocity
        - rho_w : water density 
        
    The function returns : 
        - k_QS : wave vectors of the flexural mode
        - k_QS0 : wave vecotrs of the acoustic mode
        - k_SH0 : wave vectors of the shear mode
        - cphQS : phase velocity of the flexural mode"""
    
    
    g = 9.81
    G = E/(2*(1+nu))
    cS0 = np.sqrt(E/(rho_ice*(1-nu**2))) # celerity of longitudinal wave
    cSH0 = np.sqrt(G/rho_ice) # celerity of shear wave 
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    k = np.linspace(1e-6,10,200000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for kf in range(len(freq)):
        omeg = 2*np.pi*freq[kf]
        if omeg == 0:
            flag = 1;
            idx_flag = kf;
        else:
            cph = omeg/k # phase velocity
            # Ludovic version
            func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            # Sebastien version (Stein 1998)
            # func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_ice/D + pow(omeg/cph,4)
            
            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[kf] = (np.where(np.diff(np.signbit(func)))[0]) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero] # wave vector associated to flexural mode 
    if flag:
        k_QS[idx_flag] = 0
        
    k_QS0 = freq/cS0*2*np.pi # wave vector associated to longitudinal wave
    k_SH0 = freq/cSH0*2*np.pi   # wave vector associated to shear wave
    cphQS = freq/k_QS*2*np.pi # phase velocity of the flexural wave

    
    return k_QS, k_QS0, k_SH0, cphQS

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
chan = channel
selected_stream = seismic_data_streams[(idx_geophone - 1)*3 + chan]
print(selected_stream[0])


fig, ax = plt.subplots(figsize = fig_size)
ax.plot(datetime_values,selected_stream[0].data / max(np.abs(selected_stream[0].data)))

#%% tests pour detecter les chocs
'''
curdir = os.getcwd()
os.chdir('Z:/vasco/Geophones/BicWin2025/python_functions/')
from sismo_analysis import find_9shocks
os.chdir(curdir)

typ_sig = selected_stream[0].data[90600:91600]
find_9shocks(selected_stream[0].data[:14000],typ_sig)

'''
"""
plt.plot(selected_stream[0].data / max(np.abs(selected_stream[0].data)))
plt.show()

liste_signaux = []
liste_signaux.append(selected_stream[0].data[90600:91600])
liste_signaux.append(selected_stream[0].data[95500:96500])
liste_signaux.append(selected_stream[0].data[100200:101200])
liste_signaux.append(selected_stream[0].data[108000:109000])
liste_signaux.append(selected_stream[0].data[108000:109000])
liste_signaux.append(selected_stream[0].data[113000:114000])

liste_signaux.append(selected_stream[0].data[161000:162000])
liste_signaux.append(selected_stream[0].data[166500:167500])
liste_signaux.append(selected_stream[0].data[171500:172500])
liste_signaux.append(selected_stream[0].data[183000:184000])
liste_signaux.append(selected_stream[0].data[189000:190000])
liste_signaux.append(selected_stream[0].data[194000:195000])

liste_couleurs = ['tab:blue','tab:blue','tab:blue','tab:orange','tab:orange','tab:orange','tab:blue','tab:blue','tab:blue','tab:orange','tab:orange','tab:orange']

plt.figure()
for i in range(len(liste_signaux)):
    plt.plot(np.fft.fftshift(np.abs(np.fft.fft(liste_signaux[i])))/np.max(np.fft.fftshift(np.abs(np.fft.fft(liste_signaux[i])))),color=liste_couleurs[i])
plt.show()
array_signaux = np.array(liste_signaux)
plt.figure()
plt.plot(np.mean(np.abs(np.fft.fftshift(np.fft.fft(array_signaux[:3]))),axis=0))
plt.plot(np.mean(np.abs(np.fft.fftshift(np.fft.fft(array_signaux[3:6]))),axis=0))
plt.show()
"""
#%% FFT of the signal 
fs = 1000

FFT = np.fft.fft(selected_stream[0].data - np.mean(selected_stream[0].data))
n = np.size(selected_stream[0].data)
frequencies = np.fft.fftfreq(n,1/fs)

x_fft = frequencies[:len(frequencies)//2]
y_fft = FFT[:len(FFT)//2]

fig,ax = plt.subplots()
ax.plot(x_fft,np.abs(y_fft))



##############################################################################################
#%% ---------------------- Plot all channels for a single geophone ---------------------------
##############################################################################################
plt.rcParams['text.usetex'] = False

idx_geophone = 2 # index of the geophone to be plotted

selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]

fig, ax = plt.subplots(3,figsize = fig_size)
for k, chan in enumerate([2,0,1]):
    current_stream = selected_streams[chan]
    print(current_stream[0])
    ax[k].plot(first_trace.times(),current_stream[0].data / max(np.abs(current_stream[0].data) ),label = channel_dic[chan])
    # ax[k].plot(first_trace.times(),current_stream[0].data*scale_factor,label = channel_correspondance[chan])
    ax[k].legend()
    
#    ax[k].set_ylabel(r'$V/V_{max}$')
    ax[k].set_ylabel('$V/V_{max}$')
    ax[k].set_ylim([-1.1, 1.1])
    # fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
    # # # Create an instance of ZoomHandler
    # # zoom_handler = ZoomHandler(ax, time_vector, data_vector)
    # # fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
    # # # Enable interactive labels using mplcursors
    # # mplcursors.cursor(hover=True)
    # # Adjust the backend to make it work better in Spyder
    # plt.ion()
    # plt.show(block=True)  # Use block=True to make it work better in Spyder

ax[2].set_xlabel('$t \; \mathrm{(s)}$')
fig.tight_layout()



#%% Set xlim for all subplots

Mydate = date(2024,2,10)
time_start = datetime.combine(Mydate,time(18,47,7))
time_end = datetime.combine(Mydate,time(18,47,35))
xlimits = [time_start, time_end]

for k in range(0,3):
    ax[k].set_xlim([time_start, time_end])
    
#%% Save current figure
fig.tight_layout()

figname = 'Stream_ZEN_geophone_' + str(idx_geophone) + '_acqu_' + acqu_numb
figname = fig_folder + figname

plt.savefig(figname + '.pdf',dpi = img_quality,bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality,bbox_inches = 'tight')


##############################################################################################
#%%----------------------- PLOTTING SELECTED CHANNEL FOR WHOLE RECORDING ------------------ 
##############################################################################################
 
# c'est avec ce plot qu'on va choisir les temps pour les sources

#channel = 1 #0 for E, 1 for N, 2 for Z. 

# selected indices of the corresponding channel
# seismic_data_G16 = seismic_data_streams[45:]
# current_stream  = seismic_data_G16[channel] 
# print(current_stream[0])
# ax.plot(datetime_values, 
#         current_stream[0].data / max(np.abs(current_stream[0].data) ), 
#         label=f"Stream {k}")
# selected_indices = range(channel, len(s²eismic_data_streams) - 3,3)

selected_indices = range(channel, len(seismic_data_streams),3)
fig, ax = plt.subplots(figsize = fig_size)
for k, idx in enumerate(selected_indices):
    current_stream = seismic_data_streams[idx]
    print(current_stream[0])
    Norm = max(np.abs(current_stream[0].data))
    Norm = np.std(current_stream[0].data)
    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) ) + geophones_spacing*k + geophones_spacing, 
#            current_stream[0].data / Norm + geophones_spacing*k**1.5 + geophones_spacing, # cellule Aurore pour rendre plus visible les pics
            label=f"Stream {k}")
    
    # ax.plot(datetime_values, 
    #         (current_stream[0].data-np.mean(current_stream[0].data))  + geophones_spacing*k + geophones_spacing, 
    #         label=f"Stream {k}")
    

#fig.suptitle(f"Seismic Data - {start_time_utc.strftime('%Y-%m-%d %H:%M:%S')}", fontsize=16)
# Create an instance of ZoomHandler
zoom_handler = ZoomHandler(ax, time_vector, data_vector)
fig.canvas.mpl_connect('button_press_event', zoom_handler.on_click)
# Enable interactive labels using mplcursors
mplcursors.cursor(hover=True)
# Adjust the backend to make it work better in Spyder
plt.ion()
plt.show(block=True)  # Use block=True to make it work better in Spyder

#%% Select xlim of the plot 
Mydate = date(2024,2,11)
time_start = datetime.combine(Mydate,time(18,47,7))
time_end = datetime.combine(Mydate,time(18,47,35))
ax.set_xlim([time_start, time_end])

#%% Once correclty zoomed, we can save the figure 
figname = 'Streams_all_geophones_' + channel_dic[channel] + acqu_numb
figname = fig_folder + figname

plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality ,bbox_inches = 'tight')


#################################################################################

#%%----------------------- SELECTING TIME RANGE FOR PROCESSING ------------------ 
#################################################################################
# Select signal segments manually for each track and each source 
# Store selected times in a dictionnary and save it as a .pkl file 


# open .pkl file and load dictionnary 
filename = 't1_to_time_' + date + '_' + year + '.pkl'
base = path2data
file2save = base + filename
if os.path.isfile(file2save):
    print('Time dictionnary already exists')
    with open(file2save,'rb') as f:
        time_dict = pickle.load(f)
else:
    print('No time dictionnary saved yet')
    time_dict = {}


composante = 'Z'
# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:14:11.20")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:17:58.20")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:22:04.00")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:28:54.50")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:29:00.30")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-21T16:30:31.00")


"""
composante = 'E' #Z , E or N -> direction de la source

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:26:36.60")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:27:56.80")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:29:03.90")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:32:49.30")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:34:02.20")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:35:22.20")
"""

"""
composante = 'N' #Z , E or N -> direction de la source

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:26:53.90")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:28:15.30")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:29:26.30")

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:33:13.90")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:34:30.80")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2025-02-27T18:35:41.20")
"""
# Save t0 dictionnary in pickle file 

savet0 = input('are you sure you want to save t0 dict ? y/n')

if savet0=='y':
    with open(file2save, 'wb') as f:
        pickle.dump(time_dict, f)
    print('Time dictionnary saved')
else:
    pass


###########################################################
#%% -------------- Compute FK data ----------------------
###########################################################

signal_length = 1 # duration in seconds

# load data of intial times 
#composante = 'Z'
#channel = 2
# ch = channel_dic[channel]
flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave
direction = 2 # 1 ou 2 
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
if ordi=='babasse': 
    base = f'E:/Data/{date}/Geophones/'
elif ordi=='dell_vasco':
    base = f'B:/Data/{date}/Geophones/'
    #base = f'D:/copie_BicWin25_geophones/Data/{date}/Geophones/'

pkl_path = base + 't1_to_time_' + date + '_' + year  + '.pkl'

with open(pkl_path, 'rb') as f:
    loaded_data = pickle.load(f)


# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)
t1_values = [t1,t2,t3]

# Create a matrix to store the seismic data
selected_indices = np.arange(channel, len(seismic_data_streams),3)
ta = t1 + signal_length 
num_samples = int(signal_length * first_trace.stats.sampling_rate) # Number of points generated in time_vector 
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

num_traces = len(selected_indices)

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((num_traces, third_dim_len, num_samples)) # matrix dim : geo_indices x nb_sources x time

# add G16 in first place 
# idx_G16 = selected_indices[-1]
# stream = seismic_data_streams[idx_G16]
# print(stream[0])
# for t1_index, t1 in enumerate(t1_values):
#     for j, trace in enumerate(stream):
#         start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
#         # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
#         end_sample = start_sample + num_samples
#         seismic_matrix[0, t1_index, :] = trace.data[start_sample:end_sample]


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
#plt.grid(True)


#%%  Calling fn_SVD

#####################################
##### IMPORTANT PARAMETERS ##########
#####################################

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
nb_stacking = 1 # number of times we want to stack the FK plot horizontally
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

fig, ax1 = plt.subplots(1, 1, figsize=fig_size)
# ax1 = axes  # No need for indexing when there's only one subplot
# c1 = ax1.contourf(F, K_uwp, FK_uwp, cmap='gnuplot2', vmin=vmin, vmax=vmax)

# ax1.set_xlabel('Frequency (Hz)', fontsize=font_size_medium)
# ax1.set_ylabel('Wavenumber (rad/m)', fontsize=30)
# ax1.tick_params(axis='x', labelsize=25, length=10)
# ax1.tick_params(axis='y', labelsize=25, length=10)
# #ax1.set_xlim(xmin, xmax)
# #ax1.set_ylim(ymin, ymax)

# ax1.set_title(f"Spectrum with SVD filter - {date}-{year}-{acqu_numb}\n"
#               f"Channel: {ch}, Source: {composante}, Direction: {direction}", fontsize=30)
# colorbar = plt.colorbar(c1, ax=ax1)
# colorbar.set_label('Spectrum Magnitude', fontsize=25)
# colorbar.ax.tick_params(labelsize=20)

#-----------------------------------------------
# Set parameters 
threshold = 0.2 # minimum relative amplitude to detect a local maximum in the dispersion relation FK plot 
precision_k = [0.03,0.025] # precision over k values (when searching for a local maximum)
prec = precision_k[flexure_wave]
# 0: horizontal_wave 1: flexure_wave
semirange_freq = 5 # semi-range of frequencies for horizontal waves
nb_points = 12 # points to select on graph FK flexural
#-----------------------------------------------

plt.show()
if flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
    #ax1.set_ylim([-1.5, 1.5])
#    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
#    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    ax1.set_ylabel('$f$ $(Hz)$',labelpad = 5)
    ax1.set_xlabel('$k$ $(m^{-1})$',labelpad = 5)
    ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')

    ax1.set_ylim([ymin, ymax])
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
#    plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
#    plt.ylabel(r'$f \; \mathrm{(Hz)}$')
    plt.ylabel('$f$ $(Hz)$')
    plt.xlabel('$k$ $(m^{-1})$')
    
    
    plt.legend()
    plt.show()

#--------------------------------------------------------------------------------
elif horizontal_wave:

    # Computes extents for imshow
    x = extents(k)
    y = extents(f)

    c1 = ax1.imshow(FK_normalized.T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

    #ax1.set_ylim([-1.5, 1.5])
#    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
#    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    ax1.set_ylabel('$f$ $(Hz)$',labelpad = 5)
    ax1.set_xlabel('$k$ $(m^{-1})$',labelpad = 5)
    
    # ax1.set_title('Spectrum with SVD filter')
#    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}|_{max}}(f,k)$')
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
    ax1.set_ylabel('$f$ (Hz)')
    ax1.set_xlabel('$k$ ($m^{-1}$)')
#    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
#    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    # ax1.set_title('Spectrum with SVD filter')
#    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}(f,k)$')
#    plt.colorbar(c1, ax=ax1, label= '$|\hat{s}|(f,k)$')
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
#    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
#    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    ax1.set_ylabel('$f$ (Hz)')
    ax1.set_xlabel('$k$ ($m^{-1}$)')
    # ax1.set_title('Spectrum with SVD filter')
#    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}(f,k)$')
    ax1.set_ylim([ymin, ymax])
    
    wave_type = 'QS'
    figname = fig_folder + 'Flexural_FKplot_' + wave_type + '_acqu_' + acqu_numb + '_dir'  + str(direction)

    plt.savefig(figname + '.pdf', dpi = 300, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = 300, bbox_inches = 'tight')
    

#%% Save data in a dictionnary  

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

if horizontal_wave :
    key_dir = 'dir' + str(direction)
    key = keys_phase_velocity[composante]
    key_uncertainty = 'u' + key
    s[key_dir][key] = C[0] # value 
    s[key_dir][key_uncertainty] = C[1] # standard deviation 

#%% write phase velocity in a csv file for a given direction 

# csv_file = path2data + 'Phase_velocity_acqu_' + acqu_numb + '_direction_' + str(direction) + '_with_std.csv'
# with open(csv_file, 'a') as csvfile: 
#     key_direction = 'dir' + str(direction)
#     writer = csv.DictWriter(csvfile, fieldnames = s[key_direction].keys()) 
#     writer.writeheader() 
#     for key in s[key_direction].keys(): 
#         csvfile.write("%.3f"% s[key_direction][key] + ',')
#     csvfile.write('\n')

#%% Save phase velocity data 

with open (pkl_file,'wb') as pfile:
    pickle.dump(s,pfile)
    
    
    
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

#%% Save (f,k) points associated to QS mode 
saveflexural = input('are you sure you want to overwrite flexural data?')
if saveflexural=='y':
    file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' + str(direction) +'.pkl'
    with open(file2save, 'wb') as file:
        pickle.dump([f_mode, k_mode], file)
     
     
#%% Invert ice thickness by hand 
if flexure_wave:
    
    print('Select 10 points on FK of flexural wave (QS mode)')
    points = plt.ginput(10, timeout=-1)
    
    # Extract x and y coordinates of the selected points
    f_mode, k_mode = zip(*points)
    # Plot a line between the two selected points
    plt.plot(f_mode, k_mode, linestyle='--', color='r', label='Line between points')
    
    file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' + str(direction) + '_hand_selection.pkl'
    with open(file2save, 'wb') as file:
         pickle.dump([f_mode, k_mode], file)
     


#%% --------- INVERT ICE THICKNESS USING FK MATRIX -------------------
#######################################################################

rho_ice = 917 
E_fit = 4.7449e9 # 2.43e9
nu_fit = 0.2821

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 
#path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones' # path 2 data points saved from FK flexural
# acqu_numb = '0002' # acquisition number 

# load selected points on the FK plot 
file2load = path2data +'/' +year +'_'+ date + '_acq' + acqu_numb + 'disp_QS_dir' + str(direction) + '.pkl'
with open(file2load, "rb") as filename:
    data = pickle.load(filename)
f_mode1 = data[1] 
k_mode1 = data[0] 
fig, ax = plt.subplots()
ax.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points')    
    
f_mode1 = f_mode1[1:]
k_mode1 = k_mode1[1:]

# file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir2.pkl'
# with open(file2load, "rb") as f:
#     data = pickle.load(f)
# f_mode2 = data[0] 
# k_mode2 = data[1] 


#plt.plot(f_mode2, k_mode2, linestyle='--', color='r', label='Line between points') 

f_mode = f_mode1
k_mode = k_mode1

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h[i], E_fit, nu_fit,f_mode,c_w,rho_w)
    error = np.sum(np.power((k_mode-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
   
# computes wavevectors k 
k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h_ice, E_fit, nu_fit,f_mode,c_w,rho_w)

ax.plot(k_mode_synthetic, f_mode, color='g', label='Line between points')    
ax.plot(k_mode, f_mode, linestyle='--', color='r', label='Line between points') 

 
ax.plot(h,l2_norm)

print(f'Ice thickness : {h_ice} m')
# %%
