# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:44:00 2024

@author: sebas

This script is the second version of ice thickness inversion, from geophones signals. 
This is the most recent version, but we here assume that WE KNOW THE ICE DENSITY 

Different steps : 
    - loading of geophones data for a given acuqisition along a line 
    - load the dates at which sources have been performed. These dates have already been selected using an other code : 
            'extract_initial_time.py'
    
    For each direction (from Geophone 1 -> 16) or (from Geophone 16 -> 1), we can perform the following steps
    - SVD decomposition for each track : E (shear waves), N (longitudinal waves) 
    - From in-plane waves velocity : computation of Young modulus and Poisson coefficient
    - SVD decomposition for Z (out-of-plane flexural waves)
    - Inversion of ice thickness 
    
N.B. : The inversion of the ice thickness can then be done using the averaged value of E and nu from both directions 


"""

%matplotlib qt

import os
import shutil
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime, time 
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate 
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import pickle
import csv 

# import seb
# from seb import pickle_m
##############################################################
#%% --------- SETTING PARAMETERS AND PATHS -------------------
##############################################################
date = '0211'
year = '2024'
path2data = 'F:/Rimouski_2024/Data/' + year + '/' + date +'/Geophones/'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'

#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0001' 
channel_correspondance = ['E','N','Z']
gain = 12 # gain in dB 
scale_factor = 10**(-gain/10)
geophones_spacing = 3 # space between geophones, in meter 

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


fig_folder = path2data + acqu_numb + '/' + 'Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
######################################################
#%% ----------------------- FUNCTIONS ------------------ 
######################################################

#warnings.filterwarnings("ignore", category=np.ComplexWarning)
def fn_svd(signals, fs, xs, rang,name, issaving, *varargin):
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
        
    # additional argument if needed 
    if varargin:
        if varargin[0] == 'threshold':
            threshold_user = varargin[1] # threshold between signal and noise in dB of a peak value
            # we keep left singular vectors, only if the associated singular values intensity is higher than this threshold 
        elif varargin[0] == 'rang':
            rang = varargin[1]
        else:
            print('varargin(1) unknown')
            return
    else:
        print('varargin empty')

    Nreceiv, Nt, Nemit = signals.shape # number of recepter, time steps and sources  

    # time domain fft
    Nf = 2048
    f = (np.arange(Nf) / Nf) * (fs if fs else 1)
    f_axename = 'f/fs' if not fs else 'f'
    SIGNALS = fft(signals, Nf, axis=1) # fft over time 
    SIGNALS = SIGNALS[:, :Nf + 1, :] 

    # svd
    U = np.zeros((Nreceiv, Nf, Nemit), dtype=complex) # left singular vectors
    S = np.zeros((Nemit, Nf, Nemit), dtype=complex) # diagonal matrix with singular values (in decreasing order)
    V = np.zeros((Nemit, Nf, Nemit), dtype=complex) # right singular vectors 
    D = np.zeros((Nemit, Nf), dtype=complex) # singular values for each frequency 

    for ii in range(Nf): # loop over frequencies
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False) 
        D[:, ii] = np.diag(S[:, ii, :]) # extract singular values for the associated frequency 

    for ne in range(Nemit): 
        titi = 20 * np.log10(np.abs(D[ne, :]) / np.max(np.abs(D[0, :])))
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
        # print(sigmacutsup)              
    # plt.plot(f, sigmacutsup)

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
    # projections = projections[:Nk + 1, :, :]
    # projections = projections[:, :-1, :]
    projections_sum = np.zeros((Nk, Nf, Nemit))

    # for kemit in range(Nemit):
    #     for ii in range(Nf + 1):
    #         projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]) ** 2


    for kemit in rang:
        for ii in range(Nf):
            max_value = 1  # np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2

    projections_sum = np.abs(np.mean(projections_sum, axis=2)) # average over number of emitters 

    return f, k, projections_sum

#--------------------------------------------------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------------------------------------------

def convert_to_utc_times(trace, time_vector):
    """ Create a UTC-time time vector for a given trace 
    Inputs :
        - trace : trace object (obspy) 
        - time_vector : array of time elapsed since beginning of the acquisition 
    Output : 
        - array of UTC time (datetime type)"""
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]

#--------------------------------------------------------------------------------------------------------------------------

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
            print(f"Warning: No entry found for {last_five_digits} in the geophones table.")

    # Sort traces based on the new station codes
    sorted_stream = sorted(stream, key=lambda trace: trace.stats.station)

    return sorted_stream

#--------------------------------------------------------------------------------------------------------------------------

def sort_key(trace):
    """ Create a sorting key based on stream stations """
    return trace[0].stats.station

#--------------------------------------------------------------------------------------------------------------------------

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

    k = np.linspace(1e-6,20,200000)
    
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

#--------------------------------------------------------------------------------------------------------------------------

def hydro_elastic(f,h,rho_ice,E,nu):
    """ Computes the hydro-elastic wave vector associated to an array frequency. 
    The function takes the following arguments : 
        - f : frequency arrray 
        - h = thickness of ice 
        - rho_ice = density of ice
        - E = Young modulus of ice
        - nu = Poisson coefficient 
    """
    
    omega = 2*np.pi*f
    k = omega**2 * (12*rho_ice*(1 - nu**2))/(E*h**3)
    k = np.power(k,1/5)
    
    return k

#--------------------------------------------------------------------------------------------------------------------------

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

####################################################
#%% Plot a single channel for a given geophone
####################################################
channel = 0
idx_geophone = 3

current_stream = seismic_data_streams[(idx_geophone - 1)*3 + channel]
print(current_stream[0])
normalized_data = current_stream[0].data * scale_factor

fig, ax  = plt.subplots(figsize = fig_size)
ax.plot(first_trace.times(),normalized_data)
ax.set_ylabel('$V \; \mathrm{(m.s^{-1})}$')
ax.set_xlabel('$t \; \mathrm{(s)}$')


#%%
figname = fig_folder + 'Scaled_signal_channel_' + channel_correspondance[channel] + 'idx_geophone_' + str(idx_geophone) +'_' + acqu_numb + '_gain_' + str(gain)
plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality, bbox_inches = 'tight')

##############################################################################################
#%% ---------------------- Plot all channels for a single geophone ---------------------------
##############################################################################################

idx_geophone = 2 # index of the geophone to be plotted

selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]
plt.rcParams.update({
    "text.usetex": True})
fig, ax = plt.subplots(3,figsize = fig_size)
for k, chan in enumerate([2,0,1]):
    current_stream = selected_streams[chan]
    print(current_stream[0])
    ax[k].plot(first_trace.times(),current_stream[0].data / max(np.abs(current_stream[0].data) ),label = channel_correspondance[chan])
    # ax[k].plot(first_trace.times(),current_stream[0].data*scale_factor,label = channel_correspondance[chan])
    ax[k].legend()
    
    ax[k].set_ylabel(r'$V/V_{max}$')
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

Mydate = date(2024,2,11)
time_start = datetime.combine(Mydate,time(18,47,7))
time_end = datetime.combine(Mydate,time(18,47,35))
xlimits = [time_start, time_end]

for k in range(0,3):
    ax[k].set_xlim([time_start, time_end])
    
#%% Save current figure
fig.tight_layout()

figname = 'Stream_ZEN_geophone_' + str(idx_geophone) + '_acqu_2'
figname = fig_folder + figname

plt.savefig(figname + '.pdf',dpi = img_quality,bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality,bbox_inches = 'tight')


#%%----------------------- PLOTTING SELECTED CHANNEL FOR WHOLE RECORDING ------------------ 
##############################################################################################

channel = 2 #0 for E, 1 for N, 2 for Z. 
# selected indices of the corresponding channel 
selected_indices = range(channel, len(seismic_data_streams),3)

fig, ax = plt.subplots(figsize = fig_size)
for k in selected_indices:
    current_stream = seismic_data_streams[k]
    print(current_stream[0])
    # Plot the data for each stream using precalculated datetime values
    ax.plot(datetime_values, 
            current_stream[0].data / max(np.abs(current_stream[0].data) ) + k, 
            label=f"Stream {k}")
    

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
figname = 'Streams_all_geophones_' + channel_correspondance[channel] + acqu_numb
figname = fig_folder + figname

plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png',dpi = img_quality ,bbox_inches = 'tight')

#%% OPTIONAL - FREQUENCY SPECTRUM VS WAVENUMBER for the selected signal
########################################################################

Nfft_t= 1024
Nfft_k = 1024
 
fs = 1000
ks = 2*np.pi/geophones_spacing
space_vector = np.arange(0,seismic_matrix.shape[0],geophones_spacing)

# Define the spatial and temporal frequencies
wavenum = np.fft.fftfreq(Nfft_k, d=1/ks)
freq = np.fft.fftfreq(Nfft_t, d=1/fs)

# Create a 2D Fourier transform of the seismic matrix
seismic_fft = np.fft.fft2(seismic_matrix, (Nfft_k, Nfft_t))

# Shift zero frequency component to the center
seismic_fft_shifted = np.fft.fftshift(seismic_fft)

fig, ax = plt.subplots()
# Plot the 2D Fourier transform
plt.imshow(np.transpose(np.abs(seismic_fft_shifted)), extent=(wavenum.min(), wavenum.max(),freq.min(), freq.max()), aspect='auto', cmap='viridis')
#plt.colorbar(label='Amplitude')
plt.ylabel(r'$f\; \mathrm{(Hz)}$')
plt.xlabel(r'$k \; \mathrm{(m^{-1})}$')
plt.title('2D Fourier Transform of Seismic')
plt.show()
ax.set_ylim([0 , fs/4])


##################################################################################
#%%----------------------- CONSTRUCTING MATRIX FOR FK WITH SVD ------------------
##################################################################################

signal_length = 1.0 # duration in seconds

# load data of intial times 
composante = 'Z'
channel = 2
# ch = channel_dic[channel]
flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave
direction = 1 # 1 ou 2 
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
 
base = 'C:/Users/sebas/Desktop/Amundsen_RA_2024/Data/2024/0909/Geophones/'
pkl_path = base + 't1_to_time_' + date + '_' + year + '.pkl'
with open(pkl_path, 'rb') as f:
    loaded_data = pickle.load(f)


# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)
t1_values = [t1,t2,t3]


ta = t1 + signal_length 
num_samples = int(signal_length * first_trace.stats.sampling_rate) # Number of points gereted in tim_vector (~1000)
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

# selected indices of the corresponding channel 
selected_indices = range(channel, len(seismic_data_streams),3)

# Create a 3D matrix to store the seismic data
num_traces = len(seismic_data_streams)
num_samples = int(signal_length * first_trace.stats.sampling_rate)  

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((len(selected_indices), third_dim_len, num_samples)) # matrix dim : geo_indices x nb_sources x time

# Extract the relevant data for each trace and populate the 3D matrix
for i, stream_index in enumerate(range(channel, len(seismic_data_streams), 3)):
    stream = seismic_data_streams[stream_index]
    
    for t1_index, t1 in enumerate(t1_values):
        for j, trace in enumerate(stream):
            start_sample = int((t1 - trace.stats.starttime) * trace.stats.sampling_rate)
            # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index, :] = trace.data[start_sample:end_sample]

print('Matrix computed for channel : ' + channel_correspondance[channel])

#%% Plot the data (superposition of three sources)

fig, ax = plt.subplots(figsize = fig_size)
for k in range(seismic_matrix.shape[0]):
    for t1_index in range(seismic_matrix.shape[1]):
        ax.plot(time_vector, seismic_matrix[k, t1_index, :] / max(np.abs(seismic_matrix[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams at Different t1 Values')
plt.show()

#######################################################################
#%% OPTIONAL ---------------------- INTERPOLATING 3D MATRIX ------------------
#######################################################################

# Interpolate along the first dimension

xr = np.arange(0, 46, 3)
xs = np.arange(5, 12, 3)
xr_new = np.arange(0, 46, 1)
t_new = time_vector


xx, tt = np.meshgrid(xr, t_new, sparse=True)


seismic_matrices_interp = []

for k in range(seismic_matrix.shape[1]):  # Loop over the second dimension of seismic_matrix
    data_to_interpolate = np.transpose(seismic_matrix[:, k, :])
    f = interpolate.interp2d(xr, t_new, data_to_interpolate, kind='cubic')
    seismic_matrix_interp = f(xr_new, t_new)
    
    # Append the interpolated matrix to the list
    seismic_matrices_interp.append(seismic_matrix_interp)

# Convert the list to a 3D numpy array
seismic_matrices_interp = np.stack(seismic_matrices_interp, axis=1)
seismic_matrices_interp = np.transpose(seismic_matrices_interp, (2, 1,0))


fig, ax = plt.subplots()
for k in range(seismic_matrices_interp.shape[0]):
    for t1_index in range(seismic_matrices_interp.shape[1]):
        ax.plot(time_vector, seismic_matrices_interp[k, t1_index, :] / max(np.abs(seismic_matrices_interp[k, t1_index, :])) + 3*k,
                label=f"Stream {k}, t1={t1_values[t1_index]}")

ax.set_xlabel('Time (UTC)')
ax.set_ylabel('Normalized Seismic Data')
ax.set_title('Seismic Data for Selected Streams at Different t1 Values')
plt.show()

#############################################################################
#%% -------------- COMPUTE FK-MATRIX USING SVD METHOD -------------------------
#############################################################################

rang = [0,1,2] # rank of singular values 
signals = np.transpose(seismic_matrix, (0, 2, 1)) # #1 = geo_indices, #2 = time, #3 = nb_sources 
#signals = np.transpose(seismic_matrices_interp, (0, 2, 1))

f, k, FK = fn_svd(signals, fs, geophones_spacing,rang,'ExampleName', 0, 'threshold',-90)
# dimension 1 : k
# dimension 2 : f
F, K = np.meshgrid(f, k)

###############################################################
#%% ------------ FK-PLOT USING SVD METHOD ---------
###############################################################

###########################################################
# /!\ /!\ NEEDS FLIPUD DEPENDING ON DIRECTION OF PROPAGATION /!\ /!\
###########################################################

if direction == 2:
    FK = np.flipud(FK/np.max(FK))
else : 
    FK = FK/np.max(FK)


vmin = 0     # Minimum value for colormap
vmax = 1 # Maximum value for colormap (adjust as needed)
fig, ax1 = plt.subplots(1, 1, figsize = fig_size)

# Computes extents for imshow
x = extents(k)
y = extents(f)

c1 = ax1.imshow(FK.T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

#ax1.set_ylim([-1.5, 1.5])
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
# ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}|_{max}}(f,k)$')
ax1.tick_params(axis='both', which='major', pad=7)
ax1.set_ylim([0, 400])

#%%
figname = fig_folder + 'FK_plot_dir' + str(direction + 1) + '_' + channel_correspondance[channel]
plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches='tight')
plt.savefig(figname + '.png',dpi = img_quality, bbox_inches='tight')

#########################################################################################
#%% ----------------------- UNWRAPPING SPECTRUM - HORIZONTAL STACKING ------------------
#########################################################################################

nb_stacking = 2 # number of times we want to stack the FK plot horizontally
idx_stacking = 1
FK_uwp = np.vstack((FK, FK))
while idx_stacking < nb_stacking :
    FK_uwp = np.vstack((FK_uwp, FK))
    idx_stacking += 1
    
k_uwp= np.linspace(0,(nb_stacking + 1)*max(k),FK_uwp.shape[0])


vmin = 0     # Minimum value for colormap
vmax = 1# Maximum value for colormap (adjust as needed)
fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))

c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
#ax1.set_ylim([-1.5, 1.5])
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')

ax1.set_ylim([0, 400])

##############################################
#%% Select points to determine shear velocity 
##############################################
 
points_shear = plt.ginput(10, timeout=-1)
k_mode, f_mode = zip(*points_shear)

# fit by a 1D polynome
deg = 1
p,V = np.polyfit(k_mode,f_mode,deg,cov = True )

C_shear = np.zeros((2,1))
C_shear[0] = p[0]*2*np.pi
C_shear[1] = np.sqrt(2*np.pi*np.diag(V)[0]) # standard deviation
print(C_shear)

######################################################
#%% Select points to determine longitudinal velocity 
######################################################

points_longi = plt.ginput(10, timeout=-1)
f_mode, k_mode = zip(*points_longi)

# fit by a 1D polynome
deg = 1
p, V = np.polyfit(f_mode,k_mode,deg,cov = True )

C_longi = np.zeros((2,1))
C_longi[0] = 2*np.pi*p[0]
C_longi[1] = np.sqrt(2*np.pi*np.diag(V)[0])
print(C_longi)

###################################################
#%% Computes Young modulus & Poisson coefficient 
###################################################

rho_ice = 917
nu = np.zeros((2,1))
E = np.zeros((2,1))
nu[0] = 1-2*(C_shear[0]/C_longi[0])**2
E[0] = rho_ice*C_longi[0]**2*(1-nu[0]**2)

# Computes standrad deviation  
nu[1] = 4*C_shear[0]/(C_longi[0]**2) * C_shear[1] + 4*C_shear[0]**2/(C_longi[0]**3) * C_longi[1]
E[1] = 2*rho_ice * (C_longi[0]*(1-nu[0]**2)*C_longi[1] + C_longi[0]**2*nu[0]*nu[1])

print(E*1e-9)
print(nu)


#%% Import Phase velocity data averaged on both directions : 

path2velocity_dict = 'C:/Users/sebas/Desktop/Amundsen_RA_2024/Data/2024/0909/Geophones/'
file_QS0 = os.path.join(path2velocity_dict,'2024_0909_acq0001_cQS0_bidir.pkl')
file_SH0 = os.path.join(path2velocity_dict,'2024_0909_acq0001_cSH0_bidir.pkl')

with open(file_QS0,'rb') as file:
    C_QS0 = pickle.load(file)
    
with open(file_SH0,'rb') as file:
    C_SH0 = pickle.load(file) 

# #####################################################################
#%% Only for Flexural mode (mode Z) Select points on the FK plot 
# #####################################################################

points = plt.ginput(10, timeout=-1)

# Extract x and y coordinates of the selected points
f_mode, k_mode = zip(*points)
# Plot a line between the two selected points
plt.plot(f_mode, k_mode, linestyle='--', color='r', label='Line between points')

file2save = path2data + 'dispersion_QS_acq_' + acqu_numb + '_dir' + str(direction) + '.pkl'
with open(file2save,'wb') as file:
    pickle.dump([f_mode, k_mode], file)
    
#######################################################################
#%% --------- INVERT ICE THICKNESS USING FK MATRIX -------------------
#######################################################################

rho_ice = 917 
E_fit = 2.51e9 # 2.43e9
nu_fit = 0.17

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 
#path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones' # path 2 data points saved from FK flexural
# acqu_numb = '0002' # acquisition number 

# load selected points on the FK plot 
file2load = path2data +'/' +year +'_'+ date + '_acq' + acq_number + 'disp_QS_dir' + str(direction) + '.pkl'
with open(file2load, "rb") as filename:
    data = pickle.load(filename)
f_mode1 = data[1] 
k_mode1 = data[0] 
plt.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points')    
    
f_mode1 = f_mode1[1:]
k_mode1 = k_mode1[1:]

# file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir2.pkl'
# with open(file2load, "rb") as f:
#     data = pickle.load(f)
# f_mode2 = data[0] 
# k_mode2 = data[1] 


plt.plot(k_mode1, f_mode1, linestyle='--', color='r', label='Line between points') 
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

plt.plot(k_mode_synthetic, f_mode, color='g', label='Line between points')    
plt.plot(k_mode, f_mode, linestyle='--', color='r', label='Line between points') 

 
plt.plot(h,l2_norm)

print(h_ice)

###############################################
#%% Save all relevant data into a file pickle
###############################################

pickle_file = path2data + 'Physical_data_acq_' + acqu_numb + '_direction_' + str(direction) + '.pkl'

s = dict()
s['C_shear']= C_shear[0]
s['C_longi'] = C_longi[0]
s['E'] = E[0]
s['nu'] = nu[0]
s['h_ice'] = h_ice
s['rho_ice'] = rho_ice
s['c_w'] = c_w
s['rho_w'] = rho_w

s['uC_shear'] = 2*C_shear[1]
s['uC_longi'] = 2*C_longi[1]
s['uE'] = 2*E[1]
s['unu'] = 2*nu[1]
s['uh_ice'] = h_precision
with open(pickle_file,'wb') as file:
    pickle.dump(s,file)

print( 'Data saved in file : ' + pickle_file)
#%% Save relevant data into a csv file 

csv_file = path2data + 'Physical_data_acqu_' + acqu_numb + '_direction_' + str(direction) + '_fk_with_errors.csv'

with open(csv_file, 'w') as csvfile: 
    writer = csv.DictWriter(csvfile, fieldnames = s.keys()) 
    writer.writeheader() 
    for key in s.keys(): 
        csvfile.write("%.3f"% s[key] + ',')
        
#%% Create an other plot with data points and fitted curve 

frequency_list = np.linspace(1,150,100)
k_flex_theory, k_acoust_theory, k_shear_theory, cphQS = wavenumbers_stein( rho_ice, h_ice, E, nu,frequency_list,c_w,rho_w)

fig, ax = plt.subplots()
FK = FK/np.max(FK)
c1 = ax.contourf(F, K, FK, cmap='gnuplot2', vmin=vmin, vmax=vmax)
ax.scatter(data[0],data[1], label = 'Detected points')
ax.plot(frequency_list,k_flex_theory,color = 'k', label = 'Flexural mode')
plt.colorbar(c1, ax=ax, label= r'Spectrum Magnitude')

ax.set_xlim([0, 100])
ax.grid()
ax.set_xlabel(r'$f \quad \mathrm{(Hz)}$')
ax.set_ylabel(r'$k \quad \mathrm{(m^{-1})}$')
# ax.scatter(f_mode,k_QS0)
# ax.scatter(f_mode,k_SH0)

#############################################################################
#%% Superposition of the Unwrapped FK spectrum and the fitted flexural law 
#############################################################################

frequency_list = np.linspace(1,300,100)
k_flex_theory, k_acoust_theory, k_shear_theory, cphQS = wavenumbers_stein( rho_ice, h_ice, E_fit, nu_fit,frequency_list,c_w,rho_w)
k_hydro_elast = hydro_elastic(frequency_list, h_ice, rho_ice, E_fit, nu_fit)
fig, ax1 = plt.subplots(1, 1, figsize=(18, 9))
vmin = 0 
vmax = 1 

c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)

#ax1.set_ylim([-1.5, 1.5])
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$', labelpad = 5)
# ax1.set_title(r'Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label = r'$\frac{|\hat{s}|}{|\hat{s}_{max}|}(f,k)$')
ax1.plot(k_flex_theory,frequency_list,linestyle = '--', 
         linewidth = 5 ,color = 'white', label = 'Flexural mode')
ax1.tick_params(axis='both', which='major', pad=7)
# ax1.plot(k_hydro_elast,frequency_list,linestyle = '--', 
#          linewidth = 3, color = 'r')
ax1.set_ylim([0 , 300])

#%%
figname = fig_folder + 'Unwrapped_flexural_with_theory_acqu_' + acqu_numb + '_dir'  + str(direction)

plt.savefig(figname + '.pdf', dpi = 1200, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = 1200, bbox_inches = 'tight')
