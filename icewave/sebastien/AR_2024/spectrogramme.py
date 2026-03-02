# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:21:18 2024

@author: sebas
"""

import os
import shutil
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime, time 
import numpy as np
from scipy import interpolate
from scipy.signal import spectrogram
from scipy.signal import ShortTimeFFT, windows, iirnotch, filtfilt
import scipy.integrate as integrate 
from scipy.interpolate import make_interp_spline
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import pickle
import csv 
import math

##############################################################
#%% --------- SETTING PARAMETERS AND PATHS -------------------
##############################################################
date = '0909'
year = '2024'
path2data = os.path.join('E:/Amundsen_RA_2024/Data/',date,'Geophones/')
geophones_table_path = 'C:/Users/sebas/Github/icewave/icewave/sebastien/geophones/geophones_table'

#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0001' 
channel_correspondance = ['E','N','Z']
gain = 12 # gain in dB 
scale_factor = 10**(-gain/10)
geophones_spacing = 4 # space between geophones, in meter 

# Parameters for plots
font_size_medium = 24
font_size_small = 22
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 600 # dpi to save images 


fig_folder = path2data + 'Figures_filtered_data/'  # folder where figures are saved 
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

####################################################
#%% FUNCTIONS DEFINITION 
####################################################
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
    
    fig, ax = plt.subplots()
    for ii in range(Nf):
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False)
        D[:, ii] = np.diag(S[:, ii, :])
    for ne in range(Nemit):
        titi = 20 * np.log10(np.abs(D[ne, :]) / np.max(np.abs(D[0, :])))
        ax.plot(f, titi, label=f'Slice {ne}')
        ax.set_xlabel(r'f \: \rm{(Hz)}')
        ax.set_ylabel('Singular values (in dB of peak value)')
        ax.set_xlim([f[0] , f[Nf//2]])
        

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

    k = np.linspace(1e-12,4.0,100000)
    
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

def sort_key(trace):
    """ Create a sorting key based on stream stations """
    return trace[0].stats.station

#--------------------------------------------------------------------------------------------------------------------------
def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#--------------------------------------------------------------------------------------------------------------------------
def E_nu_inversion(rho_ice,C_longi,C_shear,**kargs):
    """Computes Young modulus E and Poisson coefficient nu from phase velocity of modes QS0 and SH0 
    Inputs : - rho_ice, ice density 
             - C_longi, QS0 phase velocity
             - C_shear, SH0 phase velocity
             - uC_longi, standard deviation of C_longi
             - uC_shear, standard deviation of C_shear
   
    Outputs : - s, dictionnary containing keys 'E','nu','uE' and 'unu' """
    
    s = {}
    s['nu'] = 1-2*(C_shear/C_longi)**2
    s['E'] = rho_ice*C_longi**2*(1-s['nu']**2)

    uC_longi = kargs.get('uC_longi',None)
    uC_shear = kargs.get('uC_shear',None)
    
    if (uC_longi is not None) and (uC_shear is not None):
        # Computes standrad deviation  
        s['unu'] = 4*C_shear/(C_longi**2) * uC_shear + 4*C_shear**2/(C_longi**3) * uC_longi
        s['uE'] = 2*rho_ice * (C_longi*(1-s['nu']**2)*uC_longi + C_longi**2 * s['nu'] * s['unu'])

    return s 





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


#%% Compute spectrogramm of a given geophone 

channel = 2
idx_geophone = 1
channel_correspondance = ['E', 'N', 'Z']

current_stream = seismic_data_streams[(idx_geophone - 1)*3 + channel]
print(current_stream[0])
normalized_data = current_stream[0].data / max(current_stream[0].data) 

fig, ax  = plt.subplots(2,1,figsize = fig_size)
ax[0].plot(first_trace.times(),normalized_data)
ax[0].set_ylabel(r'$V/V_{max}$')
# ax[0].set_xlim([950 , 1200])
# ax[0].set_xlabel('$t \; \mathrm{(s)}$')


# fs = 1000
# [f,t,Sxx] = spectrogram(current_stream[0].data,fs,detrend = 'constant')
# # Sxx dim : #0 = freq, #1 = time 


# fig,ax = plt.subplots(figsize = fig_size)
# c1 = ax.imshow(Sxx, aspect = 'auto', cmap='viridis',origin = 'lower', extent = extents(t) + extents(f), vmin = np.min(Sxx), vmax = np.max(Sxx))
# plt.colorbar(c1, ax=ax, label= r'$S(f,t)$')
# c1.set_clim([0,3])
# ax.set_xlabel(r'$t \: \mathrm{(s)}$')
# ax.set_ylabel(r'$f \: \mathrm{(Hz)}$')
# ax.set_ylim([0, 150])

# other method 
fs = 1000
num_samples = len(normalized_data)
M = 1000 # window size 
w = windows.tukey(M)
hop = int(M*0.1)
SFT = ShortTimeFFT(w, hop , fs , mfft=1024, scale_to='magnitude')
Sx = SFT.stft(normalized_data)  # perform the STFT
t = SFT.t(num_samples)
f = SFT.f
extents_tf = SFT.extent(num_samples, axes_seq = 'tf')
c1 = ax[1].imshow(abs(Sx), aspect = 'auto', cmap='viridis',origin = 'lower', extent = extents_tf, vmin = np.min(abs(Sx)), vmax = np.max(abs(Sx)))
plt.colorbar(c1, ax=ax, label= r'$S(f,t)$')
# c1.set_clim([1e-5,0.005])
ax[1].set_xlabel(r'$t \: \mathrm{(s)}$')
ax[1].set_ylabel(r'$f \: \mathrm{(Hz)}$')
ax[1].set_ylim([0, 100])
# ax[1].set_xlim([950 , 1200])


file2save = fig_folder + 'STFT_geo' + str(idx_geophone) + '_ch' + channel_correspondance[channel] + '_window_tuckey_' + str(M)
plt.savefig(file2save + '.pdf',dpi = 400,bbox_inches = 'tight')
plt.savefig(file2save + '.png',dpi = 400, bbox_inches = 'tight')


#%% Try to filter noise from the Amundsen engine ...

channel = 1
idx_geophone = 1
channel_correspondance = ['E', 'N', 'Z']

current_stream = seismic_data_streams[(idx_geophone - 1)*3 + channel]
print(current_stream[0])
normalized_data = current_stream[0].data / max(current_stream[0].data) 

# select a section where there is only the engine noise 
t = first_trace.times()
mask = t > 0

fig, ax  = plt.subplots(figsize = fig_size)
ax.plot(t[mask],normalized_data[mask])
ax.set_ylabel(r'$V/V_{max}$')

# computes signal FFT
sig = normalized_data[mask]
N = len(sig)
FFT = np.fft.fft(sig,norm = 'forward')
FFT_pos = FFT[:N//2]
FFT_pos[1:] = 2*FFT_pos[1:]
freq = np.fft.fftfreq(N,1/fs)
freq = freq[:N//2]

fig, ax = plt.subplots()
ax.loglog(freq,abs(FFT_pos))
ax.set_xlabel(r'$f \: \rm (Hz)$')
ax.set_ylabel(r'$\hat{V}(f)$')
ax.grid()
# ax.set_xlim([0,200])

file2save = fig_folder + 'FFT_spectrum_geo' + str(idx_geophone) + '_ch' + channel_correspondance[channel]
plt.savefig(file2save + '.pdf',dpi = 400,bbox_inches = 'tight')
plt.savefig(file2save + '.png',dpi = 400, bbox_inches = 'tight')

#%% Coupe_bande filter 

fc = [2.3, 8.33, 16.66, 25, 41.66, 2*16.66, 3*16.66] #0909
# fc = [2.955, 8.33, 16.66, 25, 41.66, 2*16.66, 3*16.66] #1009
fs = 1000
Q = 30.0  # Quality factor
filtered_sig = sig
for f0 in fc:

    # Design notch filter
    b, a = iirnotch(f0, Q, fs)
    filtered_sig = filtfilt(b,a,filtered_sig)
    
    
    
fig, ax  = plt.subplots(figsize = fig_size)
ax.plot(t[mask],filtered_sig)
ax.set_ylabel(r'$V/V_{max}$')

fig, ax = plt.subplots(figsize = fig_size)
ax.plot(t[mask],sig)
ax.set_ylabel(r'$V/V_{max}$')

#%% Plot 3 channels for a single geophone and filter it

idx_geophone = 4
N = len(first_trace.times())
matrix_geo = np.zeros((3,N))

for k in range(3):
    current_stream = seismic_data_streams[(idx_geophone - 1)*3 + k]
    print(current_stream[0])
    matrix_geo[k,:] = current_stream[0].data

plt.rcParams.update({
    "text.usetex": True})
fig, ax = plt.subplots(3,figsize = fig_size)
for k, chan in enumerate([2,0,1]):
    ax[k].plot(first_trace.times(),matrix_geo[chan,:] / max(np.abs(matrix_geo[chan,:])),label = channel_correspondance[chan])
    # ax[k].plot(first_trace.times(),current_stream[0].data*scale_factor,label = channel_correspondance[chan])
    ax[k].legend(loc = 'upper left')
    
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

file2save = fig_folder + 'Channels_ZEN_' + year + '_' + date + '_G' + str(idx_geophone) + '_raw'
plt.savefig(file2save + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(file2save + '.png', dpi = img_quality, bbox_inches = 'tight')

# filter data and plot it 
file2load = path2data + 'Frequencies_Amundsen_' + year + '_' + date + '_channel_Z' + '.pkl'
with open(file2load,'rb') as pfile:
    fc = pickle.load(pfile)


fs = 1000
Q = 30.0  # Quality factor
matrix_geo_filt = matrix_geo
for f0 in fc:

    # Design notch filter
    b, a = iirnotch(f0, Q, fs)
    matrix_geo_filt = filtfilt(b,a,matrix_geo_filt)
    

fig, ax = plt.subplots(3,figsize = fig_size)
for k, chan in enumerate([2,0,1]):
    ax[k].plot(first_trace.times(),matrix_geo_filt[chan,:] / max(np.abs(matrix_geo_filt[chan,:])),label = channel_correspondance[chan])
    # ax[k].plot(first_trace.times(),current_stream[0].data*scale_factor,label = channel_correspondance[chan])
    ax[k].legend(loc = 'upper left')
    
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

file2save = fig_folder + 'Channels_ZEN_' + year + '_' + date + '_G' + str(idx_geophone) + '_filtered'
plt.savefig(file2save + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(file2save + '.png', dpi = img_quality, bbox_inches = 'tight')





#%% Check FFT of several geophones 

channel = 0
fs = 1000
N = len(first_trace.times())
selected_indices = range(channel,len(seismic_data_streams),3)

raw_matrix = np.zeros((len(selected_indices),N))
# ------------------------------------------------------------------------
# inverse position of G16
current_stream = seismic_data_streams[selected_indices[-1]]
print(current_stream[0])
raw_matrix[0,:] = current_stream[0].data


# ------------------------------------------------------------------------

for k,idx in enumerate(selected_indices[:-1]):
    current_stream = seismic_data_streams[idx]
    print(current_stream[0])
    raw_matrix[k,:] = current_stream[0].data
    
# compute FFT of all signals 
# FFT_matrix = np.fft.fft(raw_matrix, norm = 'forward')
# FFT_matrix = FFT_matrix[:,:N//2]
# freq = np.fft.fftfreq(N,1/fs)
# freq = freq[:N//2]

# fig, ax = plt.subplots()
# for k in range(np.shape(FFT_matrix)[0]):
#     ax.loglog(freq,abs(FFT_matrix[k,:]))
#     ax.set_xlabel(r'$f \: \rm (Hz)$')
#     ax.set_ylabel(r'$\hat{V}(f)$')
# ax.grid()


#%% Filter Amundsen noise for all geophones
                      
# fc = [2.3, 8.33, 16.66, 25, 41.66, 2*16.66, 3*16.66] #0909
#fc = [2.955, 8.33, 16.66, 25, 41.66, 2*16.66, 3*16.66] #1009


# Load frequencies to cut 
file2load = path2data + 'Frequencies_Amundsen_' + year + '_' + date + '_channel_Z' + '.pkl'
with open(file2load,'rb') as pfile:
    fc = pickle.load(pfile)


fs = 1000
Q = 30.0  # Quality factor
filtered_matrix = raw_matrix
for f0 in fc:

    # Design notch filter
    b, a = iirnotch(f0, Q, fs)
    filtered_matrix = filtfilt(b,a,filtered_matrix)

#%%
t_start = 141
t_end = t_start + 0.8

t = first_trace.times()
mask = [True if ((t_value > t_start) and (t_value < t_end)) else False for t_value in t]

fig, ax = plt.subplots(figsize = fig_size)
for k in range(np.shape(filtered_matrix)[0]):
    current_stream = filtered_matrix[k,mask]
    # norm_data = (current_stream - np.min(current_stream)) / (np.max(current_stream) - np.min(current_stream))
    norm_data = current_stream / np.max(np.abs(current_stream))
    ax.plot(t[mask] - t_start,norm_data + k * geophones_spacing)

ax.set_xlabel(r'$t \: \rm (s)$')                      
ax.set_ylabel(r'$x \: \rm (m)$')                          


file2save = fig_folder + 'zoomed_channel_Z_filtered_dir1'
plt.savefig(file2save + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(file2save + '.png', dpi = img_quality, bbox_inches = 'tight')

# fig, ax = plt.subplots()
# for k in range(np.shape(raw_matrix)[0]):
#     current_stream = raw_matrix[k,:]
#     # norm_data = (current_stream - np.min(current_stream)) / (np.max(current_stream) - np.min(current_stream))
#     norm_data = current_stream / np.max(current_stream)
#     ax.plot(datetime_values,norm_data + k * geophones_spacing)

# ax.set_xlabel(r'$t \: \rm (s)$')                      
# ax.set_ylabel(r'$x \: \rm (m)$')  


#%%               
file2save = path2data + 'Frequencies_Amundsen_' + year + '_' + date + '_channel_' + channel_correspondance[channel] + '.pkl'
with open(file2save,'wb') as pfile:
    pickle.dump(fc,pfile)
                      
#%% Set t0 dictionnary 

# open .pkl file and load dictionnary 
filename = 't1_to_time_' + date + '_' + year + '_short_length.pkl'
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

# S101, S102, S103
key = 'd' + date + 'a' + acqu_numb + 'tS' + '101' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:51:14.23")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:52:37.08")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:54:01.85")

# S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:57:00.82")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:58:03.75")
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = UTCDateTime("2024-09-10T18:59:19.40")

#%% Save t0 dictionnary in pickle file 

with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)

print('Time dictionnary saved')





#%% Load t0 dictionnary 

signal_length = 0.3 # duration in seconds
select_filtered = 1 # boolean, 1 if we want to select filtered signal 

# load data of intial times 
composante = 'E'
channel = 0
# ch = channel_dic[channel]
flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave
direction = 1 # 1 ou 2 

#-------------------------------------------------------------------------------
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
 
base = path2data
pkl_path = base + 't1_to_time_' + date + '_' + year  + '_short_length.pkl'
with open(pkl_path, 'rb') as f:
    loaded_data = pickle.load(f)


# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)
t1_values = [t1,t2,t3]

ta = t1 + signal_length 
num_samples = int(signal_length * first_trace.stats.sampling_rate) # Number of points generated in time_vector 
time_vector = np.linspace(t1.timestamp, ta.timestamp, num_samples) # timestamp gets the number of sec in ta, t1

num_traces = np.shape(filtered_matrix)[0]

# Dynamically determine the length of the third dimension
third_dim_len = len(t1_values) if len(t1_values) > 1 else 1
seismic_matrix = np.zeros((num_traces, third_dim_len, num_samples)) # matrix dim : geo_indices x nb_sources x time

if select_filtered:
    signal_matrix = filtered_matrix
    
else:
    signal_matrix = raw_matrix
    
for i, stream_index in enumerate(selected_indices):
    stream = signal_matrix[i,:]
    # print(stream[0])
    for t1_index, t1 in enumerate(t1_values):
            start_sample = int((t1 - first_trace.stats.starttime) * fs)
            # start_sample = int((t1 - trace.stats.starttime.timestamp) * trace.stats.sampling_rate)
            end_sample = start_sample + num_samples
            seismic_matrix[i, t1_index, :] = signal_matrix[i,start_sample:end_sample]

print('Matrix computed for channel : ' + composante)


# Plot the data
fig, ax = plt.subplots()
#print(seismic_matrix.shape[0]) # Seismic_matrix 1000 values [2], 3 sources [1], 16 geophones [0]
for k in range(seismic_matrix.shape[0]): # size of the array(nb of rows) = nb of geophones
    for t1_index in range(seismic_matrix.shape[1]): #size of the array(nb of collumns) = nb of sources
        ax.plot(time_vector, seismic_matrix[k, t1_index, :] / max(np.abs(seismic_matrix[k, t1_index, :])) + geophones_spacing*k,
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
geophones_spacing = 4 # in meters
signals = np.transpose(seismic_matrix, (0, 2, 1))

f, k, FK = fn_svd(signals, fs, geophones_spacing , rang ,'ExampleName', 0, 'threshold',-40) #plot valeurs  singuliere/freq
F, K = np.meshgrid(f, k)

# Normalize FK based on direction
if direction == 1:
    FK_normalized = FK / np.max(FK)
elif direction == 2:
    FK_normalized = np.flipud(FK / np.max(FK))

# Unwrapping Spectrum
nb_stacking = 0 # number of times we want to stack the FK plot horizontally
idx_stacking = 1
FK_uwp = FK_normalized
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
ymax = 250  # Maximum frequency
xmin = 0  # Minimum wavenumber
xmax = 4  # Maximum wavenumber

if channel == 0:
    wave_type = 'SH0'
elif channel == 1:
    wave_type = 'QS0'
elif channel == 2:
    wave_type = 'QS'

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
nb_points = 10 # points to select on graph FK flexural
#-----------------------------------------------

figname = fig_folder + 'FKplot_' + wave_type + '_' + year + '_' + date + '_acqu_' + acqu_numb + '_dir'  + str(direction)

plt.show()
if flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
    #ax1.set_ylim([-1.5, 1.5])
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    # ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label= r'$|\hat{s}|(f,k)$')

    ax1.set_ylim([ymin, ymax])
    
    plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')
    print('Plot saved')
    
    #-----------------------------------------------------------------------------------------------
    
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
    
    # save figure 
    figname = figname + '_selected_points'
    plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')
    
    
    # save FK matrix 
    FK_dic = {}
    FK_dic['FK'] = FK_uwp
    FK_dic['f'] = f
    FK_dic['k'] = k_uwp
    
    file2save = path2data + 'FK_matrix_QS_' + year + '_' + date + '_acq'+acqu_numb+ '_dir' + str(direction) + '_filtered.pkl'
    with open(file2save,'wb') as pfile:
        pickle.dump(FK_dic,pfile)
    

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
    plt.colorbar(c1, ax=ax1, label= r'$|\hat{s}|(f,k)$')
    ax1.tick_params(axis='both', which='major', pad=7)
    ax1.set_ylim([ymin, ymax])
    
    
    plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')
    print('Plot saved')
    
    
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

    figname = figname + '_selected_points'
    plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')

#%% Save data in a dictionnary  

pkl_file = path2data + 'Phase_velocity_dictionnary_acqu_' + acqu_numb + '_filtered.pkl'
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

#%% Save phase velocity data 

with open (pkl_file,'wb') as pfile:
    pickle.dump(s,pfile)
    
#%% Load phase velocity data
pkl_file = path2data + 'Phase_velocity_dictionnary_acqu_' + acqu_numb + '_filtered.pkl'
with open (pkl_file,'rb') as pfile:
    s = pickle.load(pfile)    


#%% Compute horizontal modes phase velocity in both directions

C_QS0 = 0.5*(s['dir1']['C_longi'] + s['dir2']['C_longi'])
file2save = path2data  +  year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir_filtered.pkl'
with open(file2save, 'wb') as file:
    pickle.dump(C_QS0, file)
    
C_SH0 = 0.5*(s['dir1']['C_shear'] + s['dir2']['C_shear'])
file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir_filtered.pkl'
with open(file2save, 'wb') as file:
    pickle.dump(C_SH0, file)

#%% Save (f,k) points associated to QS mode 
file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' + str(direction) +'_filtered.pkl'
with open(file2save, 'wb') as file:
     pickle.dump([f_mode, k_mode], file)
    
#######################################################
#%% Load data for thickness inversion 
#######################################################

# load phase velocity data 
pkl_file = path2data + 'Phase_velocity_dictionnary_acqu_' + acqu_numb + '_filtered.pkl'
with open (pkl_file,'rb') as pfile:
    s = pickle.load(pfile)    

# load (f,k) points from QS mode 

file2save = path2data + year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' + str(direction) +'_filtered.pkl'
with open(file2save, 'rb') as pfile:
     s_QS = pickle.load(pfile) # dim 1 = f, dim 2 =k 
     
# load FKuwp matrix 
file2save = path2data + 'FK_matrix_QS_' + year + '_' + date + '_acq'+acqu_numb+ '_dir' + str(direction) + '_filtered.pkl'
with open(file2save,'rb') as pfile:
    s_FK = pickle.load(pfile)

#%% Compute E , nu 
rho_ice = 917 # kg/m3
dic_phase = s['dir' + str(direction)]
s_Enu = E_nu_inversion(rho_ice, dic_phase['C_longi'], dic_phase['C_shear'])

print(s_Enu)

s_Enu = E_nu_inversion(rho_ice, C_QS0, C_SH0)

print(s_Enu)

#%% --------- INVERT ICE THICKNESS -------------------

rho_ice = 917 
E_fit = s_Enu['E'] # 2.43e9
nu_fit = s_Enu['nu']

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(1.0,5.0,h_precision) # array of height tested 
#path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones' # path 2 data points saved from FK flexural
# acqu_numb = '0002' # acquisition number 

fig, ax1 = plt.subplots(1,1,figsize = fig_size)
vmin = 0     # Minimum value for colormap
vmax = 1 # Maximum value for colormap (adjust as needed)

# Computes extents for imshow
x = extents(s_FK['k'])
y = extents(s_FK['f'])

c1 = ax1.imshow(s_FK['FK'].T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

#ax1.set_ylim([-1.5, 1.5])
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
# ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label= r'$|\hat{s}|(f,k)$')
ax1.tick_params(axis='both', which='major', pad=7)
ax1.set_ylim([0, 200])

# add detected points
ax1.scatter(s_QS[1], s_QS[0], color='r', label='Detected points')    
    

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h[i], E_fit, nu_fit,s_QS[0],c_w,rho_w)
    error = np.sum(np.power((s_QS[1]-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
   
# computes wavevectors k for the best fit
k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h_ice, E_fit, nu_fit,s_QS[0],c_w,rho_w)

ax1.plot(k_mode_synthetic, s_QS[0], color='g', label='Best fit')    
ax1.legend()


figname = fig_folder + 'FK_mode_QS_acqu' + acqu_numb + '_dir' + str(direction) +  '_filtered_dir2' 
plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')
print(h_ice)     


fig, ax = plt.subplots(1,1,figsize = (16,9))

x = extents(s_FK['k'])
y = extents(s_FK['f'])

c1 = ax.imshow(s_FK['FK'].T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)
#ax1.set_ylim([-1.5, 1.5])
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
# ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax, label= r'$|\hat{s}|(f,k)$')
ax.tick_params(axis='both', which='major', pad=7)
ax.set_ylim([0, 200])

ax.plot(k_mode_synthetic, s_QS[0], '-', color = 'g', linewidth = 3)

figname = fig_folder + 'FK_mode_QS_acqu' + acqu_numb + '_dir' + str(direction) + '_filtered_dir2_with_theory' 
plt.savefig(figname + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = img_quality, bbox_inches = 'tight')

#%% Save data used for inversion 
data = {}
data['E'] = E_fit
data['rho_ice'] = rho_ice
data['nu_fit'] = nu_fit
data['c_w'] = c_w
data['rho_w'] = rho_w
data['h_ice'] = h_ice

pkl_file_datafit = path2data + 'Data_used_for_thickness_inversion_acqu' + acqu_numb + '_dir' + str(direction) + '_filtered'
with open(pkl_file_datafit,'wb') as pfile :
    pickle.dump(data,pfile)
print('Data saved')
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     


