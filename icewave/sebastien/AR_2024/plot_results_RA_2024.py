# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 15:18:10 2025

@author: sebas
"""

import os
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime, time 
import numpy as np
from scipy.signal import iirnotch, filtfilt
from scipy.interpolate import make_interp_spline
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import pickle
import csv 
import math

import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

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


##############################################################
#%% --------- SETTING PARAMETERS AND PATHS -------------------
##############################################################
date = '0909'
year = '2024'
main_path = 'F:/Amundsen_RA_2024/Data/'
path2data = f'{main_path}{year}/{date}/Geophones/'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'

#----------------------- SELECTING ACQ NUMBER AND CHANNEL ------------------ 

acqu_numb = '0001' 
channel_correspondance = ['E','N','Z']
gain = 12 # gain in dB 
scale_factor = 10**(-gain/10)
geophones_spacing = 4 # space between geophones, in meter 

fig_folder = path2data + 'Figures_Rapport/'  # folder where figures are saved 
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
    
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

#%% Build a matrix of raw signals for a given channel

channel = 0 # 0 = 'E', 1 = 'N', 2 = 'Z'
composante = channel_correspondance[channel]
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
    
#%% Load t0 dictionnary 

signal_length = 0.3 # duration in seconds
select_filtered = 1 # boolean, 1 if we want to select filtered signal 

# load data of initial times 
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
             f"Channel: {channel}, Source: {composante}, Direction: {direction}")
ax.set_ylabel('Normalized Seismic Data')
ax.tick_params(axis='x')
ax.tick_params(axis='y')
plt.show()


#%%  Calling fn_SVD

#####################################
##### IMPORTANT PARAMETERS ##########
#####################################

rang = [0,1,2]
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

set_graphs.set_matplotlib_param('single')
fig, ax1 = plt.subplots()

#-----------------------------------------------
# Set parameters 
threshold = 0.2 # minimum relative amplitude to detect a local maximum in the dispersion relation FK plot 
precision_k = [0.03,0.025] # precision over k values (when searching for a local maximum)
prec = precision_k[flexure_wave]
# 0: horizontal_wave 1: flexure_wave
semirange_freq = 5 # semi-range of frequencies for horizontal waves
nb_points = 10 # points to select on graph FK flexural
#-----------------------------------------------

figname = f'{fig_folder}FKplot_{wave_type}_{year}_{date}_acqu_{acqu_numb}_dir{str(direction)}'

plt.show()
if flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap=parula_map,
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)

    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(c1,cax = cax)
    cbar.set_label(r'$|\hat{s}|$')

    ax1.set_ylim([ymin, ymax])
    
    plt.savefig(figname + '.pdf', bbox_inches = 'tight')
    plt.savefig(figname + '.png', bbox_inches = 'tight')
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
    figname = f'{figname}_selected_points'
    plt.savefig(figname + '.pdf', bbox_inches = 'tight')
    plt.savefig(figname + '.png', bbox_inches = 'tight')
    
    # save FK matrix 
    FK_dic = {}
    FK_dic['FK'] = FK_uwp
    FK_dic['f'] = f
    FK_dic['k'] = k_uwp
    
    # file2save = f'{path2data}FK_matrix_QS_{year}_{date}_acq{acqu_numb}_dir{str(direction)}_filtered.pkl'
    # with open(file2save,'wb') as pfile:
    #     pickle.dump(FK_dic,pfile)
    

#--------------------------------------------------------------------------------
elif horizontal_wave:

    # Computes extents for imshow
    x = extents(k)
    y = extents(f)

    c1 = ax1.imshow(FK_normalized.T, aspect = 'auto', cmap=parula_map,origin = 'lower',
                    extent = x + y,vmin = vmin, vmax = vmax)

    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(c1,cax = cax)
    cbar.set_label(r'$|\hat{s}|$')
    ax1.tick_params(axis='both', which='major', pad=7)
    ax1.set_ylim([ymin, ymax])
    
    
    plt.savefig(figname + '.pdf', bbox_inches = 'tight')
    plt.savefig(figname + '.png', bbox_inches = 'tight')
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

    figname = f'{figname}_selected_points'
    plt.savefig(figname + '.pdf', bbox_inches = 'tight')
    plt.savefig(figname + '.png', bbox_inches = 'tight')

#######################################################
#%% Load data for thickness inversion 
#######################################################

# load phase velocity data 
pkl_file = f'{path2data}Phase_velocity_dictionnary_acqu_{acqu_numb}_filtered.pkl'
with open (pkl_file,'rb') as pfile:
    s = pickle.load(pfile)    

# load (f,k) points from QS mode 

file2save = f'{path2data}{year}_{date}_acq{acqu_numb}disp_QS_dir{str(direction)}_filtered.pkl'
with open(file2save, 'rb') as pfile:
     s_QS = pickle.load(pfile) # dim 1 = f, dim 2 =k 
     
# load FKuwp matrix 
file2save = f'{path2data}FK_matrix_QS_{year}_{date}_acq{acqu_numb}_dir{str(direction)}_filtered.pkl'
with open(file2save,'rb') as pfile:
    s_FK = pickle.load(pfile)


#%% Compute E , nu 
rho_ice = 917 # kg/m3


# Compute mean value of E and nu over both directions
C_QS0 = 0.5*(s['dir1']['C_longi'] + s['dir2']['C_longi'])
C_SH0 = 0.5*(s['dir1']['C_shear'] + s['dir2']['C_shear'])
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

set_graphs.set_matplotlib_param('single')
fig, ax1 = plt.subplots()
vmin = 0     # Minimum value for colormap
vmax = 1 # Maximum value for colormap (adjust as needed)

# Computes extents for imshow
x = extents(s_FK['k'])
y = extents(s_FK['f'])

c1 = ax1.imshow(s_FK['FK'].T, aspect = 'auto', cmap=parula_map,
                origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(c1,cax = cax)
cbar.set_label(r'$|\hat{s}|$')
ax1.tick_params(axis='both', which='major', pad=7)
ax1.set_ylim([0, 200])

# add detected points
ax1.scatter(s_QS[1], s_QS[0], color='r', label='Detected points')    
    

# find h_ice that minimizes distance to data points 
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


figname = f'{fig_folder}FK_mode_QS_acqu{acqu_numb}_dir{str(direction)}_filtered_dir2' 
plt.savefig(figname + '.pdf', bbox_inches = 'tight')
plt.savefig(figname + '.png', bbox_inches = 'tight')
print(h_ice)     

#%% Create clean plot with fitted data and theory
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

x = extents(s_FK['k'])
y = extents(s_FK['f'])

c1 = ax.imshow(s_FK['FK'].T, aspect = 'auto', cmap=parula_map,
               origin = 'lower',extent = x + y,vmin = 0, vmax = vmax)

ax.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(c1,cax = cax)
cbar.set_label(r'$|\hat{s}|$')
ax.tick_params(axis='both', which='major', pad=7)
ax.set_ylim([0, 100])

# add detected points
ax.plot(s_QS[1], s_QS[0], '.', color='r',markersize = 8,zorder = 2)    
# add theory
f_th = np.linspace(2,60,500)
k_th, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h_ice, E_fit, nu_fit,f_th,c_w,rho_w)

ax.plot(k_th, f_th, '-', color = 'w', linewidth = 3,zorder = 1)

figname = f'{fig_folder}FK_mode_QS_acqu{acqu_numb}_dir{str(direction)}_filtered_dir2_with_theory' 
plt.savefig(figname + '.pdf', bbox_inches = 'tight')
plt.savefig(figname + '.png', bbox_inches = 'tight')

#%% Save data used for inversion 
data = {}
data['E'] = E_fit
data['rho_ice'] = rho_ice
data['nu_fit'] = nu_fit
data['c_w'] = c_w
data['rho_w'] = rho_w
data['h_ice'] = h_ice

pkl_file_datafit = f'{path2data}Data_used_for_thickness_inversion_acqu{acqu_numb}_dir{str(direction)}_filtered.h5'
rw.save_dict_to_h5(data,pkl_file_datafit)
print('Data saved')