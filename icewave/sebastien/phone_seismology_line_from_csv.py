# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 01:21:45 2025

@author: sebas
"""

#%% Import modules 
import numpy as np
import math
from datetime import datetime
import pandas as pd

import matplotlib.pyplot as plt
# from obspy.core import UTCDateTime
# from obspy import read
import shutil
import os
import re
import glob

import pickle
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import warnings
from matplotlib.patches import PathPatch
from scipy.interpolate import make_interp_spline

import time
import csv 
import h5py

import icewave.tools.matlab2python as mat2py
import icewave.tools.Fourier_tools as FT

plt.rcParams.update({
    "text.usetex": True}) # use latex

phone_spacing  = 3 # space between phones in meters 
acqu_numb = '1' # acquisition number 
fs = 400 # acquisition frequency 

#%%
def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#--------------------------------------------------------------------------------------------------

def svd_phones(seismic_matrix,fs,xs,Nreceiv,Nemit,rang,*varargin):
    """ Computes spatio-temporal Fourier Transform using SVD decomposition.
    Arguments : 
        - signals : seismic_matrix, python dictionnary containing phone keys, and for each phone key a dictionnary 
        containing a set of nb_sources 1D signals 
        - fs : sampling frequency in time 
        - xs : distance between cellphones (inverse of spatial sampling)
        - rang : row used to perform decomposition using SVD method 
        (number of singulavr values used to perform decomposition) 
        - varargin : optionnal argument
        
    Returns : 
        
        - f : array of frequencies, scaled
        - k : array of wavenumber, scaled 
        - FK : amplitude of FFT of the signal in time and space, using SVD method 
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

    # Time domain fft
    Nf = 1024
    Nk = 1024
    f = (np.arange(Nf) / Nf) * (fs if fs else 1)
    f_axename = 'f/fs' if not fs else 'f'
    SIGNALS = np.zeros((Nreceiv,Nf,Nemit),dtype = complex)
    for i,phone_key in enumerate(seismic_matrix.keys()):
        for j, time_key in enumerate(seismic_matrix[phone_key].keys()):
            current_fft = fft(seismic_matrix[phone_key][time_key] - np.mean(seismic_matrix[phone_key][time_key]),
                              n = Nf)

            SIGNALS[i,:,j] = current_fft

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
        #plt.plot(f, titi, label=f'Slice {ne}')
        ax.plot(f[1:-1],titi[1:-1],label = f'Slice {ne}')
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('Singular values (in dB of peak value)')

    if threshold_user is None: 
        nb_points = 5
        print(f'Select {nb_points} on the graph to define noise minimum')
        [fcut, sigmacutsup] = plt.ginput(nb_points)
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

    # projection onto each singular vector

    k = (np.arange(Nk) / Nk) * (2 * np.pi / xs)  if xs else np.arange(Nk + 1)
    k_axename = 'k/ks' if not xs else 'k'

    projections = ifft(U, Nk, axis=0)#np.fft.fftshift(ifft(U, Nk, axis=0), axes=0)
    projections_sum = np.zeros((Nf, Nk, Nemit))

    for kemit in rang:
        for ii in range(Nf):
            max_value = 1  # np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2

    FK = np.abs(np.mean(projections_sum, axis=2))
    
    return f,k,FK


#-------------------------------------------------------------------------------------------------------------

def compute_Enu(C_longi,C_shear,rho_ice):
    
    nu = 1-2*(C_shear/C_longi)**2
    E = rho_ice*C_longi**2*(1-nu**2)
    return E,nu

#-------------------------------------------------------------------------------------------------------------

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
    "text.usetex": True}) # use latex
    
#%% Load data using Victor's code 
path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Donnees_phones_Victor/'
dossier_csv = path2data

acqu_numb = '1'
channel = 'a_Z'

data = {}
# Parcours des fichiers du dossier
keys_channel = ['ts','a_N', 'a_E', 'a_Z',  'colonne5']
print("Chargement des données")

for fichier in os.listdir(dossier_csv):

    if fichier.endswith(".csv"):  # Vérifie que c'est un fichier CSV
        chemin_fichier = os.path.join(dossier_csv, fichier)

        cle = fichier.split('_')[0]
        
        if len(cle) == 2:
            cle = 'P0' + cle[-1]
        else:
            cle = 'P' + cle[-1]
        print(cle)
        # Chargement du fichier CSV en DataFrame
        df = pd.read_csv(chemin_fichier)
        # Stockage des colonnes dans un dictionnaire
        data[cle] = {keys_channel[i]: df.iloc[:, i].tolist() for i in range(df.shape[1])}


print("Chargement terminé !")


#%%
synchro_file = f'{path2data}tsync_2025-04-15 05_20_36.txt'
# file_path = directory  + file_name
# data_sync_tab = extract_data(file_path)

data_sync = {}
for ligne in data_sync_tab:

    cle = ligne[0].split('_')[0]
    if len(cle) == 1:
        cle = 'P0' + cle
    else:
        cle = 'P' + cle
    valeurs = ligne[1]  # Le reste des colonnes comme valeurs (liste)
    
    # Ajouter au dictionnaire
    data_sync[cle] = valeurs




data = dict(sorted(data.items()))
for cle in data.keys():
    data[cle]['ts'] = np.array(data[cle]['ts'], dtype = float)*10**-6  + np.array(data_sync[cle], dtype = float)
    print(data[cle]['ts'])
for i,phone_key in enumerate(data.keys()):
    data[phone_key][channel] = np.array(data[phone_key][channel], dtype = float)
    data[phone_key]['ts'] = np.array(data[phone_key]['ts'], dtype = float)
    data[phone_key]['ts'] = data[phone_key]['ts'][~np.isnan(data[phone_key][channel])]
    data[phone_key][channel] = data[phone_key][channel][~np.isnan(data[phone_key][channel])]

    print(np.shape(data[phone_key]['ts']))
    print(np.shape(data[phone_key][channel]))
























#%% Load data
 
date = '0211'
base  = f'F:/Rimouski_2025/Data/{date}/Phone/Summary/'

filelist = glob.glob(base + '*.h5')
file2load = filelist[0]

with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
    
    a_group_key = list(f.keys())
    data_unsort = {}
    for i,key in enumerate(a_group_key): 
        data_unsort[key] = mat2py.mat_to_dict(f[key], f[key])

# sort data
sorted_keys = np.sort(np.asarray(list(data_unsort.keys())).astype(int)).astype(str)
data = {}
for key in sorted_keys :
    data[key] = data_unsort[key]
    
fig_folder = base + 'Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
#%% Create data matrix 

keys_channel = ['a_E','a_N','a_Z','t_s']

# compute minimal time between which we build a seismic matrix 
# Find the largest start time (tstart) and the smallest end time (tend) among all streams
tstart = max([data[key]['ts'][0] for key in data.keys()])
tend = min([data[key]['ts'][-1] for key in data.keys()])

# truncate all streams between tstart and tend  
for phone_key in data.keys():
    relevant_signal = np.logical_and(data[phone_key]['ts'] > tstart,data[phone_key]['ts'] < tend)
    for key in data[phone_key].keys():
        data[phone_key][key] = data[phone_key][key][relevant_signal]
        
#%% Apply time shift to some phones 

# set dictionnaryof time shift in s
time_shift = {}
for key_phone in data.keys():
    time_shift[key_phone] = 0

time_shift['4'] = 0.015
time_shift['15'] = 0.010

for key_phone in data.keys():
    data[key_phone]['ts'] = data[key_phone]['ts'] + time_shift[key_phone]

#%% Plot selected channel for all phones 

channel = 'a_Z'

fig, ax = plt.subplots()
for i,phone_key in enumerate(data.keys()):
    current_stream = data[phone_key][channel] - np.mean(data[phone_key][channel])

    ax.plot(data[phone_key]['ts'],current_stream/max(current_stream) + phone_spacing*(i+1))
    
ax.set_title(f'Signal superposition for channel {channel}')

#%% Select time range of processing 

filename = 't0_time_' + date  + '.pkl'
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
time_dict[key] = 49247.0
key = 'd' + date + 'a' + acqu_numb + 'tS' + '102' + composante 
time_dict[key] = 49192.1
# key = 'd' + date + 'a' + acqu_numb + 'tS' + '103' + composante 
# time_dict[key] = ''

# # S104, S105, S106
key = 'd' + date + 'a' + acqu_numb + 'tS' + '104' + composante 
time_dict[key] = 49394.1
key = 'd' + date + 'a' + acqu_numb + 'tS' + '105' + composante 
time_dict[key] = 49456.4
key = 'd' + date + 'a' + acqu_numb + 'tS' + '106' + composante 
time_dict[key] = 49508.7


# #S107
key = 'd' + date + 'a' + acqu_numb + 'tS' + '107' + composante 
time_dict[key] = 49327.85

with open(file2save, 'wb') as f:
    pickle.dump(time_dict, f)

print('Time dictionnary saved')

#%% Load time table

path2table = 'C:/Users/sebas/OneDrive/Bureau/'
path2file = f'{path2table}Datat0_time_0211.pkl'

with open(path2file,'rb') as pf:
    loaded_data2 = pickle.load(pf)




#%% Create a seismic matrix 

signal_length = 1 # duration in seconds 
composante = 'N'
channel = 'a_N'
direction = 2 # 1 ou 2

flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave

pkl_path = base + 't0_time_' + date + '.pkl'
with open(pkl_path, 'rb') as pfile:
    loaded_data = pickle.load(pfile)

# assign a string to S values depending on the direction
if direction == 1 :
    S1 = '101' 
    S2 = '102' 
    S3 = '103'
if direction == 2: 
    S1 = '104' 
    S2 = '105' 
    S3 = '106'
    
# get the 3 differents t0 values for the required source
t1 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + date + 'a' + acqu_numb + 'tS' + S3 + composante)
t0_values = [t1,t2,t3]


# compute minimal length of signals 
minimal_length = min([len(data[key]['ts']) for key in data.keys()])
maximal_length = max([len(data[key]['ts']) for key in data.keys()])

nb_receptors = len(data) 
time_steps = signal_length*fs
nb_sources = len(t0_values) if len(t0_values) > 1 else 1
seismic_matrix = {}

# fill seismic_matrix 
for i,phone_key in enumerate(data.keys()):
    seismic_matrix[phone_key] = {}
    for j,t0 in enumerate(t0_values) :
        # start_sample = int((t0 - data[phone_key]['ts'][0])*fs)
        # end_sample = start_sample + time_steps
        # print(f'start_sample : {start_sample} / end_sample : {end_sample}')
        # seismic_matrix[phone_key][:,j] = data[phone_key][channel][start_sample:end_sample]
        
        
        key_time = f't{j+1}'
        mask_sample = np.logical_and(data[phone_key]['ts'] > t0, data[phone_key]['ts'] < t0 + signal_length)
        seismic_matrix[phone_key][key_time] = data[phone_key][channel][mask_sample]
        
        
# Plot selected signals
fig, ax = plt.subplots()
for i,phone_key in enumerate(seismic_matrix.keys()):
    # for j in range(np.shape(seismic_matrix[phone_key])[1]):
    for j,key_time in enumerate(seismic_matrix[phone_key].keys()):
        current_stream = seismic_matrix[phone_key][key_time] - np.mean(seismic_matrix[phone_key][key_time])

        ax.plot(current_stream/max(current_stream) + phone_spacing*(i+1))
    
ax.set_title(f'Selected samples for channel : {channel}')

#%% Perform Fourier transform

Nreceiv = nb_receptors
Nemit = nb_sources
xs = phone_spacing
Nf = 1024
Nk = 512

f = fs*np.arange(0,(Nf//2))/Nf
SIGNALS = np.zeros((Nreceiv,Nf//2,Nemit),dtype = complex)

for i,phone_key in enumerate(seismic_matrix.keys()):
    for j, time_key in enumerate(seismic_matrix[phone_key].keys()):
        current_fft,f = FT.fft_1D(seismic_matrix[phone_key][time_key], Nf, fs)
        SIGNALS[i,:,j] = current_fft
    
# Fourier transform in space 

k = 2*np.pi*np.arange(0,Nk//2)/Nk/xs
FK = fft(SIGNALS,n = Nk,axis = 0) 
FK = FK[:Nk//2,:,:]/Nreceiv
FK[1:-1,:,:] = 2*FK[1:-1,:,:]

# Show Fourier transform 
i = 0
fig, ax = plt.subplots()
sub_FK = FK[:,:,i]

ax.imshow(abs(sub_FK.T),origin = 'lower', aspect = 'auto',norm = 'linear',extent = extents(k) + extents(f))


#%% Compute SVD 

rang = [0,1,2]
f,k,FK = svd_phones(seismic_matrix, fs, phone_spacing, nb_receptors, nb_sources, rang, 'threshold', -7)
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

#%% Detect points

vmin = 0     # Minimum value for colormap
vmax = 1     # Maximum value for colormap (adjust as needed)
ymin = 0  # Minimum frequency
ymax = 400  # Maximum frequency
xmin = 0  # Minimum wavenumber
xmax = 4  # Maximum wavenumber

fig, ax1 = plt.subplots(1, 1, figsize=fig_size)

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
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
    ax1.set_title('Spectrum with SVD filter')
    plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')

    ax1.set_ylim([0, 100])
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
    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}|_{max}}(f,k)$')
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
vmin = 0     # Minimum value for colormap
vmax = 1     # Maximum value for colormap (adjust as needed)
ymin = 0  # Minimum frequency
ymax = 200  # Maximum frequency
xmin = 0  # Minimum wavenumber
xmax = 4  # Maximum wavenumber

fig, ax1 = plt.subplots(1, 1, figsize=fig_size)

if horizontal_wave:
    # Computes extents for imshow
    x = extents(k)
    y = extents(f)

    c1 = ax1.imshow(FK_normalized.T, aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,
                    vmin = vmin, vmax = vmax)
    # plt.plot(k_values_list, poly_function(k_values_list), color='white', linestyle='--', linewidth = 3, 
    #          label='Linear Regression')

    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)

    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}(f,k)$')
    ax1.tick_params(axis='both', which='major', pad=7)
    ax1.set_ylim([ymin, ymax])
    
    if channel == 'a_E':
        wave_type = 'SH0'
    elif channel == 'a_N':
        wave_type = 'QS0'
        
    figname = fig_folder + 'Horizontal_FKplot_' + wave_type + '_acqu_' + acqu_numb + '_dir'  + str(direction) + '_sig_length_' + str(signal_length).replace('.','p')

    plt.savefig(figname + '.pdf', dpi = 200, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = 200, bbox_inches = 'tight')
    
    
elif flexure_wave:
    c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                    origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
    # plt.plot(k_mode, f_mode, color='white', linestyle='-', linewidth = 3, label='Selected points')

    ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
    ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')

    plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}(f,k)$')
    ax1.set_ylim([ymin, ymax])
    
    wave_type = 'QS'
    figname = fig_folder + 'Flexural_FKplot_' + wave_type + '_acqu_' + acqu_numb + '_dir'  + str(direction)

    plt.savefig(figname + '.pdf', dpi = 200, bbox_inches = 'tight')
    plt.savefig(figname + '.png', dpi = 200, bbox_inches = 'tight')
    

#%% Save data in a dictionnary  

pkl_file = base + 'Phase_velocity_dictionnary.pkl'
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

with open (pkl_file,'wb') as pfile:
    pickle.dump(s,pfile)
    
#%% Save (f,k) points associated to QS mode 
file2save = f'{base}{date}_disp_QS_dir_{direction}.pkl'
with open(file2save, 'wb') as file:
     pickle.dump([f_mode, k_mode], file)
     
     
#%% Compute Young modulus and Poisson coefficient 

direction = 2
key_dir = 'dir' + str(direction)
shear_velocity = s[key_dir]['C_shear']
acoustic_velocity = s[key_dir]['C_longi']
rho_ice = 917

E,nu = compute_Enu(acoustic_velocity, shear_velocity, rho_ice)

print(f'Young modulus, E = {E*1e-9} and Poisson coefficient, nu = {nu}')


#%% Invert Ice thickness 



#######################################################################
#%% --------- INVERT ICE THICKNESS USING FK MATRIX -------------------
#######################################################################

rho_ice = 917 
E_fit = E # 2.43e9
nu_fit = nu

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 

# load selected points on the FK plot 
file2load = f'{base}{date}_disp_QS_dir_{direction}.pkl'
with open(file2load, "rb") as filename:
    data = pickle.load(filename)
f_mode1 = data[1] 
k_mode1 = data[0] 

fig, ax = plt.subplots()
ax.plot(k_mode1,f_mode1)
# plt.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points')    

f_mode1 = f_mode1[1:]
k_mode1 = k_mode1[1:]

# plt.plot(k_mode1, f_mode1, linestyle='--', color='r', label='Line between points') 
# #plt.plot(f_mode2, k_mode2, linestyle='--', color='r', label='Line between points') 

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