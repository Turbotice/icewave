
#%% Import modules 
import numpy as np
import math
from datetime import datetime

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
import pandas as pd 
from victor_func import *
import Fourier_tools as FT


phone_spacing  = 3 # space between phones in meters 
acqu_numb = '1' # acquisition number 
fs = 400 # acquisition frequency 
fig_size = (12,9)
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
        - k : arraumber, scaled 
        - FK : amplitude of FFT y of wavenof the signal in time and space, using SVD method 
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
# fontsize of the figure title


    
#%% Load data
 
nom_mesure = input("nom_mesure: ")
base  = "Data"


data = {}
dossier_csv = base +"/" + nom_mesure
# Parcours des fichiers du dossier
keys_channel = ['ts','a_N', 'a_E', 'a_Z',  'colonne5']
print("Chargement des données")

for fichier in os.listdir(dossier_csv):

    if fichier.endswith(".csv"):  # Vérifie que c'est un fichier CSV
        chemin_fichier = os.path.join(dossier_csv, fichier)

        cle = fichier.split('-')[0]
        # Chargement du fichier CSV en DataFrame
        df = pd.read_csv(chemin_fichier)
        # Stockage des colonnes dans un dictionnaire
        data[cle] = {keys_channel[i]: np.array(df.iloc[:, i].tolist()) for i in range(df.shape[1])}


print("✅ Chargement terminé !")



directory = base + "/Tsync/"  
list_files(directory)
file_name= input("Choisir un fichier de syncro: ")

file_path = directory  + file_name
data_sync_tab = extract_data(file_path)
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

for cle in data.keys():
    data[cle]['ts'] = np.array(data[cle]['ts'], dtype = float)*10**-6  + np.array(data_sync[cle], dtype = float)

# compute minimal time between which we build a seismic matrix 
# Find the largest start time (tstart) and the smallest end time (tend) among all streams
data = dict(sorted(data.items()))
liste_tel = selection_utilisateur(list(data.keys())) 


tstart = max([data[key]['ts'][0] for key in liste_tel])
tend = min([data[key]['ts'][-1] for key in liste_tel])

# truncate all streams between tstart and tend  
for phone_key in liste_tel:
    relevant_signal = np.logical_and(data[phone_key]['ts'] > tstart,data[phone_key]['ts'] < tend)

    for key in data[phone_key].keys():
        data[phone_key][key] = data[phone_key][key][relevant_signal]
        
#%% Apply time shift to some phones 

# set dictionnaryof time shift in s
time_shift = {}
for key_phone in liste_tel:
    time_shift[key_phone] = 0

time_shift['4'] = 0.015
time_shift['15'] = 0.010

for key_phone in liste_tel:
    data[key_phone]['ts'] = data[key_phone]['ts'] + time_shift[key_phone]



#%% Create a seismic matrix 

signal_length = 1 # duration in seconds 
channel = input("Choisir la voie: \na_N \na_E \na_Z\n: ")
if channel == '':
    channel = 'a_Z'
composante = channel[-1]

direction = int(input("direction 1 ou 2:")) # 1 ou 2

flexure_wave = composante == 'Z' # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes
horizontal_wave = not flexure_wave

pkl_path = base + '/pic_time/t0_time_' + nom_mesure + '.pkl'
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

t1 = loaded_data.get('d' + nom_mesure + 'a' + acqu_numb + 'tS' + S1 + composante)
t2 = loaded_data.get('d' + nom_mesure + 'a' + acqu_numb + 'tS' + S2 + composante)
t3 = loaded_data.get('d' + nom_mesure + 'a' + acqu_numb + 'tS' + S3 + composante)
t0_values = [t1,t2,t3]

t0_values = [element for element in  t0_values if element is not None]
print(t0_values)
# compute minimal length of signals 
minimal_length = min([len(data[key]['ts']) for key in data.keys()])
maximal_length = max([len(data[key]['ts']) for key in data.keys()])

nb_receptors = len(liste_tel) 
time_steps = signal_length*fs
nb_sources = len(t0_values) if len(t0_values) > 1 else 1
seismic_matrix = {}

# fill seismic_matrix 
for i,phone_key in enumerate(liste_tel):
    seismic_matrix[phone_key] = {}
    for j,t0 in enumerate(t0_values) :
        # start_sample = int((t0 - data[phone_key]['ts'][0])*fs)
        # end_sample = start_sample + time_steps
        # print(f'start_sample : {start_sample} / end_sample : {end_sample}')
        # seismic_matrix[phone_key][:,j] = data[phone_key][channel][start_sample:end_sample]
        
        key_time = f't{j+1}'
        mask_sample = np.logical_and(data[phone_key]['ts'] > t0, data[phone_key]['ts'] < t0 + signal_length)
        seismic_matrix[phone_key][key_time] = data[phone_key][channel][mask_sample]
        
        
 #%% Compute SVD 

rang = [0,1]
f,k,FK = svd_phones(seismic_matrix, fs, phone_spacing, nb_receptors, nb_sources, rang, 'threshold', -20)
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
