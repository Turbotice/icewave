
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
# fontsize of the figure title


    
#%% Load data
 
nom_mesure = input("nom_mesure: ")
base  = '../Data'


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
t0_values = [t1,t2]


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

plt.show()