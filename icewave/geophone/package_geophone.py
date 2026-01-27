# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 10:55:55 2025

@author: sebas

Gather all useful functions for geophones data 

"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from obspy.core import UTCDateTime
from obspy import read
import os
import pickle
from scipy.fftpack import fft, ifft
from scipy.linalg import svd

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

#---------------------------------------------------------------------------------------------

def build_data_streams(path2miniseeds,geophones_table_path):
    """ Build seismic data streams matrix, datetime values as well as sampling frequency. 
    The data matrix is ordered as an increasing order of geophone index, with each of the three 
    channels following for a single geophone : [E1,N1,Z1,E2,N2,Z2,...]
    Inputs : - path2miniseeds, path to miniseeds files
             - geophones_table_path, path to geophone correspondance table
    Outputs : - seismic_data_streams, matrix of all geophones streams
              - datetime_values, array of date time 
              - fs, float, sampling frequency
        """
    
    
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
    seismic_data_streams, miniseed_files = read_data(path2miniseeds)

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
    time_vector = convert_to_utc_times(first_trace, first_trace.times())
    datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
    fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 

    
    return seismic_data_streams,datetime_values,fs

#--------------------------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------------------

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
    