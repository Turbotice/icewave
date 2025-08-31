# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 14:09:02 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors 
import matplotlib.cm as cm

import h5py 
import glob
import os 
import pickle

from datetime import datetime
import pytz

import scipy
import pywt

import sys
# sys.path.append('C:/Users/sebas/git')
# sys.path.append('/media/turbots/DATA/thiou/homes/skuchly/git/')

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

global g
g = 9.81


#%% FUNCTION SECTION 

def plot_spatio_temp(spatio,fiber_length,extents,cmap):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsh, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots()
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = cmap,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#----------------------------------------------------------------------------------------------------

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega

#-------------------------------------------------------------------------------------------------------

def extract_peaks(FK,freq,k,freq_range,min_prominence,min_width):
    """ Extract peak coordinates from a space-time spectrum 
    Inputs : - FK, numpy 2D array [nf,nk], space-time spectrum
             - freq, numpy array (nf,), frequencies
             - k, numpy array (nk,), wavevectors
             - freq_range, tuple [freq_min, freq_max], frequency range within we look for peaks
             - min_prominence, minimum prominence of peaks once a row is normalized (see scipy.signal.findpeaks)
             - min_width, minimal width of peaks once a row is normalized (see scipy.signal.findpeaks)
             
    Outputs : - k_exp, array of peaks wavevector coordinate
              - f_exp, array of peaks frequency coordinate """
    
    freq_min = freq_range[0]
    freq_max = freq_range[1]
    idx_min = np.argmin(abs(freq - freq_min))
    idx_max = np.argmin(abs(freq - freq_max))

    indices = np.arange(idx_min,idx_max + 1, step = 1)

    # create arrays for peak detection 
    k_exp = []
    f_exp = []
    
    for idx_freq in indices:
    # print(f'frequency = {freq[idx_freq]:.2f} Hz')
        current_max = np.max(abs(FK[idx_freq,:]))
        # find peaks
        peaks,properties = scipy.signal.find_peaks(abs(FK[idx_freq,:])/current_max,
                                                   prominence = min_prominence,width = min_width)
        if len(peaks) != 0:
            for peak in peaks:
                k_exp.append(k[peak])
                f_exp.append(freq[idx_freq])
    
    k_exp = np.array(k_exp)
    f_exp = np.array(f_exp)

    return k_exp,f_exp


#-------------------------------------------------------------------------------------------------------

def get_water_height(DAS_water_height,UTC_t,selected_x):
    """ Compute water height for a given UTC datetime and position along the fiber
    Input : - DAS_water_height, dictionnary containing recording of water depth evolution during the whole experiment 
            - UTC_t, datetime object, time at which we look for water height
            - selected_x, float, position along optical fiber in meter 
    Output : H, float, water height at the selected time and position 
    """

    closest_idx = np.argmin(abs(UTC_t - DAS_water_height['UTC_t']))
    closest_idx_pos = np.argmin(abs(selected_x - DAS_water_height['s']))
    H = DAS_water_height['water_height'][closest_idx_pos,closest_idx]
    return H

#-------------------------------------------------------------------------------------------------------

def all_chunks_FFT_t(stack_strain,fs):
    """ Perform Fourier Transform in time for all time chunks in a given h5 file. 
    Inputs : - stack_strain, numpy array, dimensions [nb_chunks,nt,nx], contains spatio temporal signals
             - fs, float, sampling frequency 
    Outputs : - FFT_t, numpy array, [nb_chunks,nf,nx] Time Fourier transform of each spatio-temporal chunk 
              - freq, numpy array, dimension (nf,) contains frequencies """

    N = 2**(FT.nextpow2(stack_strain.shape[1]))
    
    mean_sig = np.mean(stack_strain,axis = 1)
    mean_sig = np.tile(mean_sig,(stack_strain.shape[1],1,1)).transpose((1,0,2))
    
    FFT_t = np.fft.fft(stack_strain - mean_sig,n = N,axis = 1)
    FFT_t = FFT_t/stack_strain.shape[1]
    
    # keep only half of frequencies
    FFT_t = FFT_t[:,:N//2,:]
    FFT_t[:,1:-1,:] = 2*FFT_t[:,1:-1,:]
    
    freq = fs*np.arange(0,(N/2))/N
    
    return FFT_t, freq

#-------------------------------------------------------------------------------------------------------

def plot_time_spectrum_chunk(FFT_t,freq,s,UTC_stack,fig_folder,colormap): 
    """ Plot time spectrum for all time chunks 
    Input: - FFT_t, np array, [nb_chunks,nf,nx], Time Fourier transform array
           - freq, np array, (nf,) array of frequencies
           - s, np array, (nx,) array of curvilinear coordinates
           - fig_folder, string, folder where plots will be saved 
           - colormap, colormap object used for plot """

    for chunk in range(FFT_t.shape[0]):
        
        set_graphs.set_matplotlib_param('single')
        fig, ax = plt.subplots()
        imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
                  extent = [freq[0],freq[-1],s[0],s[-1]])
        imsh.set_clim([1e1,1e3])
        
        ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
        ax.set_ylabel(r'$s \; \mathrm{(m)}$')
        ax.set_xlim([0,20])
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh,cax = cax)
        
        title_label = f'UTC : {UTC_stack[chunk,0]}'
        ax.set_title(title_label)
        
        format_date = '%Y-%m-%dT%H-%M-%S'
        label_UTC0 = UTC_stack[chunk,0].strftime(format_date)
        figname = f'{fig_folder}Time_spectrum_UTC0_{label_UTC0}'
        
        plt.savefig(f'{figname}.png',bbox_inches = 'tight')
        plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
        plt.close(fig)
       
#-----------------------------------------------------------------------------------------------------------------

def CWT_all_chunks(FFT_t,scales,wavelet,sampling_period):
    """ Compute Continuous Wavelet Transform over space dimension for each time chunk. See pywt.cwt documentation for more
    precision. This function returns the average over all time chunks of the absolute value of each CWT
    Inputs: - FFT_t, np array, [nb_chunks,nf,nx] Time Fourier Transform for each time chunk and each position
            - scales, np array 1D, scales that will be used during CWT computation 
            - wavelet, string, wavelet used for CWT computation 
            - sampling_period, float, sampling period over space 
    Outputs: - mean_CWT, np array [nf,nk,nx], space-time spectrum averaged over all time chunks
             - k, np array (nk,), array of wavevectors """


    CWT = np.zeros((FFT_t.shape[0],FFT_t.shape[1],len(scales),FFT_t.shape[2]),dtype = 'complex')
    # 0 : chunk_idx, # 1 : freq_idx, # 2 : scales, # 3 : space (curvilinear coordinate)
    
    # loop over all chunks
    for current_chunk in range(FFT_t.shape[0]):
        print(current_chunk)
        current_FFT = FFT_t[current_chunk,:,:]
        
        for i in range(current_FFT.shape[0]): # loop over time frequencies 
            cwtmatr, sigmas = pywt.cwt(current_FFT[i,:], scales, wavelet, sampling_period = sampling_period)
            CWT[current_chunk,i,:,:] = cwtmatr
    
    # Average norm of CWT over all chunks 
    mean_CWT = np.mean(abs(CWT),axis = 0)
    
    # compute wavevector 
    k = 2*np.pi*sigmas

    print('CWT computed successfully !')
    return mean_CWT,k

#----------------------------------------------------------------------------------------------------------------

def extract_D(mean_CWT,selected_x,s,freq,k,peaks_detec,DAS_water_height,UTC_mid):
    """ Compute flexural modulus D from FK plot. First peaks of maximum amplitude are detected
    in the coordinate system (f,k). Then the detected points are fitted by a shallow hydroelastic 
    dispersion relation. Flexural modulus is computed using a least square method fit. 
    Inputs : - mean_CWT, np array [nf,nk,nx], space-time spectrum averaged over all time chunks. COmputed using 
    COntinuous Wavelet Transform. 
             - selected_x, np array 1D, curvilinear coordinates at which D is computed 
             - s, np array 1D (nx,), curvilinear coordinates
             - freq, np array (nf,), frequencies array
             - k, np array (nk,), wavevectors, array
             - peaks_detec, dictionnary, contains properties for peak detection.
             The dictionnary must contain following keys : 
                     + freq_range, tuple, (f_min,f_max), frequency range over which the hydroelastic dispersion relation
                     is fitted 
                     + min_prominence, float, minimal prominence of peaks
                     + min_width, float, minimal width of peaks 
             - DAS_water_height, dictionnary, contains data relative to water depath below the fiber
             - UTC_mid, datetime object, date at which we compute the water height profile below the fiber
    Outputs : - D, np array 1D, values of flexural modulus 
              - err_D, np array 1D, values of standard deviation for each position """


    freq_range = peaks_detec['freq_range']
    min_prominence = peaks_detec['min_prominence']
    min_width = peaks_detec['min_width']
    
    list_idx_x = [np.argmin(abs(s - x)) for x in selected_x]
    
    D = np.zeros(len(selected_x))
    err_D = np.zeros(len(selected_x))
    
    rho_w = 1027

    for i,idx_x in enumerate(list_idx_x):
        FK = abs(mean_CWT[:,:,idx_x])
    
        k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
    
        H = get_water_height(DAS_water_height,UTC_mid,selected_x[i])
        popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                             bounds = (1e5,1e8))
        err_coeff = np.sqrt(np.diag(pcov))
        # print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
        D[i] = popt[0]
        err_D[i] = err_coeff[0]

    return D, err_D


#%% Load DAS parameter and Data

def main():
    date = '0211'
    
    # Load parameters for DAS
    main_path = '/media/turbots/DATA/Backup25/'
    path2DAS_param = f'{main_path}/Data/parameters_Febus_2025.pkl'
    with open(path2DAS_param,'rb') as pf:
        param = pickle.load(pf)
    print('Parameters file loaded')
    
    # Set parameters
    # fs = np.shape(strain_rate)[1] # time sampling frequency 
    # facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
    fs = param[date]['fs']
    fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
    facq_x = param[date]['facq_x'] 
    
    path2data = f'{main_path}/Data/{date}/DAS/'
    
    # Create folder for saving graphs
    fig_folder = f'{path2data}Figures/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
    
    filelist = glob.glob(path2data + '*.h5')
    Nb_minutes = 1 # duration of each stack
    
    idx_file = 0
    file2load = filelist[idx_file]
    stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    
    format_date = '%Y-%m-%dT%H-%M-%S'
    label_UTC0 = UTC_stack[0,0].strftime(format_date)
    current_fig_folder = f'{fig_folder}format_date/'
    if not os.path.isdir(current_fig_folder):
        os.mkdir(current_fig_folder)
        
    # Perform Fourier transform in time for all time chunks
    
    FFT_t, freq = all_chunks_FFT_t(stack_strain,fs)
     
    plot_time_spectrum_chunk(FFT_t,freq,s,UTC_stack,current_fig_folder,parula_map)
    
    
    # Define wavelet 
    
    wavelet = 'cmor1.5-1.5' # mother wavelet 
    [psi, x] = pywt.ContinuousWavelet(wavelet).wavefun(10)
    
    fig, axs = plt.subplots(figsize = (6,4))
    axs.plot(x,np.real(psi))
    axs.plot(x,np.imag(psi))
    axs.set_xlim([-5,5])
    
    figname = f'{current_fig_folder}wavelet_{wavelet}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    plt.close(fig)
    
    # Perform CWT for all time chunks 
    
    # parameters for CWT
    sampling_period = 1/facq_x #spatial sampling wavelength
    #logarithmic scales
    scales = np.geomspace(2,1024,num = 100)
    
    # compute normalized wavenumber (sigma)
    norm_sigma = pywt.scale2frequency(wavelet,scales)
    wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 
    print(f'Wavelengths sampled : min = {wavelength_sampled[0]:.2f}, max = {wavelength_sampled[-1]:.2f}')
    
    mean_CWT,k = CWT_all_chunks(FFT_t,scales,wavelet,sampling_period)
    
    # Load water height data 
    
    path2water_DAS = f'{main_path}/Data/{date}/DAS/'
    file2load = glob.glob(f'{path2water_DAS}fiber_water_height_GPS_structure_{date}.pkl')[0]
    with open(file2load,'rb') as pf:
        DAS_water_height = pickle.load(pf)
    
    
    # Plot (f,k) for different positions 
    
    selected_x = np.array([80,150,300,400,500,550])
    idx_x = [np.argmin(abs(s - x)) for x in selected_x]
    
    # plot FK without theory
    fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
    for ax, idx in zip(axs.flatten(),idx_x):
        imsh = ax.pcolormesh(k,freq,abs(mean_CWT[:,:,idx]),cmap = parula_map)
        ax.set_ylim([0,1])
        ax.set_xlim([0.01,0.6])
        # ax.set_ylim([0,10])
        # ax.set_xlim([0,1.0])
        ax.set_title(f's = {s[idx]:.2f} m')
    
    for j in range(axs.shape[1]):
        axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
    
    for i in range(axs.shape[0]):
        axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')
    
    figname = f'{current_fig_folder}{date}_FK_subplots_spectrum_spectrum_avg_{label_UTC0}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    
    # plot FK with theory 
    freq_range = [0.02,0.6]
    min_prominence = 0.9
    min_width = 10
    
    rho_w = 1027
    UTC_mid = UTC_stack[UTC_stack.shape[0]//2,0]
    
    fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
    for i,(ax,idx) in enumerate(zip(axs.flatten(),idx_x)):
        FK = mean_CWT[:,:,idx]
        imsh = ax.pcolormesh(k,freq,abs(FK),cmap = parula_map)
        ax.set_ylim([0,1])
        ax.set_xlim([0.01,0.6])
        # ax.set_ylim([0,8])
        # ax.set_xlim([0,0.7])
        ax.set_title(f's = {s[idx_x]:.2f} m')
        
        k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
        ax.plot(k_exp,f_exp,'r.')
    
        H = get_water_height(DAS_water_height,UTC_mid,selected_x[i])
        popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                             bounds = (1e5,1e8))
        err_coeff = np.sqrt(np.diag(pcov))
        print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
    
        k_th = np.linspace(0,1,100)
        omega_th = shallow_hydroelastic(k_th, popt[0], rho_w, H)
        label_th = r'$D = ' + f'{popt[0]:.2e}' + r'$'
        ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
        ax.legend(loc = 'lower right',fontsize = 10)
    
    for j in range(axs.shape[1]):
        axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
    
    for i in range(axs.shape[0]):
        axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')   
    
    figname = f'{fig_folder}{date}_subplots_D_extraction_spectrum_avg_{label_UTC0}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    
    # Compute flexural modulus D 
    
    zone1 = np.arange(40,210,step = 20)
    zone2 = np.arange(280,570,step = 20)
    selected_x = np.concatenate((zone1,zone2))
    
    peaks_detec = {}
    peaks_detec['freq_range'] = [0.02,0.6]
    peaks_detec['min_prominence'] = 0.9
    peaks_detec['min_width'] = 5
    
    UTC_mid = UTC_stack[UTC_stack.shape[0]//2,0]
    D,err_D = extract_D(mean_CWT, selected_x, s, freq, k, peaks_detec,DAS_water_height,UTC_mid)
    # save data 
    results_dict = {}
    results_dict['D'] = D
    results_dict['err_D'] = err_D
    results_dict['x'] = selected_x
    
    file2save = f'{path2data}{date}_wavelet_flexural_modulus_file_{label_UTC0}.pkl'
    with open(file2save, 'wb') as pf: 
        pickle.dump(results_dict,pf)
        
if __name__ == '__main__':
    main()
