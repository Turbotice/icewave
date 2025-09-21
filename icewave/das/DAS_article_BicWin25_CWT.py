# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 09:24:24 2025

@author: sebas

Create figures for DAS article 
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 

import h5py

import scipy 
import pywt


import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs
import icewave.das.DAS_package as DS
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather

# PARULA COLORMAP 
parula_map = matcmaps.parula()

global g
g = 9.81
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10

#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/CWT/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
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
    ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
    
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
        imsh.set_clim([1e1,3e2])
        
        ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
        ax.set_ylabel(r'$s \; \mathrm{(m)}$')
        ax.set_xlim([0,10])
        
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
        
        if len(k_exp) != 0:
            
            H = get_water_height(DAS_water_height,UTC_mid,selected_x[i])
            popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                                 bounds = (1e5,1e8))
            err_coeff = np.sqrt(np.diag(pcov))
            # print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
            D[i] = popt[0]
            err_D[i] = err_coeff[0]
            
        else:
            D[i] = None
            err_D[i] = None

    return D, err_D

#%% Load DAS parameters and data

# Load parameters for DAS
main_path = 'U:/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

date = '0211'
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 5
file2load = filelist[idx_file]
print(file2load)

Nb_minutes = 1 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
format_date = '%Y-%m-%dT%H-%M-%S'
label_UTC0 = UTC_stack[0,0].strftime(format_date)

# decimation for 0210 and 0212
if date in date_downsampling:
    fs = fs/down_sampling_factor # new value of sampling frequency
    stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
                                                                        stack_time,UTC_stack,down_sampling_factor)
    print(f'New value of sampling frequency, fs = {fs:.2f}')


#%% Show spatio-temporal 

chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
imsh.set_clim([-1e4,1e4])
ax.set_xlabel(r'UTC time (s)')

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

# figname = f'{fig_folder}spatio_temporal_{label_UTC0}_chunk_{chunk}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Perform Fourier transform in time for all time chunks

FFT_t, freq = all_chunks_FFT_t(stack_strain,fs)

fig, ax = plt.subplots()
imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
          extent = [freq[0],freq[-1],s[0],s[-1]])
imsh.set_clim([0,2e3])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$x \; \mathrm{(m)}$')
ax.set_xlim([0,10])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)


#%% Show Demodulated signal 

selected_freq = 0.20
idx_freq = np.argmin(abs(freq - selected_freq))
print(idx_freq)
profile = FFT_t[0,idx_freq,:]

fig, ax = plt.subplots()
ax.plot(s,np.real(profile))
ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$\hat{\dot{\epsilon}} \; \mathrm{(u.a.)}$')

ax.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))

freq_txt = f'f_{selected_freq:.2f}'.replace('.','p')

figname = f'{fig_folder}Filtered_signal_{freq_txt}_{label_UTC0}_chunk_{chunk}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Define wavelet

wavelet = 'cmor1.5-1.5' # mother wavelet 
[psi, x] = pywt.ContinuousWavelet(wavelet).wavefun(10)

fig, ax = plt.subplots()
ax.plot(x,np.real(psi), label = r'$\mathrm{Re}(\psi)$')
ax.plot(x,np.imag(psi), label = r'$\mathrm{Im}(\psi)$')
ax.set_xlim([-5,5])

ax.set_ylabel(r'$\psi(u)$')
ax.set_xlabel(r'$u$')

ax.legend()

figname = f'{fig_folder}Wavelet_morlet_{label_UTC0}_chunk_{chunk}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
#%% Compute scaleogram for this specific frequency 

sampling_period = 1/facq_x #spatial sampling wavelength
#logarithmic scales
scales = np.geomspace(2,1024,num = 100)

# compute normalized wavenumber (sigma)
norm_sigma = pywt.scale2frequency(wavelet,scales)
wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 
print(f'Wavelengths sampled : min = {wavelength_sampled[0]:.2f} m, max = {wavelength_sampled[-1]:.2f} m')

# compute cwt
cwtmatr, freqs = pywt.cwt(np.real(profile), scales, wavelet, sampling_period = sampling_period)

#%% Plot scaleogram

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
k = 2*np.pi*freqs
imsh = ax.pcolormesh(s,k,abs(cwtmatr[:,:]),cmap = 'gist_yarg',shading = 'auto',norm = 'linear')
ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$')
imsh.set_clim([0,6000])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_ylim([0.07,0.5])
# ax.set_yscale('log')

cbar.set_label(r'$\hat{\dot{\epsilon}} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

figname = f'{fig_folder}Scaleogram_{freq_txt}_{label_UTC0}_chunk_{chunk}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compile 4 figures into one 
set_graphs.set_matplotlib_param('square')
fig, ax = plt.subplots(2,2,figsize = (12,9), layout = 'constrained')
# other technique = use gridspec
# fig.subplots_adjust(wspace = 0.1)

# Plot spatio 
quad_spatio = (0,0)
imsh = ax[quad_spatio].imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = 'linear', cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax[quad_spatio].set_ylim([0,fiber_length])

# divider = make_axes_locatable(ax[0,0])
# cax = divider.append_axes("right", size="2%", pad=0.1)
# cbar = plt.colorbar(imsh,cax = cax)
# cbar = plt.colorbar(imsh,ax = ax[quad_spatio],shrink = 0.51)

ax[quad_spatio].set_xlabel(r'UTC time')
ax[quad_spatio].set_ylabel(r'$x \; \mathrm{(m)}$')

imsh.set_clim([-1e4,1e4])

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

# Plot filtered signal
quad_filtered = (0,1)
selected_freq = 0.20
idx_freq = np.argmin(abs(freq - selected_freq))
print(idx_freq)
profile = FFT_t[0,idx_freq,:]

ax[quad_filtered].plot(s,np.real(profile))
ax[quad_filtered].set_xlabel(r'$x \; \mathrm{(m)}$')
ax[quad_filtered].set_ylabel(r'$\hat{\dot{\epsilon}} \; \mathrm{(u.a.)}$')

ax[quad_filtered].ticklabel_format(axis='y', style='sci', scilimits=(3, 3))

# Plot wavelet
quad_wavelet = (1,0)
wavelet = 'cmor1.5-1.5' # mother wavelet 
[psi, x] = pywt.ContinuousWavelet(wavelet).wavefun(10)

ax[quad_wavelet].plot(x,np.real(psi), label = r'$\mathrm{Re}(\psi)$')
ax[quad_wavelet].plot(x,np.imag(psi), label = r'$\mathrm{Im}(\psi)$')
ax[quad_wavelet].set_xlim([-5,5])

ax[quad_wavelet].set_ylabel(r'$\psi(u)$')
ax[quad_wavelet].set_xlabel(r'$u$')

ax[quad_wavelet].legend()

# Plot scaleogram
quad_scaleo = (1,1)
imsh = ax[quad_scaleo].pcolormesh(s,k,abs(cwtmatr[:,:]),cmap = 'gist_yarg',shading = 'auto',norm = 'linear')
ax[quad_scaleo].set_xlabel(r'$x \; \mathrm{(m)}$')
ax[quad_scaleo].set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$')
imsh.set_clim([0,6000])

ax[quad_scaleo].set_ylim([k.min(),1.0])
ax[quad_scaleo].set_yscale('log')

# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="2%", pad=0.1)
# cbar = plt.colorbar(imsh,cax = cax)

cbar = plt.colorbar(imsh,ax = ax[quad_scaleo], shrink = 1)

cbar.set_label(r'$\hat{\dot{\epsilon}} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

# figname = f'{fig_folder}Quadrant_spatio_wavelet_scaleo_f_{freq_txt}_{label_UTC0}_chunk_{chunk}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute Continuous Wavelet Transform for all frequencies and all chunks










