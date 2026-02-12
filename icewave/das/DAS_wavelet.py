# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 16:08:03 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import h5py 
import glob
import os 
import pickle

from datetime import datetime
import pytz

import scipy
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pywt

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
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = cmap,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar


#%% load DAS parameters and data 

date = '0211'

# Load parameters for DAS
path2DAS_param = 'U:/Data/parameters_Febus_2025.pkl'
with open(path2DAS_param,'rb') as pf:
    param = pickle.load(pf)
print('Parameters file loaded')

# Set parameters
# fs = np.shape(strain_rate)[1] # time sampling frequency 
# facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
fs = param[date]['fs']
fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
facq_x = param[date]['facq_x'] 

path2data = f'U:/Data/{date}/DAS/'

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(path2data + '*.h5')
Nb_minutes = 1 # duration of each stack

idx_file = 0
file2load = filelist[idx_file]
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)

#%% Plot spatio-temporal 

chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, parula_map)
imsh.set_clim([-1e4,1e4])


#%% Compute FFT in time 

N = 2**(FT.nextpow2(stack_strain.shape[1]))

mean_sig = np.mean(stack_strain,axis = 1)
mean_sig = np.tile(mean_sig,(stack_strain.shape[1],1,1)).transpose((1,0,2))

FFT_t = np.fft.fft(stack_strain - mean_sig,n = N,axis = 1)
FFT_t = FFT_t/stack_strain.shape[1]

# keep only half of frequencies
FFT_t = FFT_t[:,:N//2,:]
FFT_t[:,1:-1,:] = 2*FFT_t[:,1:-1,:]

freq = fs*np.arange(0,(N/2))/N

fig, ax = plt.subplots()
imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
          extent = [freq[0],freq[-1],s[0],s[-1]])
imsh.set_clim([1e1,2e3])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$s \; \mathrm{(m)}$')
ax.set_xlim([0,10])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

#%% Plot profile along s for a given frequency

selected_freq = 0.2
idx_freq = np.argmin(abs(freq - selected_freq))

profile = FFT_t[0,idx_freq,:]

fig, ax = plt.subplots()
ax.plot(s,np.real(profile))
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$\hat{\dot{\epsilon}}(s)$')

#%% Perform continuous wavelet transform (cwt)

sampling_period = 1/facq_x #spatial sampling wavelength
wavelet = 'cmor1.5-1.0'

#logarithmic scales
scales = np.geomspace(2,512,num = 100)

# compute normalized wavenumber (sigma)
norm_sigma = pywt.scale2frequency(wavelet,scales)
wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 
# print(wavelength_sampled)

# compute cwt
cwtmatr, freqs = pywt.cwt(np.real(profile), scales, wavelet, sampling_period = sampling_period)

#%% Plot scaleogram

fig, ax = plt.subplots()
imsh = ax.pcolormesh(s,freqs,abs(cwtmatr[:-1,:-1]),cmap = parula_map)
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
ax.set_title(r'Continuous Wavelet Transform (scaleogram)')

#%% Check influence of wavelet used 

wavelet2 = 'cmor1.5-1.0'

# compute normalized wavenumber (sigma)
scales2 = scales
norm_sigma = pywt.scale2frequency(wavelet,scales2)
wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 

# compute cwt
cwtmatr, freqs = pywt.cwt(np.real(profile), scales2, wavelet2, sampling_period = sampling_period)

fig, ax = plt.subplots()
imsh = ax.pcolormesh(s,freqs,abs(cwtmatr[:-1,:-1]),cmap = parula_map)
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
ax.set_title(r'Continuous Wavelet Transform (scaleogram)')

#%% 

fig, ax = plt.subplots()
ax.plot(s,np.real(cwtmatr[44,:]))




