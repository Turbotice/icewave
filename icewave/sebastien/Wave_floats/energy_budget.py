# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 09:16:15 2025

@author: sebas
"""

import numpy as np
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
import pickle
import scipy.signal
from scipy.signal import find_peaks


# import modules 
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs

parula_map = matcmaps.parula()

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/Referee_answer/'

#%% Load data for object 

base = 'F:/Waves_reconstruction_wilson/9_wave_field_2024_01_02/'
exp_ID = 'D4cm_h10mm_fps_100'
test_ID = 'f6.0Hz_amplitude15mm'

path2data = f'{base}{exp_ID}/{test_ID}/'
filelist = glob.glob(f'{path2data}*Height_in_cm_surf.mat')
file2load = filelist[0]

data = {'obj':{},'wave':{}}
with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    for i in range(1,len(list_keys) - 1):
        key = list_keys[i]
        data['obj'][key] = mat2py.mat_to_dict(fmat[key],fmat[key])
        
    data['obj']['table'] = mat2py.mat_to_dict(fmat['table'],fmat['#refs#'])
    
    data['obj']['H_ww'] = np.flip(data['obj']['H_ww'],axis = 2)

#%% Load data for associated wave field 

# base_wavefield = 'W:/Banquise/Sebastien/Test_Disks_Wilson_Seb/9_wavefield_2024_01_03/wave_field/'
base_wavefield = 'F:/Waves_reconstruction_wilson/9_wave_field_2024_01_02/wave_field/'
test_ID = 'f6p0Hz_amplitude_15mm'
path2wavefield = f'{base_wavefield}{test_ID}'
file2load = f'{path2wavefield}/Height_in_cm_surf.mat'

with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    for i in range(1,len(list_keys) - 1):
        key = list_keys[i]
        data['wave'][key] = mat2py.mat_to_dict(fmat[key],fmat[key])
        
    data['wave']['table'] = mat2py.mat_to_dict(fmat['table'],fmat['#refs#'])
    
    data['wave']['H_ww'] = np.flip(data['wave']['H_ww'],axis = 2)
    
    

#%% Show a frame with object

frame = 250
field = data['obj']['H_ww'][frame,:,:]

extents = np.array([data['obj']['X'][0] , data['obj']['X'][-1], 
                    data['obj']['Y'][0], data['obj']['Y'][-1]])/data['obj']['fx']*1e2 # in centimeters 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(field.T,cmap = parula_map, aspect = 'equal', 
                 origin = 'lower', vmin = -0.2, vmax = 0.2) 

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$')


ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)


#%% Show same frame without object 

frame = 250
field = data['wave']['H_ww'][frame,:,:]

extents = np.array([data['wave']['X'][0] , data['wave']['X'][-1], 
                    data['wave']['Y'][0], data['wave']['Y'][-1]])/data['wave']['fx']*1e2 # in centimeters 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(field.T,cmap = parula_map, aspect = 'equal', 
                 origin = 'lower', extent = extents, vmin = -0.2, vmax = 0.2) 

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$')


ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)

###########################################################################################
#%% -------------------------- X-PROFILE -------------------------------------------------
###########################################################################################

#Show profile along x of the pure wavefield 

idx_mid = data['wave']['H_ww'].shape[2]//2
# spatiotemporal
spatio = data['wave']['H_ww'][:,:,idx_mid]

facq_x = data['wave']['fx']
fs = data['wave']['table'][1,4]
x = np.arange(0,spatio.shape[1])/facq_x
t = np.arange(0,spatio.shape[0])/fs
extents = np.array([t.min(),t.max(),x.min(),x.max()])

fig, ax = plt.subplots()
ax.imshow(spatio.T, origin = 'lower',cmap = parula_map, extent = extents, aspect = 'auto')
ax.set_xlabel(r'$t \; \mathrm{(s)}$')
ax.set_ylabel(r'$x \; \mathrm{(m)}$')


#%% perform FFT
# FFT 2D of the spatio
facq= [fs,facq_x]
FFT, omega, k = FT.fft_2D(spatio,facq,add_pow2 = [0,0])
k = -k

#%% Show Fourier spectrum 

# set low frequencies to zeros 
# fc = 0.5 
# idx_max = np.argmin(abs(omega - 2*np.pi*fc))
# idx_min = np.argmin(abs(omega + 2*np.pi*fc))
# FFT[idx_min:idx_max,:] = 0

# set low wavenumbers to zeros
kc = 30
idx_min = np.argmin(abs(k - kc))
idx_max = np.argmin(abs(k + kc))
FFT[:,idx_min:idx_max] = 0

# set wave propagating along -kx to zeros

FFT[FFT.shape[0]//2:,FFT.shape[1]//2:] = 0
FFT[:FFT.shape[0]//2,:FFT.shape[1]//2] = 0

extents_fourier = np.array([k[0],k[-1],omega[0],omega[-1]])
fig, ax = plt.subplots()
ax.imshow(abs(FFT),origin  = 'lower',cmap = parula_map, extent = extents_fourier, aspect = 'auto')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')
ax.set_xlim([-400,400])
ax.set_ylim([-60,60])


#%% Inverse Fourier transform 

inv_FFT2 = FT.inverse_FFT2(FFT,spatio.shape)

frame = 0
profile = spatio[frame,:]
profile2 = np.real(inv_FFT2[frame,:])
fig, ax = plt.subplots()
ax.plot(x,profile)
ax.plot(x,profile2)

#%% Demodulate 

f_ex = 6.0 # excitation frequency 

modul = np.exp(1j*2*np.pi*f_ex*t)
modul = np.tile(modul,(spatio.shape[1],1)).T

demod = np.mean(spatio * modul,axis = 0)
demod2 = np.mean(np.real(inv_FFT2) * modul, axis = 0)

fig, ax = plt.subplots()
# ax.plot(x,np.real(demod),'-')
ax.plot(x,np.real(demod2),'-')
ax.plot(x,abs(demod2),'r-')
ax.set_ylim([-0.06,0.06])

#############################################################################################
#%% ------------------------------ Y-PROFILE WITH OBJECT --------------------------------- 
#############################################################################################

# load cropped matrix
data['cropped'] = {}
file2load = f'{path2data}croppedMatrixobject.mat'

with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    
    data['cropped'] = mat2py.mat_to_dict(fmat['croppedMatrix'],fmat['croppedMatrix'])

#%% Plot a frame 
idx_center = [218,248] # indices for object center

frame = 250
field = data['cropped'][frame,:,:]
x_crop = np.arange(0,data['cropped'].shape[1])/facq_x
y_crop = np.arange(0,data['cropped'].shape[2])/facq_x
extents_crop = np.array([x_crop[0],x_crop[-1],y_crop[0],y_crop[-1]])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(field.T,cmap = parula_map, aspect = 'equal', 
                 origin = 'lower', extent = extents_crop, vmin = -0.2, vmax = 0.2) 

ax.plot(x_crop[idx_center[0]],y_crop[idx_center[1]],'ro')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$')


ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)

#%% Select profile going trhough object center

spatio = data['cropped'][:,idx_center[0],:]

extents = [t[0],t[-1],y_crop[0],y_crop[-1]]
fig, ax = plt.subplots()
ax.imshow(spatio.T, origin = 'lower',cmap = parula_map, extent = extents, aspect = 'auto')
ax.set_xlabel(r'$t \; \mathrm{(s)}$')
ax.set_ylabel(r'$x \; \mathrm{(m)}$')

#%% perform FFT
# FFT 2D of the spatio
facq= [fs,facq_x]
FFT, omega, k = FT.fft_2D(spatio,facq,add_pow2 = [0,0])
k = -k

#%% Show Fourier spectrum 

# set low frequencies to zeros 
# fc = 0.5 
# idx_max = np.argmin(abs(omega - 2*np.pi*fc))
# idx_min = np.argmin(abs(omega + 2*np.pi*fc))
# FFT[idx_min:idx_max,:] = 0

# set low wavenumbers to zeros
# kc = 30
# idx_min = np.argmin(abs(k - kc))
# idx_max = np.argmin(abs(k + kc))
# FFT[:,idx_min:idx_max] = 0

# set wave propagating along -kx to zeros

# FFT[FFT.shape[0]//2:,FFT.shape[1]//2:] = 0
# FFT[:FFT.shape[0]//2,:FFT.shape[1]//2] = 0

extents_fourier = np.array([k[0],k[-1],omega[0],omega[-1]])
fig, ax = plt.subplots()
ax.imshow(abs(FFT),origin  = 'lower',cmap = parula_map, extent = extents_fourier, aspect = 'auto')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')
ax.set_xlim([-400,400])
ax.set_ylim([-60,60])


#%% ------------------- DRAFT ---------------------------
#%% demodulate profile at excitation frequency 

f_ex = 6.0 # excitation frequency 
product = np.exp(1j*2*np.pi*f_ex*t) @ spatio 

#%%
modul = np.exp(1j*2*np.pi*f_ex*t)
modul = np.tile(modul,(spatio.shape[1],1)).T

demod = np.mean(spatio * modul,axis = 0)

fig, ax = plt.subplots()
ax.plot(x,np.real(demod),'-')

#%%
# unshift fft
unshift = np.fft.ifftshift(FFT)
inverse = np.fft.ifft2(unshift)*np.size(spatio)

frame = 0
profile = spatio[frame,:]
profile2 = inverse[frame,:spatio.shape[1]]
fig, ax = plt.subplots()
ax.plot(x,profile)
ax.plot(x,np.real(profile2))

#%% load object positions 

import scipy


data['pos'] = {}
file2load = f'{path2data}locss.mat'

test = scipy.io.loadmat(file2load)



