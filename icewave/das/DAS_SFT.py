# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 14:51:42 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.colors as colors 
import matplotlib.cm as cm

import h5py 
import glob
import os 
from time import strftime, localtime
from datetime import datetime
import pytz
import time 
import scipy
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pickle

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))


plt.rcParams.update({
    "text.usetex": True}) # use latex

global g
g = 9.81



#%% FUNCTION SECTION 

# function for plotting spatio temporal
def plot_spatio_temp(spatio,t,s,fiber_length,extents):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - t, numpy array or list, time array 
            - s, numpy array or list, curvilinear coordinate array
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsho, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#---------------------------------------------------------------------------------------------------

def wavenumbers_hydro(freq,rho_w,rho_ice,E,h,nu,H,c_w,equation):
    """ Compute wavenumbers associated to hydroelastic waves """ 
    
    g = 9.81
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    k = np.linspace(1e-4,10,200000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for i in range(len(freq)):
        omeg = 2*np.pi*freq[i]
        if omeg == 0:
            flag = 1
            idx_flag = i
        else:
            
            if equation == 'Squire_deep':
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + 1) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Squire_shallow':
                coth = 1/np.tanh(H*k)
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + coth) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Stein':
                cph = omeg/k
                func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            else : 
                print('inappropriate equation name, choose between : Squire_deep / Squire_shallow / Stein')

            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[i] = (np.where(np.diff(np.signbit(func)))[0]) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero] # wave vector associated to flexural mode 
    if flag:
        k_QS[idx_flag] = 0

    return k_QS  

#-----------------------------------------------------------------------------------------------------------------------------

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega

#-------------------------------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------------------------------
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

#%% load DAS parameters and data 

date = '0212'
# fig_folder = f'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2025/DAS/{date}/Figures/'


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

format_date = '%Y-%m-%d %H:%M:%S.f'
local_timearea = pytz.timezone('America/Montreal')
UTC_timearea = pytz.timezone('UTC')

# path2data = f'E:/Data/{date}/DAS_h5/'
path2data = f'F:/Rimouski_2025/Data/{date}/DAS/'
# path2data = f'U:/Data/{date}/DAS/'

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(path2data + '*UTC.h5')
Nb_minutes = 1 # duration of each stack


idx_file = -2
file2load = filelist[idx_file]
Nb_minutes = 1 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
down_sampling_factor = 10
date_downsampling = ['0210','0212']

# decimation for 0210 and 0212
if date in date_downsampling:
    fs = fs/down_sampling_factor # new value of sampling frequency
    stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
                                                                        stack_time,UTC_stack,down_sampling_factor)
    print(f'New value of sampling frequency, fs = {fs:.2f}')
        
#%% Plot spatio 

set_graphs.set_matplotlib_param('single')
chunk = 2
extents = [UTC_stack[chunk,0], UTC_stack[chunk,-1],s[0],s[-1]]
fig, ax , imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:],UTC_stack[chunk,:],s,fiber_length,extents = extents)

imsh.set(clim = [1, 5e3])

# figname = f'spatio_file_{idx_file}_chunk{chunk}_{Nb_minutes}min'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

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

# figname = f'FFT_t_file_{idx_file}_chunk{chunk}_{Nb_minutes}min'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Show a profile for a frequency 

selected_freq = 0.28
idx_freq = np.argmin(abs(freq - selected_freq))

profile = FFT_t[0,idx_freq,:]

fig, ax = plt.subplots()
ax.plot(s,np.real(profile))
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$\hat{\dot{\epsilon}}(s)$')

#%% Create a spatial Short-Time Fourier Transform 

# create hann window
max_wavelength = 200 # maximum wavelength we are looking for, in meter 
M = int(max_wavelength*facq_x) # number of points of window (in samples)
hann_win = scipy.signal.windows.hann(M,sym = True)

txt_wavelength = f'wavelengths measured : min = {2/facq_x:.2f} m, max = {M/facq_x:.2f} m'
print(txt_wavelength)

# create SFT transfrom
hop_meter = 20 # hops in meter between two STFT computations
hop = int(hop_meter*facq_x) # hops between two computations (in samples)
padding = 2**(FT.nextpow2(M))
SFT = scipy.signal.ShortTimeFFT(hann_win, hop,facq_x, fft_mode = 'onesided',mfft = padding,scale_to = 'magnitude')


#%% Plot spatial spectrogram for a given frequency 

Sx = SFT.spectrogram(np.real(profile),detr = 'constant',axis = -1)

# Plot spectrogram 
spectro_log = np.log10(abs(Sx))
fig, ax = plt.subplots()
imsh = ax.imshow(spectro_log,origin = 'lower',aspect = 'auto',extent = SFT.extent(stack_strain.shape[2]),
          cmap = parula_map)

imsh.set_clim([1,6])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(|\dot{\epsilon}|)$')

ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')

# figname = f'Spatial_spectro_file_{idx_file}_chunk{chunk}_{Nb_minutes}min_fdemod_{freq[idx_freq]}_M{M}_hop{hop}'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compute Spectrogram for all frequencies 

chunk = 2
Sfx = SFT.spectrogram(np.real(FFT_t[chunk,:,:]),detr = 'constant',axis = -1)

x_SFT = SFT.t(stack_strain.shape[2])

selected_x = 440 # distance from interogator
idx_x = np.argmin(abs(selected_x - SFT.t(stack_strain.shape[2])))
print(idx_x)

SFT_extents = SFT.extent(stack_strain.shape[2])

fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(abs(Sfx[:,:,idx_x])),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([0.5,1.5])

# ax.set_xlim([0,0.6])
# ax.set_ylim([0,1])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'FK for $s = ' +  f'{x_SFT[idx_x]:.2f}' + r'\; \mathrm{m}$')

# figname = f'FK_file_{idx_file}_s_{x[idx_x]}m_chunk{chunk}_{Nb_minutes}min_M{M}_hop{hop}_long_range'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot space-time spectrum for different positions 

selected_x = np.array([70,140,200,380,450,550])
idx_x = [np.argmin(abs(x_SFT - x)) for x in selected_x]

fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
for ax, idx in zip(axs.flatten(),idx_x):
    imsh = ax.imshow(np.log10(abs(Sfx[:,:,idx])),origin = 'lower',aspect = 'auto',cmap = parula_map, 
              extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
    imsh.set_clim([1,3])
    ax.set_ylim([0,1])
    ax.set_xlim([0.01,0.6])
    # ax.set_ylim([0,10])
    # ax.set_xlim([0,1.0])
    ax.set_title(f's = {x_SFT[idx]:.2f} m')

for j in range(axs.shape[1]):
    axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

for i in range(axs.shape[0]):
    axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')

# fig.colorbar(imsh,ax = axs[:,-1],shrink = 0.8, aspect = 10)

############################################################################
#%%---------- Average space-time spectrum over all chunks ------------------
############################################################################

Sfx = SFT.spectrogram(np.real(FFT_t[0,:,:]),detr = 'constant',axis = -1)
# compute spectrogram over all chunks
Sfx_all = np.zeros((FFT_t.shape[0],Sfx.shape[0],Sfx.shape[1],Sfx.shape[2]))
for i in range(Sfx_all.shape[0]):
    Sfx_all[i,:,:,:] = SFT.spectrogram(np.real(FFT_t[i,:,:]),detr = 'constant',axis = -1)

x_SFT = SFT.t(stack_strain.shape[2])
SFT_extents = SFT.extent(stack_strain.shape[2])

# average space-time spectrum 
mean_Sfx = np.mean(abs(Sfx_all),axis = 0)


#%% Plot for a given position 

selected_x = 150 # distance from interogator
idx_x = np.argmin(abs(selected_x - SFT.t(stack_strain.shape[2])))
print(idx_x)

fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx_x]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([1,2])

# ax.set_xlim([0,0.6])
# ax.set_ylim([0,1])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'FK for $s = ' +  f'{x_SFT[idx_x]:.2f}' + r'\; \mathrm{m}$')

#%% Plot space-time spectrum average over all chunks for different positions

selected_x = np.array([100,200,350,450,500,550])
idx_x = [np.argmin(abs(x_SFT - x)) for x in selected_x]

fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
for ax, idx in zip(axs.flatten(),idx_x):
    imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
              extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
    imsh.set_clim([1,2])
    ax.set_ylim([0,1])
    ax.set_xlim([0.01,0.6])
    # ax.set_ylim([0,10])
    # ax.set_xlim([0,1.0])
    ax.set_title(f's = {x_SFT[idx]:.2f} m')

for j in range(axs.shape[1]):
    axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

for i in range(axs.shape[0]):
    axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')

#%% Plot profile for a given frequency 

selected_x = 440
idx_x = np.argmin(abs(x_SFT - selected_x))

freq_min = 0.1
freq_max = 0.6
idx_min = np.argmin(abs(freq - freq_min))
idx_max = np.argmin(abs(freq - freq_max))

indices = np.arange(idx_min,idx_max + 1, step = 1)

# define norm for colormaps
norm = colors.Normalize(vmin = freq[idx_min], vmax = freq[idx_max])

min_prominence = 0.7
min_width = 0
# min_height = 0.9

dict_peaks = {}
fig, ax = plt.subplots()

for idx_freq in indices:
# print(f'frequency = {freq[idx_freq]:.2f} Hz')
    current_max = np.max(abs(mean_Sfx[idx_freq,:,idx_x]))
    # find peaks
    peaks,properties = scipy.signal.find_peaks(abs(mean_Sfx[idx_freq,:,idx_x])/current_max,
                                               height = min_prominence)
    if len(peaks) != 0:
        dict_peaks[str(idx_freq)] = peaks
    current_color = new_blues(norm(freq[idx_freq]))
    ax.plot(2*np.pi*SFT.f,abs(mean_Sfx[idx_freq,:,idx_x])/current_max,'-o',color = current_color)
    ax.plot(2*np.pi*SFT.f[peaks],abs(mean_Sfx[idx_freq,peaks,idx_x])/current_max,'r.')

ax.set_xlim([0.01,0.6])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r'$f \; \mathrm{(Hz)} $')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

#%% Load DAS water height 
path2water_DAS = f'U:/Data/{date}/DAS/'
file2water = glob.glob(f'{path2water_DAS}fiber_water_height_GPS_structure_{date}.h5')[0]
print(file2water)
DAS_water_height = rw.load_dict_from_h5(file2water)
print(DAS_water_height.keys())
# convert UTC string to datetime 
format_date = '%Y-%m-%dT%H-%M-%S.%f'
# convert UTC string to datetime 
DAS_water_height['UTC_datetime'] = []
for date_txt in DAS_water_height['UTC_t']:
    if date_txt != 'None' :
        datetime_obj = datetime.strptime(date_txt,format_date)
        datetime_obj = datetime_obj.replace(tzinfo = pytz.timezone('UTC'))
    # else :
    #     datetime_obj = None
    DAS_water_height['UTC_datetime'].append(datetime_obj)
    
DAS_water_height['UTC_t'] = np.array(DAS_water_height['UTC_datetime'])
del DAS_water_height['UTC_datetime']

UTC_chunk = UTC_stack[0,0] #initial UTC time of the studied chunk 

#%% Extract peaks from space-time spectrum and compute D 

selected_x = np.array([100,200,350,450,500,550])
idx_x = [np.argmin(abs(x_SFT - x)) for x in selected_x]

D = np.zeros(len(selected_x))
err_D = np.zeros(len(selected_x))

freq_range = [0.1,0.6]
min_prominence = 0.6
min_width = 1

rho_w = 1027

fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
for i,(ax,idx) in enumerate(zip(axs.flatten(),idx_x)):
    FK = mean_Sfx[:,:,idx]
    imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
              extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
    imsh.set_clim([1,6])
    ax.set_ylim([0,1])
    ax.set_xlim([0.01,0.6])
    # ax.set_ylim([0,8])
    # ax.set_xlim([0,0.7])
    ax.set_title(f's = {x_SFT[idx]:.2f} m')
    
    k_exp,f_exp = extract_peaks(FK,freq,2*np.pi*SFT.f,freq_range,min_prominence,min_width)
    ax.plot(k_exp,f_exp,'r.')

    H = get_water_height(DAS_water_height,UTC_chunk,selected_x[i])
    popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                         bounds = (1e5,1e8))
    err_coeff = np.sqrt(np.diag(pcov))
    print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
    D[i] = popt[0]
    err_D[i] = err_coeff[0]

    k_th = np.linspace(0,1,100)
    omega_th = shallow_hydroelastic(k_th, D[i], rho_w, H)
    label_th = r'$D = ' + f'{popt[0]:.2e}' + r'$'
    ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
    ax.legend(loc = 'lower right',fontsize = 10)

for j in range(axs.shape[1]):
    axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

for i in range(axs.shape[0]):
    axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')   

# fig_folder = 'F:/Rimouski_2025/Data/0211/DAS/Figures/'
# figname = f'{fig_folder}{date}_STFT_subplots_D_extraction_space_time_spectrum_avg_file0'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compute D for all positions of STFT 

selected_x = x_SFT[np.logical_and(x_SFT >= 0, x_SFT < 600)]
idx_x = [np.argmin(abs(x_SFT - x)) for x in selected_x]

D = np.zeros(len(selected_x))
err_D = np.zeros(len(selected_x))

freq_range = [0.02,0.6]
min_prominence = 0.9
min_width = 2

rho_w = 1027

for i,idx in enumerate(idx_x):
    FK = mean_Sfx[:,:,idx]

    k_exp,f_exp = extract_peaks(FK,freq,2*np.pi*SFT.f,freq_range,min_prominence,min_width)
    ax.plot(k_exp,f_exp,'r.')

    H = get_water_height(DAS_water_height,UTC_chunk,selected_x[i])
    popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                         bounds = (1e5,1e8))
    err_coeff = np.sqrt(np.diag(pcov))
    print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
    D[i] = popt[0]
    err_D[i] = err_coeff[0]

# save data 
results_dict = {}
results_dict['D'] = D
results_dict['err_D'] = err_D
results_dict['x'] = selected_x

file2save = f'F:/Rimouski_2025/Data/{date}/DAS/{date}_STFT_flexural_modulus_file_{idx_file}.pkl'
with open(file2save, 'wb') as pf: 
    pickle.dump(results_dict,pf)

#%% Plot data 

fig, ax = plt.subplots()
ax.plot(selected_x,D,'o')
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')

figname = f'{fig_folder}{date}_STFT_D_VS_s_spectrum_avg_file0'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')


#%% Compare data from wavelet analysis and STFT analysis

file_STFT = f'F:/Rimouski_2025/Data/{date}/DAS/{date}_STFT_flexural_modulus_file_0.pkl'
with open(file_STFT,'rb') as pf: 
    results_STFT = pickle.load(pf)

file_wavelet = f'F:/Rimouski_2025/Data/{date}/DAS/{date}_wavelet_flexural_modulus_file_0.pkl'
with open(file_wavelet,'rb') as pf: 
    results_wavelet = pickle.load(pf)
    
#%%
fig, ax = plt.subplots()
ax.errorbar(results_STFT['x'],results_STFT['D'],yerr = results_STFT['err_D'],fmt = 'o',label = 'STFT')
ax.errorbar(results_wavelet['x'],results_wavelet['D'],yerr = results_wavelet['err_D'],fmt = 'o',label = 'Wavelet')

ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')

ax.legend()

figname = f'{fig_folder}{date}_comparison_STFT_wavelet_D_extraction_spectrum_avg_file0'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')











#%% Average time FFT over all chunks 

mean_FFT = np.mean(FFT_t,axis = 0)
mean_Sfx = SFT.spectrogram(np.real(mean_FFT),detr = 'constant',axis = -1)

#%% Plot FK averaged over time for a specific abscisse

selected_x = 80 # distance from interogator
idx_x = np.argmin(abs(selected_x - SFT.t(stack_strain.shape[2])))
print(idx_x)

SFT_extents = SFT.extent(stack_strain.shape[2])

fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx_x]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [2*np.pi*SFT_extents[2],2*np.pi*SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([1,5])

# ax.set_xlim([0,0.1])
# ax.set_ylim([0,1.0])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'Average FK for $s = ' +  f'{x[idx_x]:.2f}' + r'\; \mathrm{m}$')

#%%
# Select points on the graph
print('Select points on FK plots')
Nb_points = 10
points = plt.ginput(Nb_points)

sigma,f_exp = zip(*points)
k_exp = 2*np.pi*np.array(sigma)
ax.plot(sigma,f_exp,'ro')

# fit experimental points by shallow hydroelastic 
rho_w = 1027
H = 4
popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                     bounds = (1e5,1e8))

err_coeff = np.sqrt(np.diag(pcov))
print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')

k_th = np.linspace(0,1,100)
omega_th = shallow_hydroelastic(k_th, popt, rho_w, H)
label_th = r'$D = ' + f'{popt[0]:.2e}' + r'$ , $H = ' + f'{H:.1f}' + r'$'
ax.plot(k_th/2/np.pi,omega_th/2/np.pi,'w--',label = label_th)

ax.legend()
# figname = f'FK_file_{idx_file}_s_{x[idx_x]}m_{Nb_minutes}min_M{M}_hop{hop}_short_range'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')











        