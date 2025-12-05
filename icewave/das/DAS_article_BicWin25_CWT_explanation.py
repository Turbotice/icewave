# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 14:17:09 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib as mpl
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

full_yarg = mpl.colormaps['gist_yarg'].resampled(256)
new_yarg = colors.ListedColormap(full_yarg(np.linspace(0.1,1,256)))

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))


global g
g = 9.81
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')


#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/CWT_explanation/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% FUNCTION SECTION 

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

def subpixel_polyfit(x,y,p,deg = 2):
    """ Compute subpixel position of a local maximum / minimum 
    Inputs : - x, array of abscisse 
             - y, array of y-axis values
             - p, index at which a peak has been detected 
    Outputs : - local_argmax, subpixel value for which the maximum of degree 2 polynomial fit is reached 
              - local_max, maximum of the degree 2 polynomial """
    
    x2fit = np.array([x[p-1],x[p],x[p+1]])
    y2fit = np.array([y[p-1],y[p],y[p+1]])
    coeffs = np.polyfit(x2fit,y2fit,deg = 2)
    local_argmax = -coeffs[1]/2/coeffs[0]
    local_max = np.polyval(coeffs,local_argmax)

    return local_argmax,local_max

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
        normalized = abs(FK[idx_freq,:])/current_max
        # find peaks
        peaks,properties = scipy.signal.find_peaks(normalized,
                                                   prominence = min_prominence,width = min_width)
        if len(peaks) != 0:
            for peak in peaks:
                local_argmax,local_max = subpixel_polyfit(k,normalized,peak)
                k_exp.append(local_argmax)
                f_exp.append(freq[idx_freq])
    
    k_exp = np.array(k_exp)
    f_exp = np.array(f_exp)

    return k_exp,f_exp

#-------------------------------------------------------------------------------------------------------------------

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

def load_DAS_water_height(file2water):
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
    
    return DAS_water_height

    
#%% Load parameters for DAS
main_path = 'U:/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

date = '0211'
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 5 #5 for 0211
file2load = filelist[idx_file]
print(file2load)

Nb_minutes = 1 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
format_date = '%Y-%m-%dT%H-%M-%S'
label_UTC0 = UTC_stack[0,0].strftime(format_date)

# shift curvilinear axis
offset_fiber = 37.5
s = s - offset_fiber

# decimation for 0210 and 0212
# if date in date_downsampling:
#     fs = fs/down_sampling_factor # new value of sampling frequency
#     stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
#                                                                         stack_time,UTC_stack,down_sampling_factor)
#     print(f'New value of sampling frequency, fs = {fs:.2f}')


#%% Perform Fourier transform in time for all time chunks

chunk = 0
FFT_t, freq = all_chunks_FFT_t(stack_strain,fs)

fig, ax = plt.subplots()
imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
          extent = [freq[0],freq[-1],s[0],s[-1]])
imsh.set_clim([0,2e3]) # 2e3 for 0211

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$x \; \mathrm{(m)}$')
ax.set_xlim([0,10])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

#%% Show Demodulated signal 

freq_list = [0.2,0.22,0.25,0.3,0.4,0.5]

for selected_freq in freq_list:

    idx_freq = np.argmin(abs(freq - selected_freq))
    print(idx_freq)
    profile = FFT_t[0,idx_freq,:]
    max_profile = max(np.real(profile))
    
    set_graphs.set_matplotlib_param('triple')
    fig, ax = plt.subplots()
    ax.plot(s,np.real(profile))
    ax.set_xlabel(r'$x \; \mathrm{(m)}$')
    ax.set_ylabel(r'$\hat{\dot{\epsilon}} \; \mathrm{(a.u.)}$')
    max_y = max_profile*1.1
    ax.set_ylim([-max_y,max_y])
    ax.set_xlim([0,s[-1]])
    # ax.set_ylim([-1.1,1.1])
    label = r'$f = ' + f'{selected_freq:.2f}' + '\; \mathrm{Hz}$'
    ax.set_title(label)
    
    ax.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
    
    freq_txt = f'f_{selected_freq:.2f}'.replace('.','p')
    
    figname = f'{fig_folder}Filtered_signal_{freq_txt}_{label_UTC0}_chunk_{chunk}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
#%% Compute scaleogram for this specific frequency

wavelet = 'cmor1.5-1.5' # mother wavelet  

sampling_period = 1/facq_x #spatial sampling wavelength
#logarithmic scales
scales = np.geomspace(2,1024,num = 100)

# compute normalized wavenumber (sigma)
norm_sigma = pywt.scale2frequency(wavelet,scales)
wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 
print(f'Wavelengths sampled : min = {wavelength_sampled[0]:.2f} m, max = {wavelength_sampled[-1]:.2f} m')

#%% Plot scaleogram

for selected_freq in freq_list:
    
    idx_freq = np.argmin(abs(freq - selected_freq))
    print(idx_freq)
    profile = FFT_t[0,idx_freq,:]

    # compute cwt
    cwtmatr, freqs = pywt.cwt(np.real(profile), scales, wavelet, sampling_period = sampling_period)
    
    set_graphs.set_matplotlib_param('triple')
    fig, ax = plt.subplots()
    k = 2*np.pi*freqs
    
    max_scaleo = np.max(abs(cwtmatr))
    imsh = ax.pcolormesh(s,k,abs(cwtmatr[:,:])/max_scaleo,cmap = 'gist_yarg',shading = 'auto',norm = 'linear')
    ax.set_xlabel(r'$x \; \mathrm{(m)}$')
    ax.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$')
    imsh.set_clim([0,1])
    imsh.set_rasterized(True)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)
    cbar.set_label(r'$\hat{\dot{\epsilon}}(f)/\hat{\dot{\epsilon}}_{max}(f)$')
    # cbar.formatter.set_powerlimits((3, 3))
    # cbar.update_ticks()
    
    ax.set_ylim([k.min(),0.5])
    # ax.set_yscale('log')
    ax.set_xlim([0,s[-1]])
    
    label = r'$f = ' + f'{selected_freq:.2f}' + '\; \mathrm{Hz}$'
    ax.set_title(label)

    freq_txt = f'f_{selected_freq:.2f}'.replace('.','p')
    
    figname = f'{fig_folder}Scaleogram_{freq_txt}_{label_UTC0}_chunk_{chunk}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
#%% Plot Dispersion relation for a given x
# =============================================================================

#%% Load avg CWT

date = '0211'
# main_path = f'U:'

main_path = 'U:/'
path2data = f'{main_path}Data/{date}/DAS/'

path2CWT = f'{path2data}avg_CWT/'
filelist = glob.glob(f'{path2CWT}avg_CWT*.h5')
print(filelist)

# load water height data
file2water = glob.glob(f'{path2data}fiber_water_height_GPS_structure_{date}.h5')[0]
print(file2water)    
DAS_water_height = load_DAS_water_height(file2water)

# load swell orientation 
file_swell_orientation = f'{main_path}Data/swell_orientation_beam_forming.h5'
swell_orientation = rw.load_dict_from_h5(file_swell_orientation)

swell_bool = 1 #1 if we want to account for swell orientation correction 
swell_dict =  {'bool':swell_bool,'theta':swell_orientation[date]}

idx_file = 5
file2load = filelist[idx_file]
print(file2load)
data = rw.load_dict_from_h5(file2load)

#%% define label 
label_UTC0 = data['label']

# change k if we want to account for swell orientation 
if swell_dict['bool']:
    data['k_star'] = data['k_star']/np.cos(swell_dict['theta']*np.pi/180)
    label_swell = 'swell_corrected'
else:
    label_swell = ''

#%% Plot (f,k) and D fit for three different positions 

s = data['x'] - offset_fiber
k = data['k_star']
freq = data['f']
mean_CWT = data['mean_CWT']

selected_x = np.array([100,300,500]) - 35
idx_x = [np.argmin(abs(s - x)) for x in selected_x]

freq_range = [0.02,0.8]
min_prominence = 0.7
min_width = 3

# select range of k for which we look for a peak
if date in date_downsampling:
    k_min = 0.05
    k_max = 0.8
    idx_kmin = np.argmin(abs(k - k_min))
    idx_kmax = np.argmin(abs(k - k_max))
    mean_CWT = mean_CWT[:,idx_kmax:idx_kmin,:]
    k = k[idx_kmax:idx_kmin]

rho_w = 1027

format_date = '%Y-%m-%dT%H-%M-%S'
UTC_0 = datetime.strptime(data['label'],format_date)
UTC_0 = UTC_0.replace(tzinfo = pytz.timezone('UTC'))
label_UTC0 = data['label']

set_graphs.set_matplotlib_param(38)
cmap = 'gist_yarg'

for i,idx in enumerate(idx_x):
    FK = mean_CWT[:,:,idx]
    max_FK = np.max(abs(FK))
    
    fig, ax = plt.subplots()
    imsh = ax.pcolormesh(k,freq,abs(FK)/max_FK,cmap = cmap, shading = 'auto')
    ax.set_ylim([0,1])
    ax.set_xlim([0.01,0.5])
    
    imsh.set_clim([0,1])
    imsh.set_rasterized(True)
    
    k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
    ax.plot(k_exp,f_exp,'r.')
    print(f'Detected coordinates :{k_exp}')
    
    if len(k_exp) != 0:
        H = get_water_height(DAS_water_height,UTC_0,selected_x[i])
        popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                             bounds = (1e5,1e8))
        err_coeff = np.sqrt(np.diag(pcov))
        print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e} J')
    
        k_th = np.linspace(0,1,100)
        omega_th = shallow_hydroelastic(k_th, popt[0], rho_w, H)
        label_th = r'$D = ' + f'{popt[0]:.2e}' + r' \; \mathrm{J}$'
        ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
        # ax.legend(loc = 'upper right',fontsize = 17)
        
    ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
    ax.set_ylabel(r'$f \; \mathrm{(Hz)}$') 
    
    title = r'$x = ' + f'{selected_x[i]:.0f}' + '\; \mathrm{m}$'
    ax.set_title(title)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)
    cbar.set_label(r'$\hat{\dot{\epsilon}}$')
    
    x_txt = f'x_{selected_x[i]:.0f}'
    
    figname = f'{fig_folder}FK_{label_UTC0}_{x_txt}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')