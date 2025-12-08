# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 16:51:11 2025

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
import imageio as iio
import cv2 as cv
import h5py
import csv
import scipy

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

global g
g = 9.81
global offset_fiber
offset_fiber = 37.5
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10
#%% FUNCTION SECTION 

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega

#-----------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------

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


#%% Load active (f,k)

date = '0211'
base = 'U:/Data/'
path2data = f'{base}{date}/DAS/'
save_file = f'{path2data}FK_0211_active.pkl'

with open(save_file, 'rb') as fid:
   data = pickle.load(fid)


#%% Load DAS_water_height
file2water = glob.glob(f'{path2data}fiber_water_height_GPS_structure_{date}.h5')[0]
print(file2water)    
DAS_water_height = load_DAS_water_height(file2water)

#%% Load avg CWT 

path2CWT = f'{path2data}avg_CWT/'
filelist = glob.glob(f'{path2CWT}avg_CWT*.h5')
print(filelist)

# load swell orientation 
file_swell_orientation = 'U:/Data/swell_orientation_beam_forming.h5'
swell_orientation = rw.load_dict_from_h5(file_swell_orientation)

swell_bool = 1 #1 if we want to account for swell orientation correction 
swell_dict =  {'bool':swell_bool,'theta':swell_orientation[date]}

idx_file = 5
file2load = filelist[idx_file]
print(file2load)
passive = rw.load_dict_from_h5(file2load)

#%% Set fig_folder 

# Create folder for saving graphs
fig_folder = f'{path2data}Figures_article/FK_active_passive/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

# =============================================================================
#%% Subpixel detection active method
# =============================================================================

# Perform subpix extraction 

min_prominence = 0.4
min_width = 20

freq_range = [0.8,32]
k_exp = {}
f_exp = {}

k_exp['active'],f_exp['active'] = extract_peaks(data['FK'].T,data['freq'],data['k'],freq_range,min_prominence,min_width)

#%% Plot FK and detected points 

selected_x = 80.3 - offset_fiber
format_date = '%Y-%m-%d_%H-%M-%S'
label_UTC_0 = '2025-02-11_18-32-43'
UTC_0 = datetime.strptime(label_UTC_0,format_date)
UTC_0 = UTC_0 + timedelta(seconds = 356.1)
UTC_0 = UTC_0.replace(tzinfo = pytz.utc)

k = data['k']
freq = data['freq']
rho_w = 1027
xth = np.linspace(0,1.8,300)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.imshow(data['FK'].T,aspect='auto', extent=[k[0], k[-1], freq[-1], freq[0]], cmap='gray_r')
ax.plot(k_exp['active'],f_exp['active'],'r.',ms = 2)

# Perform fit
H = get_water_height(DAS_water_height,UTC_0,selected_x)
popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp['active'],
                                     f_exp['active'],bounds = (1e5,1e8))
err_coeff = np.sqrt(np.diag(pcov))
yth = shallow_hydroelastic(xth, popt[0], rho_w, H)/2/np.pi
ax.plot(xth,yth,'r--')

ax.invert_yaxis()
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_xlim([0,k[-1]])
ax.set_ylim([0,freq[-1]])


# =============================================================================
#%% Subpixel detection passive method
# =============================================================================

s = passive['x'] - offset_fiber
k = passive['k_star']
freq = passive['f']
mean_CWT = passive['mean_CWT']

selected_x = 80 - offset_fiber # x value 
idx = np.argmin(abs(s - selected_x))
freq_range = [0.02,0.8]
min_prominence = 0.7
min_width = 3

rho_w = 1027

format_date = '%Y-%m-%dT%H-%M-%S'
UTC_0 = datetime.strptime(passive['label'],format_date)
UTC_0 = UTC_0.replace(tzinfo = pytz.timezone('UTC'))
label_UTC0 = passive['label']

set_graphs.set_matplotlib_param('square')
cmap = parula_map
fig, ax = plt.subplots()
FK = mean_CWT[:,:,idx]
max_FK = np.max(abs(FK))
imsh = ax.pcolormesh(k,freq,abs(FK)/max_FK,cmap = cmap, shading = 'auto')
# imsh = ax.imshow(abs(FK).T/max_FK, cmap = parula_map, interpolation = 'auto', origin = 'lower',
#                  extent = [k.min(),k.max(),freq.min(),freq.max()])
ax.set_ylim([0,1])
ax.set_xlim([0.01,0.6])
# ax.set_ylim([0,8])
# ax.set_xlim([0,0.7])
# ax.set_title(f'x = {s[idx]:.0f} m')
imsh.set_clim([0,1])
imsh.set_rasterized(True)

# extract subpixel positions
k_exp['passive'],f_exp['passive'] = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
ax.plot(k_exp['passive'],f_exp['passive'],'r.')
print(f'Detected coordinates :{k_exp}')

if len(k_exp['passive']) != 0:
    H = get_water_height(DAS_water_height,UTC_0,selected_x)
    popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,
                                         k_exp['passive'],f_exp['passive'],bounds = (1e5,1e8))
    err_coeff = np.sqrt(np.diag(pcov))
    print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e} J')

    k_th = np.linspace(0,1,100)
    omega_th = shallow_hydroelastic(k_th, popt[0], rho_w, H)
    label_th = r'$D = ' + f'{popt[0]:.2e}' + r' \; \mathrm{J}$'
    ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
    ax.legend(loc = 'upper right',fontsize = 17)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')   

# fig.colorbar(imsh, ax=axs.ravel().tolist(), pad = 0.02)


#%% Create figure with 3 panels 

rho_w = 1027
scale = 'log'
corresp_axs = ['active','passive','combination']
D = {}
set_graphs.set_matplotlib_param('triple')
fig, axs = plt.subplots(nrows = 1,ncols = 3,figsize = (16,6.5),layout = 'constrained')

for i,ax in enumerate(axs):
    
    # show active FK and detection 
    if corresp_axs[i] == 'active':

        selected_x = 80.3 - offset_fiber
        format_date = '%Y-%m-%d_%H-%M-%S'
        label_UTC_0 = '2025-02-11_18-32-43'
        UTC_0 = datetime.strptime(label_UTC_0,format_date)
        UTC_0 = UTC_0 + timedelta(seconds = 356.1)
        UTC_0 = UTC_0.replace(tzinfo = pytz.utc)
        
        k = data['k']
        freq = data['freq']

        xth = np.linspace(0,1.8,300)
        
        # Perform fit
        H = get_water_height(DAS_water_height,UTC_0,selected_x)
        popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp['active'],
                                             f_exp['active'],bounds = (1e5,1e8))
        err_coeff = np.sqrt(np.diag(pcov))
        print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e} J')
        D['active'] = popt[0]
        yth = shallow_hydroelastic(xth, popt[0], rho_w, H)/2/np.pi
        
        ax.imshow(data['FK'].T/np.max(abs(data['FK'])),aspect='auto', extent=[k[0], k[-1], freq[-1], freq[0]], cmap='gray_r')
        ax.plot(k_exp['active'],f_exp['active'],'.',color = 'tab:blue',ms = 3)
        # ax.plot(xth,yth,'r--')
        
        ax.invert_yaxis()
        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
        ax.set_xlim([0,k[-1]])
        ax.set_ylim([0,freq[-1]])
    
    
    elif corresp_axs[i] == 'passive':
        
        s = passive['x'] - offset_fiber
        k = passive['k_star']
        freq = passive['f']
        mean_CWT = passive['mean_CWT']
    
        selected_x = 80 - offset_fiber # x value 
        idx = np.argmin(abs(s - selected_x))
        freq_range = [0.02,0.8]
        min_prominence = 0.7
        min_width = 3
        
        # select range of k for which we look for a peak
        k_min = 0.02
        k_max = 0.7
        idx_kmin = np.argmin(abs(k - k_min))
        idx_kmax = np.argmin(abs(k - k_max))
        mean_CWT = mean_CWT[:,idx_kmax:idx_kmin,:]
        k = k[idx_kmax:idx_kmin]
    
        format_date = '%Y-%m-%dT%H-%M-%S'
        UTC_0 = datetime.strptime(passive['label'],format_date)
        UTC_0 = UTC_0.replace(tzinfo = pytz.timezone('UTC'))
        label_UTC0 = passive['label']
    
        FK = mean_CWT[:,:,idx]
        max_FK = np.max(abs(FK))
        imsh = ax.pcolormesh(k,freq,abs(FK)/max_FK,cmap = 'gray_r', shading = 'auto')
        # imsh = ax.imshow(abs(FK).T/max_FK, cmap = parula_map, interpolation = 'auto', origin = 'lower',
        #                  extent = [k.min(),k.max(),freq.min(),freq.max()])
        ax.set_ylim([0,1])
        ax.set_xlim([0.02,0.5])
        imsh.set_clim([0,1])
        imsh.set_rasterized(True)
    
        # extract subpixel positions
        k_exp['passive'],f_exp['passive'] = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
        ax.plot(k_exp['passive'],f_exp['passive'],'.',color = 'tab:orange')
        print(f'Detected coordinates :{k_exp}')
    
        H = get_water_height(DAS_water_height,UTC_0,selected_x)
        popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,
                                             k_exp['passive'],f_exp['passive'],bounds = (1e5,1e8))
        err_coeff = np.sqrt(np.diag(pcov))
        print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e} J')
    
        D['passive'] = popt[0]
        #     k_th = np.linspace(0,1,100)
        #     omega_th = shallow_hydroelastic(k_th, popt[0], rho_w, H)
        #     label_th = r'$D = ' + f'{popt[0]:.2e}' + r' \; \mathrm{J}$'
        #     ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
        #     ax.legend(loc = 'upper right',fontsize = 17)
    
        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

    elif corresp_axs[i] == 'combination':
        
        xth = np.linspace(0,2,200)
        yth = shallow_hydroelastic(xth, D['active'], rho_w, H)/2/np.pi
        for key in k_exp.keys():
            ax.plot(k_exp[key],f_exp[key],'.',ms = 8)
        label_th = 'Hydro-elastic'
        ax.plot(xth,yth,'r--',lw = 2,label = label_th)
            
        # ax.legend()
        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

        if scale == 'log':
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([4e-2,2])
            ax.set_ylim([3e-2,60])
        
        elif scale == 'linear':
            ax.set_xlim([0,1.6])
            ax.set_ylim([-1.5,35])

figname = f'{fig_folder}FK_active_passive_x_80_scale_{scale}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')



















#%% Show subpixel extraction for a single frequency 

k_exp = []
f_exp = []

# select a frequency
selected_freq = 10
idx = np.argmin(abs(data['freq'] - selected_freq))

current_max = np.max(abs(data['FK'][:,idx]))
normalized =abs(data['FK'][:,idx])/current_max
# find peaks
peaks,properties = scipy.signal.find_peaks(normalized,
                                           prominence = min_prominence,width = min_width)

if len(peaks) != 0:
    for peak in peaks:
        local_argmax,local_max = subpixel_polyfit(data['k'],normalized,peak)
        k_exp.append(local_argmax)
        f_exp.append(data['freq'][idx])

fig, ax = plt.subplots()
ax.plot(data['k'],normalized,'o-')
ax.plot(local_argmax,local_max,'r.')



