# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 16:22:28 2026

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
import icewave.drone.drone_projection as dp
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

#%% Set fig_folder path 

fig_folder = 'F:/Rimouski_2025/Data/Summary/DAS/Uncertainties/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Function section 

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

def plot_FK_panel(data,selected_x):
    """ Plot FK for different positions into a 3 x 2 subplot panel 
    Inputs : - data, dictionnary
             - selected_x, selected positions in meter for which we want to show FK plot """
             
    # Plot (f,k) for different positions 
    s = data['x']
    k = data['k_star']
    freq = data['f']
    mean_CWT = data['mean_CWT']
    
    idx_x = [np.argmin(abs(s - x)) for x in selected_x]

    # plot FK without theory
    fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
    for ax, idx in zip(axs.flatten(),idx_x):
        imsh = ax.pcolormesh(k,freq,abs(mean_CWT[:,:,idx]),cmap = parula_map)
        ax.set_ylim([0,1])
        ax.set_xlim([0.01,0.6])
        # ax.set_ylim([0,10])
        # ax.set_xlim([0,1.0])
        ax.set_title(f'x = {s[idx]:.2f} m')
        # imsh.set_clim([0,300])

    for j in range(axs.shape[1]):
        axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

    for i in range(axs.shape[0]):
        axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')
    
    return fig, axs

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

def flexural_day(D_results):
    """ Compute mean value of flexural values for each position along optical fiber 
    Inputs : - D_results, dictionnary, containing keys 'D' and 'err_D' which are the flexural modulus and 
    associate error obtained from the fit. 
    
    Outputs : - mean_filtered, array like, mean value of flexural modulus at each position, getting rid off outliers
              - main_std, array like, standard deviation computed using both uncertainty from type A (std/np.sqrt(N))
            and type B (average std over all fits for each position)
    """
    
    keys = list(D_results.keys())
    # create an array of flexural modulus and error
    D_array = np.zeros((len(D_results[keys[0]]['D']),len(keys)))
    errD_array  = np.zeros((len(D_results[keys[0]]['D']),len(keys)))
    for i,key in enumerate(keys) :
        D_array[:,i] = D_results[key]['D']
        errD_array[:,i] = D_results[key]['err_D']
        
    # Compute average of flexural modulus for each position 
    mean_D = np.nanmean(D_array, axis = 1)
    std_D_time = np.nanstd(D_array,axis = 1)
    
    # get rid off outliers 
    filtered_D = np.zeros((D_array.shape))
    for j in range(D_array.shape[1]):
        test = abs(D_array[:,j] - mean_D) > 2*std_D_time 
        for i in range(D_array.shape[0]):
            if test[i] == 1:
                filtered_D[i,j] = None
                print('Outliers detected')
            else : 
                filtered_D[i,j] = D_array[i,j]
        
    filtered_errD = np.zeros((errD_array.shape))
    for j in range(D_array.shape[1]):
        test = abs(D_array[:,j] - mean_D) > 2*std_D_time 
        for i in range(D_array.shape[0]):
            if test[i] == 1:
                filtered_errD[i,j] = None
            else : 
                filtered_errD[i,j] = errD_array[i,j]
            
    mean_filtered = np.nanmean(filtered_D,axis = 1)
    std_filtered = np.nanstd(filtered_D,axis = 1) # standard deviation type A
    
    std_b = np.nanmean(filtered_errD, axis = 1) # standard deviation type b (mean std over each fit for a given position)
    
    main_std = np.sqrt((std_filtered/np.sqrt(filtered_D.shape[0]))**2 + std_b**2)
    
    return mean_filtered, main_std 


# =============================================================================
# %% Load averaged CWT 
# =============================================================================

date = '0211'
# main_path = f'U:'

main_path = 'F:/Rimouski_2025/'
path2data = f'{main_path}Data/{date}/DAS/'

path2CWT = f'{path2data}avg_CWT/'
filelist = glob.glob(f'{path2CWT}avg_CWT*.h5')
print(filelist)

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

# load water height data
file2water = glob.glob(f'{path2data}fiber_water_height_GPS_structure_{date}.h5')[0]
print(file2water)    
DAS_water_height = load_DAS_water_height(file2water)

# # load swell orientation 
# file_swell_orientation = f'{main_path}Data/swell_orientation_beam_forming.h5'
# swell_orientation = rw.load_dict_from_h5(file_swell_orientation)
swell_orientation = 0

#%% Load data from a single file 

file2load = filelist[0]
data = rw.load_dict_from_h5(file2load)

selected_x = np.array([80,150,300,400,500,550]) # 0211
set_graphs.set_matplotlib_param('single')
fig, axs = plot_FK_panel(data,selected_x)

# extract peaks from CWT
mean_CWT = data['mean_CWT']
s = data['x']
freq = data['f']
k = data['k_star']

format_date = '%Y-%m-%dT%H-%M-%S'
UTC_0 = datetime.strptime(data['label'],format_date)
UTC_mid = UTC_0.replace(tzinfo = pytz.timezone('UTC'))

peaks_detec = {}
peaks_detec['freq_range'] = [0.02,0.8]
peaks_detec['min_prominence'] = 0.7
peaks_detec['min_width'] = 3
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

#%% Evaluate model sensitivity to D

def error_model(f_exp,k_exp,best_D,rho_w,H,DELTA_D = 1e7):
    
    # estimate error 
    N = len(f_exp)
    f_model_best = shallow_hydroelastic(k_exp,best_D, rho_w, H)/2/np.pi
    residuals = f_model_best - f_exp
    # écart-type des résidus (1 paramètre ajusté -> dof = N-1)
    dof = max(N - 1, 1)
    sigma_res = np.sqrt(np.sum(residuals**2) / dof)
    
    D_plus = best_D + DELTA_D
    D_minus = best_D - DELTA_D
    actual_delta = (D_plus - D_minus) / 2
    f_model_plus = shallow_hydroelastic(k_exp, D_plus , rho_w, H)/2/np.pi
    f_model_minus = shallow_hydroelastic(k_exp, D_minus , rho_w, H)/2/np.pi
    df_dD = (f_model_plus - f_model_minus) / (2 * actual_delta)
    denom = np.sum(df_dD**2)
    if denom > 0:
        sigma_D_res = sigma_res / np.sqrt(denom)
    else:
        sigma_D_res = 0.0
    
    return sigma_D_res

#%% Plot results for a given x
DELTA_D = 1e7
# estimate error 
N = len(f_exp)
f_model_best = shallow_hydroelastic(k_exp, popt[0] , rho_w, H)/2/np.pi
residuals = f_model_best - f_exp
# écart-type des résidus (1 paramètre ajusté -> dof = N-1)
dof = max(N - 1, 1)
sigma_res = np.sqrt(np.sum(residuals**2) / dof)

D_plus = popt[0] + DELTA_D
D_minus = popt[0] - DELTA_D
actual_delta = (D_plus - D_minus) / 2
f_model_plus = shallow_hydroelastic(k_exp, D_plus , rho_w, H)/2/np.pi
f_model_minus = shallow_hydroelastic(k_exp, D_minus , rho_w, H)/2/np.pi
df_dD = (f_model_plus - f_model_minus) / (2 * actual_delta)
denom = np.sum(df_dD**2)
if denom > 0:
    sigma_D_res = sigma_res / np.sqrt(denom)
else:
    sigma_D_res = 0.0

sigma_D = np.sqrt(sigma_D_res**2 + err_coeff[0]**2)

xth = np.linspace(1e-2,2e-1,100)
yth = shallow_hydroelastic(xth, popt[0], rho_w, H)/2/np.pi
y_low = shallow_hydroelastic(xth, popt[0] - sigma_D, rho_w, H)/2/np.pi
y_up = shallow_hydroelastic(xth, popt[0] + sigma_D, rho_w, H)/2/np.pi


fig, ax = plt.subplots()
ax.plot(k_exp,f_exp,'.',color = 'tab:red')
ax.plot(xth,yth)
ax.plot(xth,y_low)
ax.plot(xth,y_up)



#%% Load flexural modulus data from CWT analysis #0211

main_path = 'F:/Rimouski_2025/Data/'
# date_DAS = ['0211', '0212']
date_DAS = ['0211']

# set offset fiber
offset_fiber = 37.5 # in meters

main_D = {}
main_D['corrected'] = {}
for date in date_DAS : 
    filepath = f'{main_path}{date}/DAS/Figures/Wavelet_study*/{date}_wavelet_flexural_modulus_subpix_swell_corrected_file*.h5'
    filelist = glob.glob(filepath, recursive = True)
    print(filelist)
    
    D_results = {}
    for file2load in filelist:
        data = rw.load_dict_from_h5(file2load)
        UTC_0 = data['UTC_0']
        D_results[UTC_0] = data
    
    keys = list(D_results.keys())
    mean_filtered, main_std = flexural_day(D_results)
    
    main_D['corrected'][date] = {'D':mean_filtered, 'err_D': main_std, 'x' : D_results[keys[0]]['x'] - offset_fiber}

# main_D['uncorrected'] = {}
# for date in date_DAS : 
#     filepath = f'{main_path}{date}/DAS/Figures/Wavelet_study*/{date}_wavelet_flexural_modulus_subpix__file*.h5'
#     filelist = glob.glob(filepath, recursive = True)
#     print(filelist)
    
#     D_results = {}
#     for file2load in filelist:
#         data = rw.load_dict_from_h5(file2load)
#         UTC_0 = data['UTC_0']
#         D_results[UTC_0] = data
    
#     keys = list(D_results.keys())
#     mean_filtered, main_std = flexural_day(D_results)
    
#     main_D['uncorrected'][date] = {'D':mean_filtered, 'err_D': main_std, 'x' : D_results[keys[0]]['x'] - offset_fiber}

#%% 

# fig, ax = plt.subplots()
# for key in D_results.keys():
#     D = D_results[key]['D']
#     err_D = D_results[key]['err_D']
#     x = D_results[key]['x']
    
#     ax.plot(x,D,'o-')
#     ax.fill_between(x,D - err_D,D + err_D, alpha = 0.2)
    
#%%

keys = list(D_results.keys())
# create an array of flexural modulus and error
D_array = np.zeros((len(D_results[keys[0]]['D']),len(keys)))
errD_array  = np.zeros((len(D_results[keys[0]]['D']),len(keys)))

fig, ax = plt.subplots()
for i,key in enumerate(keys) :
    D = D_results[key]['D']
    err_D = D_results[key]['err_D']
    x = D_results[key]['x']
    D_array[:,i] = D_results[key]['D']
    errD_array[:,i] = D_results[key]['err_D']
        
    # Compute average of flexural modulus for each position 
    mean_D = np.nanmean(D_array, axis = 1)
    std_D_time = np.nanstd(D_array,axis = 1)
    
    # # get rid off outliers 
    # filtered_D = np.zeros((D_array.shape))
    # for j in range(D_array.shape[1]):
    #     test = abs(D_array[:,j] - mean_D) > 2*std_D_time 
    #     for i in range(D_array.shape[0]):
    #         if test[i] == 1:
    #             filtered_D[i,j] = None
    #             print('Outliers detected')
    #         else : 
    #             filtered_D[i,j] = D_array[i,j]
        
    # filtered_errD = np.zeros((errD_array.shape))
    # for j in range(D_array.shape[1]):
    #     test = abs(D_array[:,j] - mean_D) > 2*std_D_time 
    #     for i in range(D_array.shape[0]):
    #         if test[i] == 1:
    #             filtered_errD[i,j] = None
    #         else : 
    #             filtered_errD[i,j] = errD_array[i,j]
            
    # mean_filtered = np.nanmean(filtered_D,axis = 1)
    # std_filtered = np.nanstd(filtered_D,axis = 1) # standard deviation type A
    
    std_b = np.nanmean(errD_array, axis = 1) # standard deviation type b (mean std over each fit for a given position)
    
    main_std = np.sqrt((std_D_time/np.sqrt(D_array.shape[1]))**2 + std_b**2)
    
    
    ax.plot(x,D_array[:,i],'o-')
    ax.fill_between(x,D - err_D,D + err_D, alpha = 0.2)

ax.plot(x,mean_D,'o-',color = 'k')
ax.fill_between(x,mean_D - main_std, mean_D + main_std, color = 'k', alpha = 0.2)

