# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 15:42:50 2025

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

fig_folder = 'U:/Data/0211/DAS/Figures_article/CWT/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Function section 

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

#%% Load flexural modulus data from CWT analysis #0211

main_path = 'U:/Data/'
date_DAS = ['0211', '0212']


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
    
    main_D['corrected'][date] = {'D':mean_filtered, 'err_D': main_std, 'x' : D_results[keys[0]]['x']}

main_D['uncorrected'] = {}
for date in date_DAS : 
    filepath = f'{main_path}{date}/DAS/Figures/Wavelet_study*/{date}_wavelet_flexural_modulus_subpix__file*.h5'
    filelist = glob.glob(filepath, recursive = True)
    print(filelist)
    
    D_results = {}
    for file2load in filelist:
        data = rw.load_dict_from_h5(file2load)
        UTC_0 = data['UTC_0']
        D_results[UTC_0] = data
    
    keys = list(D_results.keys())
    mean_filtered, main_std = flexural_day(D_results)
    
    main_D['uncorrected'][date] = {'D':mean_filtered, 'err_D': main_std, 'x' : D_results[keys[0]]['x']}

#%% Plot mean flexural modulus

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots(figsize = (12,9))
for key in main_D.keys():
    for date in date_DAS :
        label = f'{key}_{date}'
        ax.errorbar(main_D[key][date]['x'],main_D[key][date]['D'],yerr = main_D[key][date]['err_D'],fmt = '-o',
                    label = label)
ax.legend()

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')

figname = f'{fig_folder}D_VS_x_swell_correction_comparison'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Load Ludo's results 

path2active_results = f'{main_path}Summary/DAS/Results_active_sources/'

# load '0211' data
date = '0211'
path2values_h = os.path.join(path2active_results,f'thicknesses{date}.pkl')
path2values_E = os.path.join(path2active_results, 'young_modulus.pkl')

with open(path2values_h, 'rb') as f:
    thicknesses = pickle.load(f)
    print("Contents of thicknesses.pkl:")
    print(thicknesses)

with open(path2values_E, 'rb') as f:
    young_modulus = pickle.load(f)
    print("\nContents of young_modulus.pkl:")
    print(young_modulus)

active = {date:{}}
active[date]['h'] = thicknesses['thickness']
active[date]['x'] = thicknesses['position']

#%% load 0212 data
date = '0212'

path2values_h = os.path.join(path2active_results, f'thicknesses{date}.pkl')
path2values_E = os.path.join(path2active_results, 'young_modulus.pkl')

with open(path2values_h, 'rb') as f:
    thicknesses = pickle.load(f)
    print("Contents of thicknesses.pkl:")
    print(thicknesses)

with open(path2values_E, 'rb') as f:
    young_modulus = pickle.load(f)
    print("\nContents of young_modulus.pkl:")
    print(young_modulus)


# --- Prepare h data ---
h_raw = np.array(thicknesses['h'])
h_xpos_raw = np.array(list(thicknesses['xposition'].values()))

# Compute moving average (window = 4)
window = 2
h_ma = np.convolve(h_raw, np.ones(window)/window, mode='valid')

# Adjust x-position for moving average (centered)
h_xpos_ma = h_xpos_raw[(window - 1)//2 : -(window//2)] if len(h_xpos_raw) >= window else []

# --- Prepare E data ---
E_raw = np.array(young_modulus['E'][:-1])  # drop last value
E_xpos_raw = np.array(list(young_modulus['xposition'].values())[:-1])

E_ma = np.convolve(E_raw, np.ones(window)/window, mode='valid')
E_xpos_ma = E_xpos_raw[(window - 1)//2 : -(window//2)] if len(E_xpos_raw) >= window else []

nu = 0.3

# --- Interpolate E onto h_xpos_raw ---
E_interp = np.interp(h_xpos_raw, E_xpos_raw, E_raw)

# Compute D on h_xpos_raw grid
D_raw = E_interp * (h_raw**3) / (12 * (1 - nu**2))

# eventual outliers 
outliers = np.array([])
mask = []
for pos in h_xpos_raw :
    mask.append(pos not in outliers)

# save data in dictionnary
active[date] = {}
active[date]['E'] = E_interp[mask]
active[date]['h'] = h_raw[mask]
active[date]['x'] = h_xpos_raw[mask]
active[date]['D'] = D_raw[mask]

#%% Plot Ludo's results and passive results 

# create mask for hidding data
# outliers = np.array([236,256])
outliers = np.array([])
mask = []
for pos in h_xpos_raw :
    mask.append(pos not in outliers)


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

color_date = {'0211':new_blues(0.5), '0212': new_blues(0.8)}
key_swell = 'corrected'
for date in date_DAS :
    ax.errorbar(main_D[key_swell][date]['x'],main_D[key_swell][date]['D'],yerr = main_D[key_swell][date]['err_D'],
                fmt = '-o',label = date, color = color_date[date])

ax.plot(h_xpos_raw[mask],D_raw[mask],'o-',color = 'tab:red',label = 'Active 0212')

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')
ax.set_ylim([4e6,2e8])
ax.set_yscale('log')

ax.legend(loc = 'lower right')

figname = f'{fig_folder}D_VS_x_comparison_with_active_sources_swell_corrected_log_scale'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Create a dictionnary that gathers all results

# load DAS gps positions
file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

# GPS coordinates of beginning and end of fiber
Lat0 = DAS_water_height['lat'][0]
Long0 = DAS_water_height['long'][0]

Lat1 = DAS_water_height['lat'][-1]
Long1 = DAS_water_height['long'][-1]

# azimuth angle of fiber 
psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(Long1 - Long0)*np.pi/180,(Lat1 - Lat0)*np.pi/180)*180/np.pi

# create a dictionnary
results = {}
results['passive'] = main_D
    
# store active results
results['active'] = {}
for key_date in ['0211','0212']:
    results['active'][key_date] = active[key_date]
    
# compute GPS positions
for key_corr in results['passive'].keys():
    for key_date in results['passive'][key_corr].keys():
        
        Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0, Long0, psi, results['passive'][key_corr][key_date]['x'])
        results['passive'][key_corr][key_date]['latitude'] = Lat
        results['passive'][key_corr][key_date]['longitude'] = Long


for key_date in results['active'].keys():
    
    Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0, Long0, psi, results['active'][key_date]['x'])
    results['active'][key_date]['latitude'] = Lat
    results['active'][key_date]['longitude'] = Long
    

filename = 'U:/Data/Summary/DAS/main_results_active_passive_V2.h5'
rw.save_dict_to_h5(results, filename)

filename = filename.replace('.h5','.pkl')
with open(filename,'wb') as pf:
    pickle.dump(results,pf)

#%%

file_swell_orientation = f'{main_path}swell_orientation_beam_forming.h5'
swell_orientation = rw.load_dict_from_h5(file_swell_orientation)



