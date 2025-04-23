# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 09:28:57 2025

@author: sebas

This script aims at gathering all results collected from drones Data

"""

#%% Import modules

import numpy as np
import os
import glob
import pickle

import matplotlib.pyplot as plt

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.Fourier_tools as FT
import icewave.multi_instruments.results as res
import icewave.tools.datafolders as df

#%% Function section 

def get_amp_freq_from_drone_results(drone_results):
    TF_spectrum = drone_results['FFT_spectrum']['TF_spectrum']
    f = drone_results['FFT_spectrum']['f']

    idx_max = np.argmax(TF_spectrum)
    Amax,f0 = FT.subpix_precision_array(TF_spectrum, f, idx_max)
    
    return Amax,f0

#%% Load an exemple of results file 


date = '0226'
typ = 'phones'
main_path = 'K:/Share_hublot/Data/'
path2results = f'{main_path}{date}/Summary/results_{typ}_{date}.pkl'

with open(path2results,'rb') as pf:
    results = pickle.load(pf)
    
#%% Load results associated to a single drone movie
date = '0226'
typ = 'Drones'
drone_ID = 'mesange'
exp_ID = '23-waves_012'
path2drone_results = f'{main_path}{date}/{typ}/{drone_ID}/matData/{exp_ID}/Results/'

file2load = f'{path2drone_results}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2load,'rb') as pf:
    drone_results = pickle.load(pf)

#%% Detect maximum amplitude and associated frequency of TF spectrum 

Amax,f0 = get_amp_freq_from_drone_results(drone_results)

#%% Load records for corresponding date 

file2load = f'{main_path}{date}/Summary/records_{date}.pkl'
with open(file2load,'rb') as pf:
    records = pickle.load(pf)
    
#%%
# write Amax and f0 in a results dictionnary
main_results = {}
result = res.make_result(date, 'drones', drone_ID, exp_ID, 'f0', f0, comments = 'unit : Hz')
main_results.update(result)
result = res.make_result(date, 'drones', drone_ID, exp_ID, 'Amax', Amax)
main_results.update(result)

# write A,f,k,alpha 
keys2write = ['A','f','k','alpha','hw']
for key in keys2write:
    value = drone_results['attenuation'][key]
    result = res.make_result(date, 'drones', drone_ID, exp_ID, key, value)
    main_results.update(result)

# write powerlaw attenuation alpha = B*f**beta
keys2write = ['B','beta']
for key in keys2write:
    value = drone_results['attenuation']['power_law'][key]
    result = res.make_result(date, 'drones', drone_ID, exp_ID, key, value)
    main_results.update(result)





















