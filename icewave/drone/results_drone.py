# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 09:28:57 2025

@author: sebas

This script aims at gathering all results collected from drones Data

"""

# Import modules

import numpy as np
import os
import glob
import pickle

import icewave.tools.Fourier_tools as FT
import icewave.multi_instruments.results as res
import icewave.tools.datafolders as df
import argparse


# Function section 
# def gen_parser():
#     """ Generate parser to use this script in command line (see doc of argparse module for more details) 
#     Output : - args, object containing a single attribute 'date', which is the date for which we want to gather all results """
#     parser = argparse.ArgumentParser(description="Gather Drone results data")
#     parser.add_argument('-date', dest='date', type=str, default='0226',help='select date for which drone results are gathered')
    
#     args = parser.parse_args()
#     print(args)
#     return args

def get_amp_freq_from_drone_results(drone_results):
    TF_spectrum = drone_results['FFT_spectrum']['TF_spectrum']
    f = drone_results['FFT_spectrum']['f']

    idx_max = np.argmax(TF_spectrum)
    Amax,f0 = FT.subpix_precision_array(TF_spectrum, f, idx_max)
    
    return Amax,f0

def get_drone_results_from_dict(drone_results):
    
    date = drone_results['date']
    drone_ID = drone_results['drone_ID']
    exp_ID = drone_results['exp_ID']

    main_results = {}
    
    # Detect maximum amplitude and associated frequency of TF spectrum
    Amax,f0 = get_amp_freq_from_drone_results(drone_results)

    # write Amax and f0 in a results dictionnary
    result = res.make_result(date, 'drones', drone_ID, exp_ID, 'f0', f0)
    main_results.update(result)
    result = res.make_result(date, 'drones', drone_ID, exp_ID, 'Amax', Amax)
    main_results.update(result)

    # write A,f,k,alpha 
    attenuation_dict = drone_results['attenuation']
    keys2write = ['A','f','k','alpha','hw']
    for key in keys2write:
        value = attenuation_dict[key]
        result = res.make_result(date, 'drones', drone_ID, exp_ID, key, value)
        main_results.update(result)

    # write powerlaw attenuation alpha = B*f**beta
    keys2write = ['B','beta']
    for key in keys2write:
        value = attenuation_dict['power_law'][key]
        result = res.make_result(date, 'drones', drone_ID, exp_ID, key, value)
        main_results.update(result)
    
    return main_results

def collect_drone_results(date,disk = 'Elements',year = '2024'):
    
    # search all results files for given date
    main_path = df.find_path(disk,year)
    path2drone = f'{main_path}{date}/Drones/'
    filelist = glob.glob(path2drone + '**/Results/main_results*.pkl',recursive = True)
    print(filelist)
    
    main_results = {}
    
    for file2load in filelist:
        with open(file2load,'rb') as pf:
            drone_results = pickle.load(pf)
        
        current_results = get_drone_results_from_dict(drone_results)
        main_results.update(current_results)
        
        print(f'Drone results updated with file {file2load}')
    
        res.save_result(main_results)
        print('Drone results saved')

# if code is used a main script from command line 
if __name__ == '__main__':
    # args = gen_parser()
    date = '0226'
    collect_drone_results(date)

