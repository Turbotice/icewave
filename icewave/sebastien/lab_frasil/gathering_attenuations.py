# -*- coding: utf-8 -*-
"""
Created on Tue May  6 10:59:19 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

import pickle
import os
import glob
import re
import cv2 as cv 

from concurrent.futures import ProcessPoolExecutor

import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% FUNCTION SECTION 
def collect_all_experiments(filelist):
    data_list = []
    for file in filelist:
        with open(file,'rb') as pf:
            data = pickle.load(pf)
        data_list.append(data)
        
    return data_list

def affine(x,a,b):
    y = a*x + b
    return y

def powerlaw_fit(x,y,err_y):
    """ Fit data using a power law, taking into account standard deviation of y """
    log_x = np.log(x)
    log_y = np.log(y)
    err_log_y = err_y/y
    
    popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y,sigma = err_log_y,absolute_sigma = True)
    err_affine = np.sqrt(np.diag(pcov))
    beta = popt[0]
    err_beta = err_affine[0]
    B = np.exp(popt[1])
    err_B = B*err_affine[1]
    
    coeffs = (beta,B)
    err_coeffs = (err_beta,err_B)
    return coeffs,err_coeffs

#%% Load data 
h = 7.5
date = '2024_07_11'
main_path = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/'

# data from laser 
path2data_laser = f'{main_path}Laser_attenuation/'
filelist = glob.glob(f'{path2data_laser}**/laser_attenuation*.pkl')
# data from PIV
path2data_PIV = f'{main_path}matData/'
filelist_PIV = glob.glob(f'{path2data_PIV}**/PIV_attenuation_results*.pkl')

data_laser = collect_all_experiments(filelist)
data_PIV = collect_all_experiments(filelist_PIV)


#%% Compare attenuation coefficient

main_dict = {}
main_dict['laser'] = {}
keys = ['alpha','err_alpha','k0','f_demod','f_ex']
for key in keys:
    main_dict['laser'][key] = np.array([data[key] for data in data_laser ])

main_dict['PIV'] = {}
keys = ['alpha','err_alpha','k','f_demod','f_ex']
for key in keys:
    main_dict['PIV'][key] = np.array([data[key] for data in data_PIV ])

#%% Plot data 

power_laws = {}
coeffs,err_coeffs = powerlaw_fit(main_dict['laser']['f_demod'], main_dict['laser']['alpha'], main_dict['laser']['err_alpha'])
power_laws['laser'] = {}
power_laws['laser']['coeffs'] = coeffs
power_laws['laser']['err_coeffs'] = err_coeffs
coeffs,err_coeffs = powerlaw_fit(main_dict['PIV']['f_demod'], main_dict['PIV']['alpha'], main_dict['PIV']['err_alpha'])
power_laws['PIV'] = {}
power_laws['PIV']['coeffs'] = coeffs
power_laws['PIV']['err_coeffs'] = err_coeffs

yth = {}
xth = np.linspace(1e0,1e1,200)
for key in ['laser','PIV']:
    yth[key] = power_laws[key]['coeffs'][1]*xth**power_laws[key]['coeffs'][0]


fig, ax = plt.subplots(figsize = (12,9))
ax.plot(main_dict['laser']['f_demod'],main_dict['laser']['alpha'],'o',color = 'tab:blue',label = 'laser')
coeffs = power_laws['laser']['coeffs']
label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
ax.plot(xth,yth['laser'],'-',color = 'tab:blue',label = label_th)
ax.plot(main_dict['PIV']['f_demod'],main_dict['PIV']['alpha'],'d',color = 'tab:orange',label = 'PIV')
coeffs = power_laws['PIV']['coeffs']
label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
ax.plot(xth,yth['PIV'],'-',color = 'tab:orange',label = label_th)
ax.legend()

ax.set_xlabel(r'$f_{\mathrm{demod}} \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2e0,8e0])
ax.set_ylim([1e0,1e2])






