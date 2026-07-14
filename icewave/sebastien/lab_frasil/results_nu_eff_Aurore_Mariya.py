# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:49:38 2026

@author: sebas
"""

#%%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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
import icewave.sebastien.set_graphs as set_graphs
import icewave.sebastien.theory.module_bilayer_viscous_dimensionless as theory
import icewave.tools.rw_data as rw

global g
g = 9.81

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 

def affine(x,a,b):
    y = a*x + b
    return y

def powerlaw_fit(x,y,err_y = None):
    """ Fit data using a power law, taking into account standard deviation of y """
    log_x = np.log(x)
    log_y = np.log(y)
    
    if err_y is None :
        popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y)
    
    else : 
        
        err_log_y = err_y/y
        popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y,sigma = err_log_y,
                                             absolute_sigma = True)
        
    err_affine = np.sqrt(np.diag(pcov))
    beta = popt[0]
    err_beta = err_affine[0]
    B = np.exp(popt[1])
    err_B = B*err_affine[1]
    
    coeffs = (beta,B)
    err_coeffs = (err_beta,err_B)
    return coeffs,err_coeffs

def get_wavevector_theory(freq,h,rho_1,rho_2,nu_1,nu_2):
    """ Compute wavector using viscous bilayers, for a given set of frequencies and 
    fluid physical properties. """
    
    params_array = np.array([(f_ex,h,rho_1,rho_2,nu_1,nu_2) for f_ex in freq])
    dimensionless_array = np.array([theory.dim_numbers(params) for params in params_array]) 
    omega_array = 2*np.pi*freq
    k0_array = omega_array**2/g
    
    initial_guess = [1 , 0.1] # initial guess for [Re(kappa), Im(kappa)]
    k_array = np.zeros(len(freq),dtype = 'complex')
    for i,dimensionless in enumerate(dimensionless_array): # loop over dimensionless parameters
        
        solution = scipy.optimize.fsolve(theory.func_real_imag,initial_guess,(dimensionless,))
        root = solution[0] + 1j*solution[1]
        print(f'Root of det(M) = {root}')
        initial_guess = solution
        k_array[i] = root*k0_array[i]
        
    return k_array

#%% Set fig folder

fig_folder = 'F:/PhD_Manuscript/ch3/Attenuation/'   
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)
    
#%% Load effective viscosities from Aurore and Mariya experiments

main_table = {'aurore':{},'mariya':{}}

file2effective_table = 'U:/Aurore_frasil/effective_nu_Lamb.h5'
print(f'{file2effective_table} already exists, loading...')
main_table['aurore']['lamb'] = rw.load_dict_from_h5(file2effective_table)

file2effective_bilayer = 'U:/Aurore_frasil/effective_nu_viscous_bilayer_imag.h5'
print(f'{file2effective_bilayer} already exists, loading...')
main_table['aurore']['bilayer'] = rw.load_dict_from_h5(file2effective_bilayer)

file2effective_table = 'F:/Stagiaires/Mariya/effective_nu_Lamb.h5'
main_table['mariya']['lamb'] = rw.load_dict_from_h5(file2effective_table)

file2effective_bilayer = 'F:/Stagiaires/Mariya/effective_nu_viscous_bilayer_imag.h5'
main_table['mariya']['bilayer'] = rw.load_dict_from_h5(file2effective_bilayer)


#%% Plot results

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
nu_w = 1e-6

color_table = {'aurore':'tab:blue','mariya':'tab:red'}
marker_table = {'lamb':'o','bilayer':'^'}
ms = 8

method = 'laser'
for key in main_table.keys():
    for key_model in main_table[key].keys():
        
        if key == 'mariya' and key_model == 'bilayer':
            table = main_table[key][key_model]

        else:
            table = main_table[key][key_model][method]

        label = f'{key} - {key_model}'
        marker = marker_table[key_model]
        ax.errorbar(table['h'],table['nu']/nu_w,
                    yerr = table['err_nu']/nu_w,fmt = marker,color = color_table[key],
                    label = label, markersize = ms)
                
    
ax.set_xlabel(r'$h \; \mathrm{(mm)}$')
ax.set_ylabel(r'$\nu / \nu_w$')
ax.legend()

figname = f'{fig_folder}nu_eff_VS_h_rho_1p07_rho_1p0'
# plt.savefig(f'{figname}.png', bbox_inches = 'tight')
# plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
# plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Plot only Aurore results 


set_graphs.set_matplotlib_param('powerpoint')
fig, ax = plt.subplots()
nu_w = 1e-6

color_table = {'aurore':'tab:blue','mariya':'tab:red'}
marker_table = {'lamb':'o','bilayer':'^'}
ms = 12

method = 'laser'
key = 'aurore'
key_model = 'lamb'
table = main_table[key][key_model][method]

marker = marker_table[key_model]
ax.errorbar(table['h'],table['nu']/nu_w,
            yerr = table['err_nu']/nu_w,fmt = marker,color = color_table[key],
            markersize = ms)
    
ax.set_xlabel(r'$h \; \mathrm{(mm)}$')
ax.set_ylabel(r'$\nu / \nu_w$')

ax.set_xlim([0,16])
# ax.set_yscale('log')

figname = f'{fig_folder}nu_eff_VS_h_Aurore_data_Lamb_model'
# plt.savefig(f'{figname}.png', bbox_inches = 'tight')
# plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
# plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

