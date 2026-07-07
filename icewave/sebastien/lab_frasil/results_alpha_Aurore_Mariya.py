# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 09:32:46 2026

@author: sebas

Gather attenuation results from Aurore and Mariya experiments
"""


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

def func_real_imag(z,dimensionless):
    x,y = z # z = [Re(z),Im(z)]
    
    kappa = x + 1j*y
    determinant = theory.det_M(kappa,dimensionless)
    return [np.real(determinant),np.imag(determinant)] # system of 2 real equations 

def get_f_list(data,h,amplitude):
    """ Get frequency array for a given set of experiments"""
    
    valid_keys = []
    for key in data.keys():
        if data[key]['h'] == h and data[key]['amplitude'] == amplitude:
            valid_keys.append(key)

    f_list = [data[key]['f_ex'] for key in valid_keys]
    f_list = np.array(f_list)
    return f_list

#%% Load data 

main_dict = {}

file2load = 'U:/Aurore_frasil/attenuation_results.pkl'
with open(file2load,'rb') as pf:
    main_dict['aurore'] = pickle.load(pf)
    
file2load = 'F:/Stagiaires/Mariya/attenuation_results.pkl'
with open(file2load,'rb') as pf:
    main_dict['mariya'] = pickle.load(pf)
    
#%% Set fig_folder 

fig_folder = 'F:/PhD_Manuscript/ch3/Attenuation/'   
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)

#%% Plot Aurore results 

data = main_dict['aurore']
mksize = 12

set_graphs.set_matplotlib_param('powerpoint')
method = 'laser'
thickness_list = [5.0,7.5,10.0,12.5,15.0]

bounds = np.array([3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots(figsize = (12,9))

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])

    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,
                markeredgecolor = 'k',ecolor = 'k',
                markersize = mksize)
    
    xth = np.linspace(2.0,6.0,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    ax.plot(xth,yth,'-',color = current_color,label = label_th)
    
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2.7,5.3])
ax.set_ylim([2e0,70])


# =============================================================================
# %% Plot Mariya results    
# =============================================================================
    
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

data = main_dict['mariya']

full_reds = mpl.colormaps['Reds'].resampled(256)
new_reds = mcolors.ListedColormap(full_reds(np.linspace(0.2,0.8,256)))
cnorm = mcolors.Normalize(vmin = 8.0,vmax = 22.0)
cmap = new_reds

thickness_list = [8.0,17.0,22.0]

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    
    current_color = new_reds(cnorm(h))
    marker = 's'
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,
                markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
    xth = np.linspace(1,1e1,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    # ax.plot(xth,yth,'-',color = current_color,label = label_th)
    
# ax.legend()

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
# sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1,1e1])
ax.set_ylim([1,1e2])


figname = f'{fig_folder}alpha_VS_freq_Mariya_rhow_1p07_laser_lorentz'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')


# =============================================================================
# %% Average Mariya's results over motor amplitude + superposition with Aurore results
# =============================================================================

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

# set markersize
ms = 6

data = main_dict['aurore']
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
cnorm = mcolors.Normalize(2.5,22.0)
cmap = new_blues
# marker_list = ['o','s','D','^','h','p']

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    
    current_color = cmap(cnorm(h))
    marker = 'o'
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = ms)
    
    xth = np.linspace(1.0,1e1,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    # ax.plot(xth,yth,'-',color = current_color,label = label_th)
    

data = main_dict['mariya']

a = 15.0
thickness_list = [8.0,17.0,22.0]

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['h'] == h:
            valid_keys.append(key)
    
    f_list = get_f_list(data,h,a)
    
    alpha = []
    err_alpha = []
    f_demod = []
    for f_ex in f_list:
        keys_fex = []
        for key in valid_keys:
            if data[key]['f_ex'] == f_ex:
                keys_fex.append(key)
        
        alpha.append(np.mean([data[key]['alpha'] for key in keys_fex]))
        err_alpha.append(np.std([data[key]['err_alpha'] for key in keys_fex]))
        f_demod.append(np.mean([data[key]['f_demod'] for key in keys_fex]))
    
    alpha = np.array(alpha)
    err_alpha = np.array(err_alpha)
    f_demod = np.array(f_demod)
    
    current_color = cmap(cnorm(h))
    marker = 's'
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'firebrick',
                ecolor = 'k',markeredgewidth = 1.1,
                markersize = ms)
    
    xth = np.linspace(1,1e1,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    # ax.plot(xth,yth,'-',color = current_color,label = label_th)


divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
# sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1.8,6])
ax.set_ylim([3,1e2])

figname = f'{fig_folder}alpha_VS_freq_Mariya_Aurore_superposition'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')


#%% Superpose results 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots(figsize = (12,9))

# plot aurore's data
data = main_dict['aurore']
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
# marker_list = ['o','s','D','^','h','p']

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    
    current_color = cmap(norm(h))
    marker = '-'
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
    xth = np.linspace(2.0,6.0,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    # ax.plot(xth,yth,'-',color = current_color,label = label_th)
    
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')


# plot Mariya's data

color_list = ['lightcoral','indianred','firebrick']
data = main_dict['mariya']
thickness_list = [8.0,17.0,22.0]

for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    
    current_color = color_list[i]
    marker = '-'
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
    xth = np.linspace(2.0,6.0,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    # ax.plot(xth,yth,'-',color = current_color,label = label_th)


# ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2.7,5.3])
ax.set_ylim([2e0,70])




