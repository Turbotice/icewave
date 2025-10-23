# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 10:11:37 2025

@author: sebas

This script gathers all experimental results from laboratory experiment. 
An effective viscosity is extracted from the experimental points using either Lamb viscous surface dissipation model
or bi-layer viscous dissipation model (De Carolis et al.)

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

#%% Set figure folder 



#%% Load experimental results

file2load = 'U:/Aurore_frasil/attenuation_results.pkl'
with open(file2load,'rb') as pf:
    data = pickle.load(pf)


#%% Plot attenuation coeff for a single method and all thicknesses

set_graphs.set_matplotlib_param('single')
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots(figsize = (12,9))

effective_table = {}
initial_array = np.zeros((len(thickness_list),1))
effective_table[method] = {'h':initial_array,'nu':initial_array,'prefactor':initial_array,
                           'err_nu':initial_array,'err_prefactor':initial_array}

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
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
    xth = np.linspace(2.0,6.0,100)
    
    # power law fit 
    # coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    # yth = coeffs[1]*xth**coeffs[0]
    # label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    # print(label_th)
    
    # fit by power law omega**3.5
    beta = 3.5
    popt,pcov = scipy.optimize.curve_fit(lambda x,b : affine(x, beta, b),np.log(f_demod*2*np.pi),np.log(alpha),
                                         sigma = err_alpha/alpha,absolute_sigma = True,
                                         bounds = (np.log(1e-5),np.log(1e-1)))
    print(popt)
    coeff = popt[0]
    err_coeff = np.sqrt(np.diag(pcov))[0]
    yth = np.exp(affine(np.log(xth*2*np.pi),beta,coeff))
    
    effective_table[method]['h'][i] = h
    effective_table[method]['prefactor'][i] = np.exp(coeff)
    effective_table[method]['err_prefactor'][i] = np.exp(coeff)*err_coeff
    
    nu = 2*pow(np.exp(coeff),2)*pow(g,4)
    err_nu = 4*pow(g,4)*np.exp(coeff)*np.exp(coeff)*err_coeff
    effective_table[method]['nu'][i] = nu
    effective_table[method]['err_nu'][i] = err_nu
    
    # label_th = f'nu = {nu}'
    label_th = r'$\nu_{eff} = ' + f'{nu*1e6:.1f}' + r'\, \nu_w $'
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

ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2.7,5.3])
ax.set_ylim([2e0,70])



#%% Fit effective viscosity using bilayer viscous model 





