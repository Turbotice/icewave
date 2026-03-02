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


#%% Set figure folder 

fig_folder = 'U:/Aurore_frasil/Results_Seb/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load experimental results

file2load = 'U:/Aurore_frasil/attenuation_results.pkl'
with open(file2load,'rb') as pf:
    data = pickle.load(pf)

#%% Load table of effective viscosity

file2effective_table = 'U:/Aurore_frasil/effective_nu_Lamb.h5'
if os.path.isfile(file2effective_table):
    print(f'{file2load} already exists, loading...')
    effective_table = rw.load_dict_from_h5(file2effective_table)
else:
    print(f'{file2load} does not exist')
    effective_table = {}

#%% Plot attenuation coeff for a single method and all thicknesses


set_graphs.set_matplotlib_param('single')
method = 'PIV'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots(figsize = (12,9))

effective_table[method] = {'h':[],'nu':[],'prefactor':[],
                           'err_nu':[],'err_prefactor':[]}

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
    
    effective_table[method]['h'].append(h)
    effective_table[method]['prefactor'].append(np.exp(coeff))
    effective_table[method]['err_prefactor'].append(np.exp(coeff)*err_coeff)
    
    nu = 2*pow(np.exp(coeff),2)*pow(g,4)
    err_nu = 4*pow(g,4)*np.exp(coeff)*np.exp(coeff)*err_coeff
    effective_table[method]['nu'].append(nu)
    effective_table[method]['err_nu'].append(err_nu)
    
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


figname = f'{fig_folder}alpha_VS_freq_surface_Lamb_attenuation_{method}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

for key in effective_table[method].keys():
    effective_table[method][key] = np.array(effective_table[method][key])

rw.save_dict_to_h5(effective_table,file2effective_table)

#%% Plot effective viscosity as a function of h

nu_w = 1e-6 # kinematik viscosity of water 
set_graphs.set_matplotlib_param('single')
color_table = {'laser':'tab:blue','PIV':'tab:orange'}
fig, ax = plt.subplots()
for method in effective_table.keys():
    
    ax.errorbar(effective_table[method]['h'],effective_table[method]['nu']/nu_w,
                yerr = effective_table[method]['err_nu'],fmt = 'o',color = color_table[method],
                label = method,markersize = 8)

ax.set_xlabel(r'$h \; \mathrm{(mm)}$')
ax.set_ylabel(r'$\nu_{eff} / \nu_w$')
ax.legend()

figname = f'{fig_folder}nu_eff_VS_h_Lamb_attenuation'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


# =============================================================================
#%% Fit effective viscosity using bilayer viscous model - imaginary part
# =============================================================================

#%% Load table of effective viscosity
file2effective_bilayer = 'U:/Aurore_frasil/effective_nu_viscous_bilayer_imag.h5'
if os.path.isfile(file2effective_bilayer):
    print(f'{file2load} already exists, loading...')
    effective_bilayer = rw.load_dict_from_h5(file2effective_bilayer)
else:
    print(f'{file2load} does not exist')
    effective_bilayer = {}


#%% Fit effective viscosity using bilayer model 

method = 'laser'
effective_bilayer[method] = {'h':[],'nu':[],'err_nu':[]}


# physical properties
f_th = np.linspace(1e0,1e1,500)
rho_1 = 1e3 # water density
r = 1
rho_2 = rho_1/r
nu_2 = 1e-6 # water viscosity
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0] # in mm

# properties for plot
zoomed_bool = False
bounds = np.array([3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256) 
marker_list = ['s','D','^','h','p'] # list of markers 
mksize = 10

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for i in range(1,len(thickness_list)):
    h_label = thickness_list[i]
    h_value = thickness_list[i]*1e-3
    
    # get experimental values
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h_label:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    # curve_fit
    popt,pcov = scipy.optimize.curve_fit(lambda freq,nu_1 : 
                                         np.imag(get_wavevector_theory(freq, h_value, rho_1, rho_2, nu_1, nu_2)),
                                         f_demod,alpha,p0 = 1e-3, bounds = (1e-5,1e-1))
                
    nu = popt[0]
    err_nu = np.sqrt(np.diag(pcov))[0]
    
    effective_bilayer[method]['h'].append(h_label)
    effective_bilayer[method]['nu'].append(nu)
    effective_bilayer[method]['err_nu'].append(err_nu)
    
    label_th = r'$\nu_{1} = ' + f'{nu*1e6:.1f}' + r'\, \nu_w $'
    
    # plot theory 
    current_color = new_blues(norm(h_label))
    current_marker = marker_list[i-1]
    k_th = get_wavevector_theory(f_th, h_value, rho_1, rho_2, nu, nu_2)
    ax.plot(f_th,np.imag(k_th),'-',color= current_color,label = label_th)
    
    # plot experimental results
    ax.plot(f_demod,alpha,current_marker,color = current_color,markeredgecolor = 'k',markersize = mksize)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\mathrm{Im}(k) \; \mathrm{(m^{-1})}$')
ax.legend()

# set colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')

if zoomed_bool:
    ax.set_xlim([2.7,5.3])
    # ax.set_ylim([2e0,70])
else:
    ax.set_xlim([2,10])
    ax.set_ylim([1e-2,2e2])

figname = f'{fig_folder}alpha_VS_freq_bilayer_viscous_viscous_{method}_zoomed_{zoomed_bool}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# convert list in dictionnary to array and save it
for key in effective_bilayer[method].keys():
    effective_bilayer[method][key] = np.array(effective_bilayer[method][key])
rw.save_dict_to_h5(effective_bilayer,file2effective_bilayer)


# =============================================================================
#%% Plot top viscosity nu_1 as a function of h
# =============================================================================

file2effective_bilayer = 'U:/Aurore_frasil/effective_nu_viscous_bilayer_imag.h5'
effective_bilayer = rw.load_dict_from_h5(file2effective_bilayer)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
nu_w = 1e-6

color_table = {'laser':'tab:blue','PIV':'tab:orange'}

for method in effective_bilayer.keys():
    
    ax.errorbar(effective_bilayer[method]['h'],effective_bilayer[method]['nu']/nu_w,
                yerr = effective_bilayer[method]['err_nu'],fmt = 'o',color = color_table[method],
                label = method, markersize = 8)

ax.set_xlabel(r'$h \; \mathrm{(mm)}$')
ax.set_ylabel(r'$\nu_{1} / \nu_w$')
ax.legend()

figname = f'{fig_folder}nu_1_VS_h_bilayer_viscous_viscous'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


# =============================================================================
#%% Plot viscosity using both Lamb model and bilayer viscous / viscous
# =============================================================================

file2effective_bilayer = 'U:/Aurore_frasil/effective_nu_viscous_bilayer_imag.h5'
effective_bilayer = rw.load_dict_from_h5(file2effective_bilayer)

file2effective_table = 'U:/Aurore_frasil/effective_nu_Lamb.h5'
effective_table = rw.load_dict_from_h5(file2effective_table)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
nu_w = 1e-6

color_table = {'laser':'tab:blue','PIV':'tab:orange'}
mksize = 8

for method in effective_bilayer.keys():
    
    # plot bilayer model
    label = f'bilayer - {method}'
    ax.errorbar(effective_bilayer[method]['h'],effective_bilayer[method]['nu']/nu_w,
                yerr = effective_bilayer[method]['err_nu'],fmt = 'o',color = color_table[method],
                label = label, markersize = mksize)
    
    # plot Lamb model
    label = f'Lamb - {method}'
    ax.errorbar(effective_table[method]['h'],effective_table[method]['nu']/nu_w,
                yerr = effective_table[method]['err_nu'],fmt = '^',color = color_table[method],
                label = label, markersize = 8)

ax.set_xlabel(r'$h \; \mathrm{(mm)}$')
ax.set_ylabel(r'$\nu / \nu_w$')
ax.legend()

figname = f'{fig_folder}nu_VS_h_model_comparison_bilayer_Lamb'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')



