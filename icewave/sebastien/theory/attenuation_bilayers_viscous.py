# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 12:02:25 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors 
import scipy

import os

import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.sebastien.theory.module_bilayers_viscous as theory


# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# Set global parameters
global g
g = 9.81

main_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/waves_attenuation/'
fig_folder = f'{main_path}attenuation/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
 
#%% FUNCTION SECTION 

def plot_det_M(x,y,params,figname):
    """ Compute determinant of matrix M for different values of k
    Inputs : - x, numpy array, real part of wavevector
             - y, numpy array, imaginary part of wavevector
             - params, tuple (f_ex,h,rho_1,rho_2,nu_1,nu_2)
             - figname, name of the figure to be saved """
    
    f_ex = params[0]
    omega = 2*np.pi*f_ex
    
    k_linear = omega**2/g
    
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    # Evaluate |det(M(k))| over the grid

    det_values = np.zeros_like(Z, dtype=complex)

    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            k = Z[i, j]

            try:
                det_values[i, j] = theory.det_M(k,params)
            except np.linalg.LinAlgError:
                det_values[i, j] = np.nan  # handle singularities or bad evals
                
    
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots()
    extents = np.array((x.min(),x.max(),y.min(),y.max()))/k_linear
    c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
                  interpolation = 'gaussian',extent = extents)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
    ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    ax.set_ylabel(r'$\mathrm{Im}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    
    ax.set_title('abs(det(M(k))) in complex plane')
    fig_label = f'det_M_{suffixe}'
    figname = f'{sub_fig_folder}{fig_label}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    print(f'figure saved under name: {figname}')
    
def get_wavevector(params,k0):
    """ Compute wavevector that is a root of det(M) 
    Inputs: - params, tuple (f_ex,h,rho_1,rho_2,nu_1,nu_2)
            - k0, complex, initial guess for which we look for a root using Newton-Raphson algorithm
    Output: - new_k, complex, wavector which is a root of det(M) after performing Newton-Raphson algorithm """

    tol = 1e-2 # tolerence below which we stop Newton Algorithm
    max_iter = 100 # maximum iteration 
    step_derivative = 1e0 # step in wavevector space to compute derivative
    threshold_derivative = 1e-100 # threshold below which derivative is considered to be too small 
    new_k = theory.newton_raphson(params, k0,tol,max_iter,step_derivative,threshold_derivative)
    print(f'Wavevector from newton method : {new_k}')
    
    return new_k

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

 
#%% Set parameters

h = 1e-2 # thickness of frasil in meters
rho_1 = 1000 # in kg/m3
rho_2 = 1000 # in kg/m3
nu_1 = 1e-4 # in m2/s
nu_2 = 1e-6 # in m2/s
f_ex = 4.0 # in Hz

r = rho_1/rho_2

omega = 2*np.pi*f_ex
k_linear = omega**2/g # wave vector for non dissipative linear waves 

xi_0 = 1.5e-3 # wave amplitude in meter

params = np.array([f_ex,h,rho_1,rho_2,nu_1,nu_2])

# Order of magnitudes
U_0 = xi_0*omega
p_0 = rho_1*g*xi_0
vort_0 = U_0/h

suffixe = f'h_{h:.3f}_r_{r:.1f}_nu1_{nu_1:.2e}_nu2_{nu_2:.2e}'
sub_fig_folder = f'{fig_folder}/{suffixe}/'
if not os.path.isdir(sub_fig_folder):
    os.mkdir(sub_fig_folder)
    
    
    
###################### COMPUTE ATTENUATION FOR SEVERAL NU_1 ########################
#%%
h = 1e-2 # thickness of frasil in meters
rho_1 = 1000 # in kg/m3
rho_2 = 1000 # in kg/m3
nu_2 = 1e-5 # in m2/s

nu_1_array = np.logspace(-6,-3,10)

f_ex_array = np.arange(2.0,6.1,0.1)
omega_array = 2*np.pi*f_ex_array

k_array = np.zeros((len(nu_1_array),len(f_ex_array)),dtype = 'complex')

x = np.linspace(1e-1,200,100)
y = np.linspace(1e-1,200,100)
for i,nu_1 in enumerate(nu_1_array):
    
    suffixe = f'h_{h:.3f}_r_{r:.1f}_nu1_{nu_1:.2e}_nu2_{nu_2:.2e}'
    sub_fig_folder = f'{fig_folder}/{suffixe}/'
    if not os.path.isdir(sub_fig_folder):
        os.mkdir(sub_fig_folder)
        
    for j,f_ex in enumerate(f_ex_array):
        params = (f_ex,h,rho_1,rho_2,nu_1,nu_2)
        
        if j == len(f_ex_array)//2:
            
            fig_label = f'det_M_{suffixe}'
            figname = f'{sub_fig_folder}{fig_label}'
            plot_det_M(x, y, params, figname)
            
        plt.close('all')
        
        omega = 2*np.pi*params[0]
        k_linear = omega**2/g # wave vector for non dissipative linear waves 
        k0 = k_linear*(1*0.9 + 1j*0.1)  # initial wavevector to start Newton algorithm 
        print(f'k = omega2/g leads to : {k_linear}')

        new_k = get_wavevector(params, k0)
        k_array[i,j] = new_k

#%% Create colorbar for values of nu_1

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))
norm = mcolors.LogNorm(nu_1_array.min(),nu_1_array.max()) 

#%% Plot dispersion relation 

k_gravity = omega_array**2/g
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
for i in range (len(nu_1_array)):
    nu_1 = nu_1_array[i]
    current_color = new_blues(norm(nu_1))
    ax.plot(np.real(k_array[i,:]),omega_array,'o',color = current_color)
    
ax.plot(k_gravity,omega_array,'r--')
ax.set_xlabel(r'$\mathrm{Re}(k) \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r'$\nu_1 \; \mathrm{(m^2.s^{-1})}$')
suffixe = f'h_{h:.1e}_r_{r:.2f}_nu_2_{nu_2:.2e}'
figname = f'{fig_folder}disp_relation_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot attenuation 

k_gravity = omega_array**2/g
set_graphs.set_matplotlib_param('single')
f_ex_th = np.linspace(1e0,1e1,100)

fig, ax = plt.subplots()
for i in range (len(nu_1_array)):
    nu_1 = nu_1_array[i]
    current_color = new_blues(norm(nu_1))
    ax.plot(f_ex_array,np.imag(k_array[i,:]),'o',color = current_color)
    
    pfit,err_pfit = powerlaw_fit(f_ex_array,np.imag(k_array[i,:]))
    beta = pfit[0]
    B = pfit[1]
    yth = B*f_ex_th**beta
    print(beta,B)
    label_fit = r'$q = ' + f'{B:.2f}' + 'f^{' + f'{beta:.1f}'  + '}$'
    ax.plot(f_ex_th,yth,color = current_color)
    
    
ax.set_ylabel(r'$\mathrm{Im}(k) \; \mathrm{(rad.m^{-1})}$')
ax.set_xlabel(r'$f_{ex} \; \mathrm{(Hz)}$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r'$\nu_1 \; \mathrm{(m^2.s^{-1})}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1e0,1e1])
ax.set_ylim([1e-4,1e3])

suffixe = f'h_{h:.1e}_r_{r:.2f}_nu_2_{nu_2:.2e}'
figname = f'{fig_folder}attenuation_VS_freq_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

###################################################################################
##################### COMPUTE ATTENUATION FOR A SINGLE VISCOSITY ##################
#%% Compute attenuation for various frequencies 

f_ex_array = np.arange(2.0,7.1,0.1)
omega_array = 2*np.pi*f_ex_array
params_array = [(f_ex,h,rho_1,rho_2,nu_1,nu_2) for f_ex in f_ex_array]

k = np.zeros(f_ex_array.shape,dtype = 'complex')
for i,params in enumerate(params_array):
    
    omega = 2*np.pi*params[0]
    k_linear = omega**2/g # wave vector for non dissipative linear waves 
    k0 = k_linear*(1 + 1j*0.1)  # initial wavevector to start Newton algorithm 
    print(f'k = omega2/g leads to : {k_linear}')

    new_k = get_wavevector(params, k0)
    k[i] = new_k


#%% Plot dispersion relation 

k_gravity = omega_array**2/g

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(np.real(k),omega_array,'o')
ax.plot(k_gravity,omega_array,'r--')
ax.set_xlabel(r'$\mathrm{Re}(k) \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')

#%% Plot attenuation 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(f_ex_array,np.imag(k),'o')

ax.set_xlabel(r'$f_{ex} \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\mathrm{Im}(k) \; \mathrm{(rad.m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')

#%% Compute determinant for different values of k 
x = np.linspace(1e-1,200,100)
y = np.linspace(1e-1,200,100)
fig_label = f'det_M_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'

plot_det_M(x, y, params, figname)


