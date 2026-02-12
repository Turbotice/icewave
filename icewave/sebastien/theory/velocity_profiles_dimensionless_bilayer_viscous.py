# -*- coding: utf-8 -*-
"""
Created on Fri Jul  4 10:20:52 2025

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
import icewave.sebastien.theory.module_bilayer_viscous_dimensionless as theory

import pickle

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# Set global parameters
global g
g = 9.81

main_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/waves_attenuation/'
fig_folder = f'{main_path}velocity_dimensionless/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% FUNCTION SECTION

def u_dimensionless(kappa,coeffs,dimensionless,z,fluid_idx):
    
    [A_1,B_1,C_1,D_1,A_2,C_2] = coeffs
    gamma, delta_1, delta_2, r = dimensionless
    
    m_1 = np.sqrt(kappa**2 - 1j*1/delta_1**2)
    m_2 = np.sqrt(kappa**2 - 1j*1/delta_2**2)
    
    if fluid_idx == 1:
        pressure_term = kappa*(C_1*np.exp(kappa*z) + D_1*np.exp(-kappa*z))
        vort_term = 1j*m_1*(A_1*np.exp(m_1*z) - B_1*np.exp(-m_1*z))
        
        u = pressure_term + vort_term
    
    else:
        pressure_term = r*kappa*C_2*np.exp(kappa*z) 
        vort_term = 1j*m_2*A_2*np.exp(m_2*z)
        
        u = pressure_term + vort_term
        
    return u
    
def w_dimensionless(kappa,coeffs,dimensionless,z,fluid_idx):
    
    [A_1,B_1,C_1,D_1,A_2,C_2] = coeffs
    gamma, delta_1, delta_2, r = dimensionless
    
    m_1 = np.sqrt(kappa**2 - 1j*1/delta_1**2)
    m_2 = np.sqrt(kappa**2 - 1j*1/delta_2**2)
    
    if fluid_idx == 1:
        pressure_term = -1j*kappa*(C_1*np.exp(kappa*z) - D_1*np.exp(-kappa*z))
        vort_term = kappa*(A_1*np.exp(m_1*z) + B_1*np.exp(-m_1*z))

        u = pressure_term + vort_term
        
    else:
        
        pressure_term = -1j*kappa*r*C_2*np.exp(kappa*z) 
        vort_term = kappa*A_2*np.exp(m_2*z)
        
        u = pressure_term + vort_term
        
    return u 

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























###################################################################################################################
#%% Plot velocity fields for a set of (nu_1, h)
###################################################################################################################

# load experimental data 
path2data =  'C:/Users/sebas/OneDrive/Bureau/These PMMH/attenuation_results.pkl'
with open(path2data,'rb') as pf:
    data = pickle.load(pf)
    
results_folder = f'{fig_folder}Experiments_regime/'
if not os.path.isdir(results_folder):
    os.mkdir(results_folder)    

#%% Compute wavevectors for different values of thickness and a given nu_1

f_array = np.linspace(1e0,1e1,500) # range of frequencies
omega_array = 2*np.pi*f_array
k0_array = omega_array**2/g
nu_2 = 1e-6 # kinematik viscosity of water

rho_1 = 1e3
r = 1
rho_2 = rho_1/r

nu_1_array = np.logspace(np.log10(nu_2),0,10)
# h_array = np.logspace(-3,np.log10(0.1),10)
h_array = np.array([5.0,7.5,10.0,12.5,15.0])*1e-3

# select a value of nu_1
# idx_nu1 = 2
# nu_1 = nu_1_array[idx_nu1]
nu_1 = 9e-3
print(f'nu_1 = {nu_1:.2e}')

kappa_array = np.zeros((len(h_array),len(f_array)),dtype = 'complex')
k_array = np.zeros((len(h_array),len(f_array)),dtype = 'complex')

for idx_h,h in enumerate(h_array):

    params_array = np.array([(f_ex,h,rho_1,rho_2,nu_1,nu_2) for f_ex in f_array])
    dimensionless_array = np.array([dim_numbers(params) for params in params_array])
    
    initial_guess = [1 , 0.1] # initial guess for [Re(kappa), Im(kappa)]
    for i,dimensionless in enumerate(dimensionless_array):
        
        solution = scipy.optimize.fsolve(func_real_imag,initial_guess,(dimensionless,))
        root = solution[0] + 1j*solution[1]
        print(f'Root of det(M) = {root}')
        initial_guess = solution
        kappa_array[idx_h,i] = root
        k_array[idx_h,i] = root*k0_array[i]



#%% Compute velocity field associated to a set of dimensionless numbers

coeffs = theory.get_coeffs(kappa, dimensionless)

#%% Compute velocity profile 

Nstep = 100 # spatial step in fluid 1
dz = -gamma/Nstep
H = 5*gamma # total depth 
z1 = np.arange(0,-gamma,dz)
z2 = np.arange(-gamma,-H,dz)
z = np.concatenate((z1,z2))

u1 = u_dimensionless(kappa, coeffs, dimensionless, z1, 1)
u2 = u_dimensionless(kappa, coeffs, dimensionless, z2, 2)
u = np.concatenate((u1,u2))

w1 = w_dimensionless(kappa, coeffs, dimensionless, z1, 1)
w2 = w_dimensionless(kappa, coeffs, dimensionless, z2, 2)
w = np.concatenate((w1,w2))

fig, ax = plt.subplots()
ax.plot(z,np.real(u),label = r'$\mathrm{Re}(u^*)$')
ax.plot(z,np.imag(u),label = r'$\mathrm{Im}(u^*)$')
ax.plot(z,np.real(w),label = r'$\mathrm{Re}(w^*)$')
ax.plot(z,np.imag(w),label = r'$\mathrm{Im}(w^*)$')


ax.set_xlabel(r'$z^*$',labelpad = 5)
ax.set_ylabel(r'$u^*$',labelpad = 5)
ax.legend()






