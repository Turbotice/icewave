# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 13:53:30 2025

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

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# Set global parameters
global g
g = 9.81

main_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/waves_attenuation/'
fig_folder = f'{main_path}attenuation_dimensionless/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
   
    
#%% FUNCTION SECTION 

# dimensionless numbers
def dim_numbers(params,g = 9.81):
    """ Return dimensionless numbers from set of parameters
    Inputs: - params, tuple, (f_ex,h,rho_1,rho_2,nu_1,nu_2)
    Output: - dim_numbers, tuple, (gamma,delta_1,delta_2,r) """
    
    omega = 2*np.pi*params[0]
    k_0 = omega**2/g
    
    gamma = params[1]*k_0
    delta_1 = np.sqrt(nu_1/omega)*k_0
    delta_2 = np.sqrt(nu_2/omega)*k_0
    r = rho_1/rho_2
    dimensionless = (gamma,delta_1,delta_2,r)
    return dimensionless


def plot_det_M(x,y,dimensionless):
    """ Compute determinant of matrix M for different values of k
    Inputs : - x, numpy array, real part of wavevector
             - y, numpy array, imaginary part of wavevector
             - dimensionless, tuple (gamma,delta_1,delta_2,r)
             - figname, name of the figure to be saved """
        
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    # Evaluate |det(M(k))| over the grid

    det_values = np.zeros_like(Z, dtype=complex)

    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            k = Z[i, j]

            try:
                det_values[i, j] = theory.det_M(k,dimensionless)
            except np.linalg.LinAlgError:
                det_values[i, j] = np.nan  # handle singularities or bad evals
                
    fig, ax = plt.subplots()
    extents = np.array((x.min(),x.max(),y.min(),y.max()))
    c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
                  interpolation = 'gaussian',extent = extents)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
    ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    ax.set_ylabel(r'$\mathrm{Im}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    
    return fig, ax 
    

def print_dimensionless(dimensionless):
    print(f'gamma = {dimensionless[0]:.2e}')
    print(f'delta_1 = {dimensionless[1]:.2e}')
    print(f'delta_2 = {dimensionless[2]:.2e}')
    print(f'r = {dimensionless[3]:.2e}')
    
def get_wavevector(dimensionless,k0):
    """ Compute wavevector that is a root of det(M) 
    Inputs: - dimensionless, tuple (gamma,delta_1,delta_2,r)
            - k0, complex, initial guess for which we look for a root using Newton-Raphson algorithm
    Output: - new_k, complex, wavector which is a root of det(M) after performing Newton-Raphson algorithm """

    tol = 1e-4 # tolerence below which we stop Newton Algorithm
    max_iter = 200 # maximum iteration 
    step_derivative = 1e-4 # step in wavevector space to compute derivative
    threshold_derivative = 1e-80 # threshold below which derivative is considered to be too small 
    new_k = theory.newton_raphson(dimensionless, k0,tol,max_iter,step_derivative,threshold_derivative)
    
    return new_k

def func_real_imag(z,dimensionless):
    x,y = z # z = [Re(z),Im(z)]
    
    kappa = x + 1j*y
    determinant = theory.det_M(kappa,dimensionless)
    return [np.real(determinant),np.imag(determinant)] # system of 2 real equations 


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
        vort_term = 1j*m_2*A_2*np.exp(m_2*kappa*z)
        
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

#%% Define input parameters 

h = 1e-2 # thickness of frasil in meters
rho_1 = 1000 # in kg/m3
rho_2 = 1000 # in kg/m3
nu_1 = 1e-6 # in m2/s
nu_2 = 1e-6 # in m2/s
f_ex = 3.0 # in Hz

xi_0 = 1.5e-3 # wave amplitude in meter

r = rho_1/rho_2

omega = 2*np.pi*f_ex
k_0 = omega**2/g # wave vector for non dissipative linear waves 

params = np.array([f_ex,h,rho_1,rho_2,nu_1,nu_2])

# Order of magnitudes
U_0 = xi_0*omega
p_0 = rho_1*g*xi_0
vort_0 = U_0*k_0

dimensionless = dim_numbers(params)
print_dimensionless(dimensionless)

suffixe = f'h_{h:.3f}_r_{r:.1f}_nu1_{nu_1:.2e}_nu2_{nu_2:.2e}'
sub_fig_folder = f'{fig_folder}/{suffixe}/'
if not os.path.isdir(sub_fig_folder):
    os.mkdir(sub_fig_folder)
    
#%% Plot determinant in complex space

x = np.linspace(1e-3,2,100)
y = np.linspace(-0.5,2,100)

figname = f'{sub_fig_folder}det_M_f_ex_{f_ex:.2f}_Hz'
plot_det_M(x, y, dimensionless, figname)

#%% Compute wavevector for a given set of dimensionless numbers

# set dimensionless numbers
gamma = 1.0
delta_1 = 5e-3
delta_2 = delta_1*1e-1
r = 1.0
dimensionless = (gamma,delta_1,delta_2,r)

x = np.linspace(1e-3,2,100)
y = np.linspace(1e-3,2,100)

suffixe = f'gamma_{gamma:.2e}_delta1_{delta_1:.2e}_delta2_{delta_2:.2e}_r_{r:.2f}'
sub_fig_folder = f'{fig_folder}/dimensionless_numbers/'
if not os.path.isdir(sub_fig_folder):
    os.mkdir(sub_fig_folder)

figname = f'{sub_fig_folder}suffixe'
set_graphs.set_matplotlib_param('single')
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# Evaluate |det(M(k))| over the grid

det_values = np.zeros_like(Z, dtype=complex)

for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        k = Z[i, j]

        try:
            det_values[i, j] = theory.det_M(k,dimensionless)
        except np.linalg.LinAlgError:
            det_values[i, j] = np.nan  # handle singularities or bad evals
            
fig, ax = plt.subplots()
extents = np.array((x.min(),x.max(),y.min(),y.max()))
c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
              interpolation = 'gaussian',extent = extents)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\mathrm{Im}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)

# compute wavevector for a set of dimensionless numbers 
initial_guess = 1 + 1j*0.1
kappa = get_wavevector(dimensionless, initial_guess)
ax.plot(np.real(initial_guess),np.imag(initial_guess),'ro')
ax.plot(np.real(kappa),np.imag(kappa),'rd')
print(f'Root of det(M) = {kappa}')

# Try to compute wavevector using fsolve function 

initial_guess = [1 , 0.1] # initial guess for [Re(z), Im(z)]
solution = scipy.optimize.fsolve(func_real_imag,initial_guess,(dimensionless,))
root = solution[0] + 1j*solution[1]

print(f'Root of det(M) = {root}')

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

######################################################################
#%% Compute attenuation for different dimensionless parameters
######################################################################

gamma_array = np.logspace(-2,0,10)
delta_1_array = np.logspace(-2,1,100)
ratio_delta = 1e-3
r = 1
# compute kappa for different delta_1 and different values of delta_1/gamma, keeping ratio delta_1/delta_2 = 1

# folder_study = f'{fig_folder}/lab_influence_delta_1_gamma_delta_ratio_{ratio_delta:.1e}'
# if not os.path.isdir(folder_study):
#     os.mkdir(folder_study)
        
x = np.linspace(1e-3,2,100)
y = np.linspace(-0.5,2,100)

kappa_array = np.zeros((len(gamma_array),len(delta_1_array)),dtype = 'complex')


for i,gamma in enumerate(gamma_array):
    # sub_fig_folder = f'{folder_study}/gamma_{gamma:.2e}/'
    # if not os.path.isdir(sub_fig_folder):
    #     os.mkdir(sub_fig_folder)
    
    initial_guess = [1 , 0.1] # initial guess for [Re(kappa), Im(kappa)]
    previous_root = initial_guess
    
    for k, delta_1 in enumerate(delta_1_array):
        delta_2 = delta_1*ratio_delta
        dimensionless = (gamma,delta_1,delta_2,r)
        
        suffixe = f'gamma_{gamma:.2e}_delta1_{delta_1:.2e}_delta2_{delta_2:.2e}_r_{r:.2f}'
        
        figname = f'{sub_fig_folder}{suffixe}'
        # set_graphs.set_matplotlib_param('single')
        # fig, ax = plot_det_M(x, y, dimensionless)
    
        # compute wavevector for a set of dimensionless numbers 
        solution = scipy.optimize.fsolve(func_real_imag,initial_guess,(dimensionless,))
        root = solution[0] + 1j*solution[1]
        print(f'Root of det(M) = {root}')
        initial_guess = solution
        kappa_array[i,k] = root
        
        # initial_guess = 0.9 + 1j*0.1
        # kappa = get_wavevector(dimensionless, initial_guess)
        # print(f'Root of det(M) = {kappa}')

        
        # ax.plot(np.real(initial_guess),np.imag(initial_guess),'ro')
        # ax.plot(np.real(kappa),np.imag(kappa),'rd')

        
        # plt.savefig(figname + '.pdf', bbox_inches='tight')
        # plt.savefig(figname + '.png', bbox_inches='tight')
        
        # plt.close('all')


#%% Plot dispersion relation 

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))
norm = mcolors.LogNorm(gamma_array.min(),gamma_array.max()) 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
for i in range(kappa_array.shape[0]):
    
    current_kappa = kappa_array[i,:]
    current_color = new_blues(norm(gamma_array[i]))
    ax.plot(delta_1_array,np.real(current_kappa),'-',color = current_color)
    
ax.set_xscale('log')
ax.set_xlabel(r'$\delta_1^*$')
ax.set_ylabel(r'$\mathrm{Re}(\kappa)$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r'$\gamma $')

label = r'$\delta_2^* / \delta_1^* = ' + f'{ratio_delta:.2e}' + r'$'
ax.set_title(label)

ax.set_xlim([7e-3,1.5e1])
ax.set_ylim([-0.1,1.5])

# results_folder = f'{fig_folder}lab_influence_delta_1_gamma/Results/'
# if not os.path.isdir(results_folder):
#     os.mkdir(results_folder)
# figname = f'{results_folder}dispersion_relation_delta_ratio_{ratio_delta:.2e}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot attenuation coefficient 

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))
norm = mcolors.LogNorm(gamma_array.min(),gamma_array.max()) 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
for i in range(kappa_array.shape[0]):
    
    current_kappa = kappa_array[i,:]
    current_color = new_blues(norm(gamma_array[i]))
    ax.plot(delta_1_array,np.imag(current_kappa),'-',color = current_color)
    
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\delta_1^*$')
ax.set_ylabel(r'$\mathrm{Im}(\kappa)$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=new_blues, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r'$\gamma $')
label = r'$\delta_2^* / \delta_1^* = ' + f'{ratio_delta:.2e}' + r'$'
ax.set_title(label)

ax.set_xlim([7e-4,1.5e1])
ax.set_ylim([6e-7,1e0])

# coeffs,err_coeffs = powerlaw_fit(delta_1_array[:-1],np.imag(kappa_array)[:-1])
# print(f'beta = {coeffs[0]:.2f} ± {err_coeffs[0]:.2f}')

# figname = f'{results_folder}Im_k_delta_ratio_{ratio_delta:.2e}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')





