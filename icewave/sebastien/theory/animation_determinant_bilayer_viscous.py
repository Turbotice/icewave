# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 15:52:24 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

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
fig_folder = f'{main_path}animation/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% FUNCTION SECTION 

def dim_numbers(params,g = 9.81):
    """ Return dimensionless numbers from set of parameters
    Inputs: - params, tuple, (f_ex,h,rho_1,rho_2,nu_1,nu_2)
    Output: - dim_numbers, tuple, (gamma,delta_1,delta_2,r) """
    
    omega = 2*np.pi*params[0]
    k_0 = omega**2/g
    
    rho_1 = params[2]
    rho_2 = params[3]
    
    nu_1 = params[4]
    nu_2 = params[5]
    
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

def get_det_values(x,y,dimensionless):
    
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
    return det_values

def get_title(dimensionless):
    
    label = r'$\gamma = ' + f'{dimensionless[0]:.2e}'  + r'\; \delta_1 = ' 
    label = label + f'{dimensionless[1]:.2e}' + r'\; \delta_2 = ' + f'{dimensionless[2]:.2e}'  
    # label = label + r'\; r = ' + f'{dimensionless[3]:.2f}' + r'$'
    
    return label


def animation_det(dimensionless_array,x,y,nb_frames,time_interval):
    
    """ Create and return animation of a line plot using matplotlib.animation
    Inputs : - fig: matplotlib figure
         - ax: axis object
         - data: numpy array, data to plot #dim0 : space, #dim1 : time
         - x: numpy array, x-axis (space)
         - t: numpy array, (time)
         - nb_frames : int, number of frames to show
         - time_interval : time between two consecutive frames
         
    Output : - ani: matplotlib animation object"""
    
    extents = np.array((x.min(),x.max(),y.min(),y.max()))
    # Evaluate |det(M(k))| over the grid
    j0 = 0
    dimensionless = dimensionless_array[j0]
    det_values = get_det_values(x, y, dimensionless)
                
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots()
    c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
                  interpolation = 'gaussian',extent = extents)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
    ax.set_xlabel(r'$\mathrm{Re}(\kappa)$', labelpad = 5)
    ax.set_ylabel(r'$\mathrm{Im}(\kappa)$', labelpad = 5)
    
    label = get_title(dimensionless)
    ax.set_title(label)
    
    # small function to update the current figure 
    def update_plot(frame):
        dimensionless = dimensionless_array[frame]
        det_values = get_det_values(x, y, dimensionless)
        log_values = np.log10(abs(det_values))
        c.set_data(log_values)
        c.set_clim(np.min(log_values),np.max(log_values))
        
        label = get_title(dimensionless)
        ax.set_title(label)
        return c
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update_plot, frames=nb_frames, interval=time_interval)
    plt.show()
    print('Animation computed')
    
    return ani


#%%
gamma_array = np.logspace(-2,0,10)
for gamma in gamma_array:

    delta_1_array = np.logspace(-2,1,50)
    r = 1.0
    ratio_delta = 1e-1
    
    dimensionless_array = [(gamma,delta_1,delta_1*ratio_delta,r) for delta_1 in delta_1_array]
    
    x = np.linspace(1e-3,2,100)
    y = np.linspace(-0.5,2,100)
    
    # set parameters for animation
    nb_frames = len(delta_1_array)
    time_interval = 200
    
    ani = animation_det(dimensionless_array,x,y,nb_frames,time_interval)
    animation_title = f'det_M_gamma_{gamma:.2e}_ratio_delta_{ratio_delta:.2e}'
    file2save = f'{fig_folder}{animation_title}.mp4'
    ani.save(file2save)


#%%

label = get_title(dimensionless_array[0])


