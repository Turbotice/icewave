# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 11:48:27 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy
import math
import cv2 as cv

import h5py
import pickle
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.datafolders as df
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.drone.drone_projection as dp
import icewave.tools.interactive_scatter_filter as scatter_filter
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')


#%% FUNCTION SECTION 

def matrix_peak_correlation_detection(M,x,y,x_range,y_range,model_signal,detec_param):
    """ Perform peak detection using gaussian correlation 
    Inputs : - M, numpy.ndarray [ny,nx]: matrix over which we perform peak detection
             - x, array type (nx,) : detection is performed at fixed x
             - y, array type (ny,): correlation is performed with respect to y coordinate
             - x_range,y_range, array (2,) : range of x and y over which we look for maxima
             - model_signal, array type (ny,): reference array used to perform correlation
             - detec_param, dict : parameters dictionnary to be used for find_peaks """
    
    delta_y = np.diff(y)[0]
    
    mask_y = np.where(np.logical_and(y < y_range[1],y > y_range[0]))[0]
    y = y[mask_y]
    M = M[mask_y,:]

    idx_min = np.argmin(np.abs(x_range[0] - x))
    idx_max = np.argmin(np.abs(x_range[1] - x))
    indices = np.arange(idx_min,idx_max,step = 1)

    # parameters for find peaks
    min_prominence = detec_param['prominence']
    rel_height = detec_param['rel_height']
    
    S = {'idx':[],'peaks':[],'A':[],'y':[],'x':[]}
    # create dictionnary in which results are saved 
    
    fig, ax = plt.subplots()
    for i in range(len(indices)):
        ax.clear()
        
        idx = indices[i]
        
        detected_x = x[idx]
        print(detected_x)
        section = M[:,idx]
        label_x = f'x = {detected_x:.2f}'
        
        norm = section/section.max()
        
        ax.plot(y,norm,'-o',label = label_x)
        x_model = np.arange(0,len(model_signal),step = 1)*delta_y
        ax.plot(x_model,model_signal,'-k')
        
        # Apply gaussian filter 
        correlated = scipy.signal.fftconvolve(norm,model_signal,mode = 'same')
        norm_correlated = correlated/correlated.max()
        ax.plot(y,norm_correlated)
        # find peaks of correlation
        peaks,properties = scipy.signal.find_peaks(norm_correlated,prominence = min_prominence)
        prominences = scipy.signal.peak_prominences(norm_correlated,peaks)
        widths = scipy.signal.peak_widths(norm_correlated,peaks,rel_height = rel_height)
        ax.vlines(x=y[peaks], ymin=norm_correlated[peaks] - prominences[0],
                ymax = norm_correlated[peaks],color = 'g')
        ax.hlines(y=widths[1], xmin=y_range[0] + widths[2]*delta_y,
                xmax=y_range[0] + widths[3]*delta_y,color = 'g')
        
        ax.plot(y[peaks],norm_correlated[peaks],'d')
        
        ax.set_xlabel(r'$y \; \mathrm{(Hz)}$')
        ax.set_ylabel(r'$A_{norm}$')
        ax.legend()
        
        plt.pause(1.0)
        # save parameters in dictionnary
        for i,p in enumerate(peaks) :
            S['idx'].append(idx)
            S['peaks'].append(p)
            S['A'].append(norm[p]*section.max())
            S['x'].append(detected_x)
            S['y'].append(y[p])
        
        detected_y = y[p]
        print(f'y = {detected_y}')
        # f_txt = f'{detected_f:.2f}'.replace('.','p')
        # figname = f'{fig_folder}Section_f{f_txt}_corr_gaussian_sigma{gaussian_width}'
        # plt.savefig(figname + '.pdf', bbox_inches='tight')
        # plt.savefig(figname + '.png', bbox_inches='tight')
        
    for key in S.keys():
        S[key] = np.array(S[key])

    return S

#-----------------------------------------------------------------------------------------------------

def indices2fit(y,x,p,rel_height):
    """ Compute indices of points to be fitted, based on index of a peak p and the y-values to fit
    Inputs : y, numpy.ndarray, values to fit
             x, numpy.ndarray, abscisse coordinates
             p, index of the detected peak, y[p] is therefore a local maximum
             rel_height, float between 0 and 1, relative prominence at which peak width should be computed if left or right limits
             are too large
    Outputs : points2fit, indices of point to be fitted"""

    # compute left/right limits
    decreasing_curve = 1
    left_limit = p
    right_limit = p
    while decreasing_curve :
        next_left = left_limit-1
        next_right = right_limit+1
        if np.logical_or(y[next_left] >= y[left_limit],
                         y[next_right] >= y[right_limit]) : 
            decreasing_curve = 0
        else : 
            left_limit = next_left
            right_limit = next_right

        if np.logical_or(left_limit == 0,right_limit == len(x)-1):
            decreasing_curve = 0
            widths = scipy.signal.peak_widths(y,[p],rel_height = rel_height)
            left_limit = math.floor(widths[2][0])-1
            right_limit = math.ceil(widths[3][0])+1
            
    # compute left/right distance to peak
    left_ips = x[p] - x[left_limit]
    right_ips = x[right_limit] - x[p] 
    # min_ips = min(left_ips,right_ips)
    points2fit = np.where(np.logical_and(x > x[p] - left_ips,x < x[p] + right_ips))[0]
    # points2fit = np.where(np.logical_and(Afk['k'] > Afk['k'][p] - min_ips,Afk['k'] < Afk['k'][p] + min_ips))[0]
    
    return points2fit

#----------------------------------------------------------------------------------------------------------------

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

#----------------------------------------------------------------------------------------------------------------

def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN

#----------------------------------------------------------------------------------------------------------------

def phase_velocity_shallow(k,hw,g = 9.81):
    """ Compute phase velocity based on shallow water equation"""
    c_phi = np.sqrt(g*np.tanh(k*hw)/k)
    return c_phi

def group_velocity_shallow(k,hw,g = 9.81):
    """ Compute group velocity based on shallow water equation """
    c_phi = phase_velocity_shallow(k,hw)
    c_g =0.5*c_phi*(1 + 2*k*hw/np.sinh(2*k*hw))
    
    return c_g

#----------------------------------------------------------------------------------------------------------------

def time2space_attenuation(time_att,k,water_depth,error_bar = 1):
    """ Compute space attenuation from time attenuation using shallow water dispersion equation
    Inputs : - time_att, numpy.ndarray : N x 2 array, column 0 : time attenuation coefficient,
    column 1 : error on this coefficient
             - k, wavevector array
             - water_depth, tuple or array of size (2,), first value : water depth, second value : error on water depth
             - error_bar, boolean : default is 1, computes error_bar or not
    Outputs : - alpha, numpy.ndarray, spatial attenuation coefficients
              - err_alpha, numpy.ndarray, error on spatial attenuation coefficient
    """
    
    lambda_time = time_att[:,0]
    hw = water_depth[0]
    c_g = group_velocity_shallow(k,hw)
    alpha = lambda_time/c_g
    
    # compute errorbars
    if error_bar :
        print('Computing error bars, be sure that errorbars are second column of time_att and water_depth')
        err_lambda = time_att[:,1]
        err_hw = water_depth[1]
        err_khw = k*err_hw
        c_phi = phase_velocity_shallow(k, hw)
        error_c_phi = c_phi*err_khw/np.sinh(2*k*hw)
        error_c_g = np.sqrt((0.5*np.sinh(2*k*hw) + k*hw)**2 * error_c_phi**2 + 
                            c_phi**2 * (1 - 2*k*hw/np.tanh(2*k*hw))**2 * err_khw**2)/np.sinh(2*k*hw)
        err_alpha = np.sqrt(err_lambda**2 + lambda_time**2 * (error_c_g/c_g)**2)/c_g   
    
    return alpha,err_alpha

#-----------------------------------------------------------------------------------------------------------------------------
# Create a function for fitting peaks by Lorentzian

def lorentzian_fit(signal,x,bounds_fit):
    """ Perform Lorentzian fit over a signal 
    Inputs : - signal, array like 
             - x, array like 
             - bounds_fit, tuple, bounds_fit = ([x0_min,alpha_min],[x0_max,alpha_max]) boundaries of possible parameters
             values, see documentation of scipy.optimize.curve_fit 
    Output :  - results, dictionnary of results """
    
    # peform Lorentzian fitting 
    y_exp = (signal - signal.min())/(signal.max() - signal.min())
    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                         bounds = bounds_fit)
    
    err_coeff = np.sqrt(np.diag(pcov))
    print(f'x0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')
     
    # compute distance to fit
    xfit = np.linspace(x[0],x[-1],len(x)) #pulsation
    yth = lorentzian(xfit, popt[0], popt[1])
    dist2fit = np.sum((y_exp - yth)**2)/np.sum(y_exp**2)

    # create a small dictionnary of results
    results = {}
    results['d'] = dist2fit
    results['A'] = signal.max()
    results['x0'] = popt[0]
    results['alpha'] = popt[1]
    results['err_x0'] = err_coeff[0]
    results['err_alpha'] = err_coeff[1]
    
    return results

#---------------------------------------------------------------------------------------------------------------------

def fit_peaks_by_lorentzian(signal,x,peaks,signal_model,rel_height, figname):
    """ Fit each peaks detected from a signal by a lorentzian curve
    Inputs : - signal, 1D array like 
             - x, 1D array like 
             - peaks, array like, contains indices of signal array at which peaks have been detected 
             - signal_model, array like, signal model with which points to be fitted by Lorentzian are computed 
             - rel_height, float, relative height of peaks used to compute points to be fitted by Lorentzian 
    Outputs : - results, dictionnary of results """
    
    # find points to fit
    norm = signal/signal.max()
    correlated = scipy.signal.fftconvolve(norm,signal_model,mode = 'same')
    norm_correlated = correlated/correlated.max()
    
    # initialize a dictionnary of results 
    results = {'x0':[],'alpha':[],'err_x0':[],'err_alpha':[],'A':[],'d':[]}
    fig, ax = plt.subplots()
    
    ax.plot(x,signal)
    
    for p in peaks:
        points2fit = indices2fit(norm_correlated,x,p,rel_height) 
        if signal[p] > signal[points2fit].min():
            
            # peform Lorentzian fitting 
            x2fit = x[points2fit]
            signal2fit = signal[points2fit]
            bounds_fit = ([x2fit[0],1e-4],[x2fit[-1],1])
            results_fit = lorentzian_fit(signal2fit,x2fit,bounds_fit)

            # create a small dictionnary of results
            for key in results.keys():
                results[key].append(results_fit[key])
                
            # plot lorentzian fit 
            ax.plot(x2fit,signal2fit,'.', color = 'r')
            x_fit = np.linspace(x2fit[0],x2fit[-1],len(x2fit))
            yth = lorentzian(x_fit,results_fit['x0'],results_fit['alpha'])*(signal2fit.max() - signal2fit.min()) + signal2fit.min()
            ax.plot(x_fit,yth,'r-')
            
    # save plot
    ax.set_xlabel(r'$x$',labelpad = 5)
    ax.set_ylabel(r'$A$',labelpad = 5)
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
        
    plt.close('all')
    
    for key in results.keys():
        results[key] = np.array(results[key])
            
    return results

#------------------------------------------------------------------------------------------------------------------

def temporal_attenuation_from_FK_spectrum(Efk,peaks_properties,f_range,k_range,rel_height,sub_folder):
    """ Compute temporal attenuation of detected peaks from a FK spectrum. 
    Inputs : - Efk, dictionnary, contains follwoing keys :+ E, array like [nf,nk], space-time spectrum 
                                                          + f, array like (nf,), frequency array 
                                                          + k, array like (nk,), wavevector array 
             - peaks_properties, dictionnary, containing at least keys : 'f','k' and 'peaks'. Contains peaks indices 
             and their associate coordinates (f,k). 
             - rel_height, float, relative height used to compute width, see function indices2fit 
             - sub_folder, string, folder where figures of each fit will be saved 
    Output : - m, dictionnary, contains all results """
    
    
    # select area of FK where we want to compute attenuation
    mask_f = np.where(np.logical_and(Efk['f'] < frequency_range[1],Efk['f'] > frequency_range[0]))[0]
    mask_k = np.where(np.logical_and(Efk['k'] > wavevector_range[0],Efk['k'] < wavevector_range[1]))[0]

    f = Efk['f'][mask_f]
    k = Efk['k'][mask_k]
    M = Efk['E'][np.ix_(mask_f,mask_k)]
    
    m = {'A':[],'f':[],'k':[],'lambda':[],'err_f':[],'err_lambda':[],'d':[]}
        
    for idx, current_k in enumerate(k):
        
        # select peaks detected at current wavevector
        mask = np.where(peaks_properties['k'] == current_k)[0]
        peaks = peaks_properties['peaks'][mask]
        
        if len(peaks) > 0:
        
            k_txt = f'{current_k:.2f}'.replace('.','p')
            figname = f'{sub_folder}Lorentzian_fit_respectto_f_fixed_k_{k_txt}'
            
            # compute attenuation 
            current_results = fit_peaks_by_lorentzian(M[:,idx], f*2*np.pi, peaks, gaussian, rel_height, figname)
            m['A'].append(current_results['A'])
            m['f'].append(current_results['x0']/2/np.pi)
            m['err_f'].append(current_results['err_x0']/2/np.pi)
            m['k'].append(current_k)
            m['lambda'].append(current_results['alpha'])
            m['err_lambda'].append(current_results['err_alpha'])
            m['d'].append(current_results['d'])
        
    for key in m.keys():
        m[key] = np.array(m[key])
        m[key] = np.squeeze(m[key])
        
    return m

#----------------------------------------------------------------------------------------------------------
def fit_water_height(f,k,fun,err_f = None):
    """ Compute water height from dispersion relation fit 
    Inputs : - f, array like, frequencies
             - k, array like, wavevectors
             - fun, function of k and other parameters, hw must be the first parameter of this function 
             - err_f, optional, error on f 
    Outputs : - hw and err_hw, water height and its standard deviation computed from fit """
    popt,pcov = scipy.optimize.curve_fit(fun,k,f,sigma = err_f,absolute_sigma = True)
    
    hw = popt[0]
    err_hw = np.sqrt(np.diag(pcov)[0])  
    return hw,err_hw

#----------------------------------------------------------------------------------------------------------
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







#%% Import data

base = 'F:/Rimouski_2024/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '10-waves_005'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
filelist = glob.glob(f'{path2data}*scaled.mat')
print(filelist)

idx_file = 0
file2load = filelist[idx_file]

# load file 
with h5py.File(file2load, 'r') as fmat:
    S = {}
    print('Top-level keys : ', list(fmat.keys()))
    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    S = mat2py.transpose_PIVmat_fields(S)
    
#%%  Supress quadratic noise

Vx = FT.supress_quadratic_noise(np.transpose(S['Vx'],(1,0,2)),S['x'],S['y'])
Vy = FT.supress_quadratic_noise(np.transpose(S['Vy'],(1,0,2)),S['x'],S['y'])
Vx = np.transpose(Vx,(1,0,2))
Vy = np.transpose(Vy,(1,0,2))
# Vx = np.flip(Vx,(0,1))
# Vy = np.flip(Vy,(0,1))

print('Quadratic field supressed')

#%% Define fig_folder

fig_folder = f'{path2data}Figures_attenuation_ux_uz/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Show apparent velocity fields Vx and Vy

frame = 0
extents_meter = np.array([S['x'].min(),S['x'].max(),S['y'].min(),S['y'].max()])
extents_pix = np.array([S['PIXEL']['x_pix'].min(),S['PIXEL']['x_pix'].max(),
                    S['PIXEL']['y_pix'].min(),S['PIXEL']['y_pix'].max()])

set_graphs.set_matplotlib_param('single')
fig, axs = plt.subplots(ncols = 2,sharey = True, figsize = (12,8))
imsh = axs[0].imshow(Vx[:,:,frame].T,cmap = parula_map,origin = 'lower',aspect = 'equal',extent = extents_meter)
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_x \; \mathrm{(m.s^{-1})}$')

imsh = axs[1].imshow(Vy[:,:,frame].T,cmap = parula_map,origin = 'lower',aspect = 'equal',extent = extents_meter)
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_y \; \mathrm{(m.s^{-1})}$')

axs[0].set_ylabel(r'$y \; \mathrm{(m)}$')
for ax in axs:
    ax.set_xlabel(r'$x \; \mathrm{(m)}$')

plt.tight_layout()

#%% Compute ux and uz 

uz,ux,err_uz = dp.get_uz_ux_from_structure(Vx,Vy,S)

#%% Show spatio for uz and err_uz

extents_spatio = [S['t'][0],S['t'][-1],S['x'][0],S['x'][-1]]
fig,ax = plt.subplots()
imsh = ax.imshow(uz,origin = 'lower',aspect = 'auto',cmap = 'seismic',
          interpolation = 'gaussian',extent = extents_spatio)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_z \; \mathrm{(m.s^{-1})}$')

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

extents_spatio = [S['t'][0],S['t'][-1],S['x'][0],S['x'][-1]]
fig,ax = plt.subplots()
imsh = ax.imshow(err_uz,origin = 'lower',aspect = 'auto',cmap = parula_map,
          interpolation = 'gaussian',extent = extents_spatio)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_z \; \mathrm{(m.s^{-1})}$')

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)


#%% Show field for ux

frame = 800
fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,9))
imsh = axs[0].imshow(ux[:,:,frame].T,cmap = parula_map,origin = 'lower',extent = extents_meter)
divider = make_axes_locatable(axs[0])
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_x \; \mathrm{(m.s^{-1})}$')

imsh = axs[1].imshow(Vx[:,:,frame].T,cmap = parula_map,origin = 'lower',extent = extents_meter)
divider = make_axes_locatable(axs[1])
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_x \; \mathrm{(m.s^{-1})}$')

axs[0].set_ylabel(r'$y \; \mathrm{(m)}$')
axs[0].set_xlabel(r'$x \; \mathrm{(m)}$')
axs[1].set_xlabel(r'$x \; \mathrm{(m)}$')

plt.tight_layout()

#%% Compute space-time spectrum for ux

N = ux.shape[2]
Efk = FT.space_time_spectrum(ux,1/S['SCALE']['fx'],S['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Plot space time spectrum 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
Amin = 6e-6 # change to adjust colormap
Amax = 2e-3 # change to adjust colormap
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
kbounds = [0.05,3.3] # bounds for k axis on Efk plot
fbounds = [0.05,1.5] # bounds for f axis on Efk plot
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{u}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}ux_spacetime_spectrum_raw_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Detect peaks from dispersion relation curve

# Build Gaussian to be convoluted (can change gaussian width to detect more or less peaks)
gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.01,2.5] # range of wavevector spanned
frequency_range = [0,1.5] # range of frequency over which we look for peaks

Speaks = matrix_peak_correlation_detection(Efk['E'], Efk['k'], Efk['f'], wavevector_range, frequency_range, 
                                      gaussian, detec_param)
# change key name
Speaks['k'] = Speaks['x']
Speaks['f'] = Speaks['y']
del Speaks['x']
del Speaks['y']

#%% Filter detected points, plot results and save filtered peaks

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(Speaks['k'], Speaks['f'],Speaks)

fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_respectto_f_fixed_k_{date}_{drone_ID}_{exp_ID}.h5'
rw.save_dict_to_h5(filtered_properties, file2save)
print('DONE.')


#%% Compute temporal attenuation from Efk using Lorentzian fit 

rel_height = detec_param['rel_height']

sub_folder = f'{fig_folder}Lorentzian_fit/'
if not os.path.isdir(sub_folder):
    os.mkdir(sub_folder)

m = temporal_attenuation_from_FK_spectrum(Efk,filtered_properties,frequency_range,wavevector_range,rel_height,sub_folder)
    
#%% Compute water depth hw

fun = lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi
hw,err_hw = fit_water_height(m['f'], m['k'], fun, err_f = m['err_f'])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.errorbar(m['k'],m['f'],yerr = m['err_f'],fmt = '.',color = 'w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{u}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'

k_fit = np.linspace(kbounds[0],kbounds[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
ax.legend()

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}ux_spacetime_spectrum_fit_hw_{hw_txt}_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw

    
#%% Compute spatial attenuation coefficient 

time_att = np.column_stack((m['lambda'],m['err_lambda'])) # time attenuation coefficient and standard deviation
water_depth = [hw , err_hw] # water depth and standard deviation
alpha,err_alpha = time2space_attenuation(time_att, m['k'], water_depth) # compute spatial attenuation

xbounds = [0.1,1]
ybounds = [1e-3,1e0]
x_fit = np.linspace(xbounds[0],xbounds[1],100)

# fit by a powerlaw
coeffs,err_coeffs = powerlaw_fit(m['f'], alpha, err_alpha)
yth = coeffs[1]*x_fit**coeffs[0]

label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'

fig, ax = plt.subplots()
ax.errorbar(m['f'],alpha,yerr = err_alpha,xerr = m['err_f'],fmt = 'o')
ax.plot(x_fit,yth,linewidth = 2,label = label_th)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_from_temporal_attenuation_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  

# Save data 
m['alpha'] = alpha
m['err_alpha'] = err_alpha
m['power_law'] = {}
m['power_law']['B'] = coeffs[1]
m['power_law']['err_B'] = err_coeffs[1]
m['power_law']['beta'] = coeffs[0]
m['power_law']['err_beta'] = err_coeffs[0]

file2save = f'{fig_folder}attenuation_from_Efk_fixed_k_respect2f_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)
    
    
#%% Create main_results structure and save it 

main_results = {}
main_results['date'] = date
main_results['drone_ID'] = drone_ID
main_results['exp_ID'] = exp_ID
main_results['DRONE'] = S['DRONE']
main_results['SCALE'] = S['SCALE']
main_results['GPS'] = S['GPS']
main_results['t0_UTC'] = S['t0_UTC']

# main_results['Efk'] = Efk
main_results['attenuation'] = m

file2save = f'{fig_folder}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(main_results,pf)

print('Main results file saved !')






###############################################################################################################
################### Compute attenuation coefficient using fixed frequency #########################################################################
###############################################################################################################

#%% Detect harmonics for a given frequency and perform peak detection
# Build Gaussian to be convoluted
gaussian_width = 2
N = 6
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.4} # parameters for find peaks function
frequency_range = [0.1,1.1] 
wavevector_range = [0.05,3.0]

Speaks = matrix_peak_correlation_detection(Efk['E'].T, Efk['f'], Efk['k'], frequency_range, wavevector_range, 
                                           gaussian, detec_param)
# change key name
Speaks['f'] = Speaks['x']
Speaks['k'] = Speaks['y']
del S['x']
del S['y']


#%% Filter detected points, plot results and save filtered peaks

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(Speaks['k'], Speaks['f'],Speaks)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_respectto_k_fixed_f_{date}_{drone_ID}_{exp_ID}.h5'
rw.save_dict_to_h5(filtered_properties, file2save)
print('DONE.')
    
#%% Show detected peaks superposed with space-time spectrum

fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'.w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}Detected_fk_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute peaks width and attenuation coefficients 

sub_folder = f'{fig_folder}Lorentzian_fit_respectto_k/'
peaks_properties = filtered_properties

# select area of FK where we want to compute attenuation
mask_f = np.where(np.logical_and(Efk['f'] < frequency_range[1],Efk['f'] > frequency_range[0]))[0]
mask_k = np.where(np.logical_and(Efk['k'] > wavevector_range[0],Efk['k'] < wavevector_range[1]))[0]

f = Efk['f'][mask_f]
k = Efk['k'][mask_k]
M = Efk['E'][np.ix_(mask_f,mask_k)]

m = {'A':[],'f':[],'k':[],'alpha':[],'err_k':[],'err_alpha':[],'d':[]}
    
for idx, current_f in enumerate(f):
    
    # select peaks detected at current wavevector
    mask = np.where(peaks_properties['f'] == current_f)[0]
    peaks = peaks_properties['peaks'][mask]
    
    if len(peaks) > 0:
    
        f_txt = f'{current_f:.2f}'.replace('.','p')
        figname = f'{sub_folder}Lorentzian_fit_respectto_k_fixed_f_{f_txt}'
        
        # compute attenuation 
        current_results = fit_peaks_by_lorentzian(M.T[:,idx], k, peaks, gaussian, rel_height, figname)
        m['A'].append(current_results['A'])
        m['k'].append(current_results['x0'])
        m['err_k'].append(current_results['err_x0'])
        m['f'].append(current_f)
        m['alpha'].append(current_results['alpha'])
        m['err_alpha'].append(current_results['err_alpha'])
        m['d'].append(current_results['d'])
    
for key in m.keys():
    m[key] = np.array(m[key])
    m[key] = np.squeeze(m[key])
    

#%% Compute water depth hw

fun = lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi
hw,err_hw = fit_water_height(m['f'], m['k'], fun)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.errorbar(m['k'],m['f'],xerr = m['err_k'],fmt = '.',color = 'w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{u}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'

k_fit = np.linspace(kbounds[0],kbounds[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
ax.legend()

hw_txt = f'{hw:.2f}'.replace('.','p')
# figname = f'{fig_folder}ux_spacetime_spectrum_fit_hw_{hw_txt}_{date}_{drone_ID}_{exp_ID}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw



#%%
xbounds = [0.1,1]
ybounds = [1e-3,1e0]
x_fit = np.linspace(xbounds[0],xbounds[1],100)

# fit by a powerlaw
coeffs,err_coeffs = powerlaw_fit(m['f'], m['alpha'], m['err_alpha'])
yth = coeffs[1]*x_fit**coeffs[0]

label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'

fig, ax = plt.subplots()
ax.errorbar(m['f'],m['alpha'],yerr = m['err_alpha'],fmt = 'o')
ax.plot(x_fit,yth,linewidth = 2,label = label_th)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_fixed_f_respectto_k_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  


m['power_law'] = {}
m['power_law']['B'] = coeffs[1]
m['power_law']['err_B'] = err_coeffs[1]
m['power_law']['beta'] = coeffs[0]
m['power_law']['err_beta'] = err_coeffs[0]

file2save = f'{fig_folder}attenuation_from_Efk_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)

#%% Compare both detections with respect to f / k 

# load both set of coordinates (f,k)
m = {}
file2save = f'{fig_folder}attenuation_from_Efk_fixed_k_respect2f_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'rb') as pf:
    m['time'] = pickle.load(pf)

file2save = f'{fig_folder}attenuation_from_Efk_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'rb') as pf:
    m['space'] = pickle.load(pf)

#%% Compare dispersion relations 

fig, ax = plt.subplots()
for key in m.keys():
    ax.plot(m[key]['k'],m[key]['f'],'.',label = key)
    
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.legend()

figname = f'{fig_folder}comparison_detection_methods_dispersion_relation_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compare attenuation laws
xbounds = np.array([1e-1,1e0])
ybounds = np.array([1e-3,1e0])
xfit = np.linspace(xbounds[0],xbounds[-1],100)

color_tab = {'time':'tab:blue','space':'tab:orange'}

fig, ax = plt.subplots()
for key in m.keys():
    ax.errorbar(m[key]['f'],m[key]['alpha'],yerr = m[key]['err_alpha'],fmt = '.',
                color = color_tab[key],label = key)
    # yth = m[key]['power_law']['B']*xfit**m[key]['power_law']['beta']
    # ax.plot(xfit,yth,'-',color = color_tab[key])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}comparison_detection_methods_attenuation_law_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')




