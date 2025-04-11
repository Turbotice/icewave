# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 13:36:48 2025

@author: sebas
"""

import numpy as np
import os
import glob 
import pandas as pd 
import h5py
import pickle
import cv2 as cv
import scipy
import math
# from skimage.filters import meijering, sato, frangi, hessian


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
from matplotlib.path import Path

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 
def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN


def unnormed_gaussian(x,x0,alpha):
    y = np.exp(-0.5*((x - x0)/alpha)**2)
    return y

def lorentzian(x,x0,alpha):
    y = 1/(1 + ((x - x0)/alpha)**2)
    return y

def phase_velocity_shallow(k,hw,g = 9.81):
    """ Compute phase velocity based on shallow water equation"""
    c_phi = np.sqrt(g*np.tanh(k*hw)/k)
    return c_phi

def group_velocity_shallow(k,hw,g = 9.81):
    """ Compute group velocity based on shallow water equation """
    c_phi = phase_velocity_shallow(k,hw)
    c_g =0.5*c_phi*(1 + 2*k*hw/np.sinh(2*k*hw))
    
    return c_g

def time2space_attenuation(time_att,k,water_depth,error_bar = 1):
    """ Compute space attenuation from time attenuation using shallow water dispersion equation """
    
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
        error_c_g = np.sqrt((0.5*np.sinh(2*k*hw) + 2*k*hw)**2 * error_c_phi**2 + 
                            c_phi**2 * (1 - 2*k*hw/np.tanh(2*k*hw))**2 * err_khw**2)/np.sinh(2*k*hw)
        err_alpha = np.sqrt(err_lambda**2 + lambda_time**2 * (error_c_g/c_g)**2)/c_g   
    
    return alpha,err_alpha

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
        ax.hlines(y=widths[1], xmin=widths[2]*delta_y,
                xmax=widths[3]*delta_y,color = 'g')
        
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

#%% 
date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

path2data = f'K:/Share_hublot/Data/{date}/Drones/{drone_ID}/matData/{exp_ID}/'

file2load = glob.glob(f'{path2data}*scaled.mat')[0]
with h5py.File(file2load,'r') as fmat:
    print('Top-level keys : ', list(fmat.keys()))
    data = mat2py.mat_to_dict(fmat['m'],fmat['m'])

data = mat2py.transpose_PIVmat_fields(data)

#%% Define fig_folder and save folder

fig_folder = f'{path2data}Plots/attenuation_from_disp_relation/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

save_folder =f'{path2data}Results/'
if not os.path.isdir(save_folder):
    os.mkdir(save_folder)
    
#%% Reduce drone noise and compute FFT spectrum 
Vs = FT.supress_quadratic_noise(np.transpose(data['Vx'],(1,0,2)),data['x'],data['y'])

#%% Compute FFT spectrum 
TF_spectrum,freq = FT.temporal_FFT(Vs,data['SCALE']['facq_t'],padding_bool = 1,add_pow2 = 1)

#%%
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle _{x,y}(f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylim([1e-3,1e1])

figname = f'{fig_folder}TF_spectrum_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# save data in a structure
FFT_spectrum = {}
FFT_spectrum['TF_spectrum'] = TF_spectrum
FFT_spectrum['f'] = freq

#%% Compute 
Vs = np.transpose(Vs,(1,0,2))
Efk = FT.space_time_spectrum(Vs,1/data['SCALE']['fx'],data['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Plot Afk 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0.05,3.3)
ax.set_ylim(0.05,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}Efk_raw_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Method for peaks detection
# step 1 - find peaks using gaussian convolution
# step 2 - filter and keep only relevant peaks
# step 3 - for each relevant peak, compute points on which gaussian will be fitted
# step 4 - save alpha coefficient for each detected peak 

#%% Select a profile for a given wavevector and perform peak detection with gaussian fit over frequency 
Afk = Efk.copy()
# Build Gaussian to be convoluted
gaussian_width = 10
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6}
wavevector_range = [0.05,2.0]
frequency_range = [0,1.2]

S = matrix_peak_correlation_detection(Afk['E'], Afk['k'], Afk['f'], wavevector_range, frequency_range, gaussian, detec_param)
# change key name
S['k'] = S['x']
S['f'] = S['y']
del S['x']
del S['y']

#%% Show detected peaks superposed with space-time spectrum

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(S['k'], S['f'],S)

fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(0.09,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_respectto_f_fixed_k.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(filtered_properties,pf)
    
#%% Select a peak and cascade down to compute width over which we perform fitting

rel_height = detec_param['rel_height']

m = {'A':[],'f':[],'k':[],'lambda':[],'err_f':[],'err_lambda':[],'d':[]}

indices = np.where(np.logical_and(Afk['k'] > wavevector_range[0],Afk['k'] < wavevector_range[1]))[0]
selected_k = Afk['k'][indices]
fig, ax = plt.subplots()
for idx in indices:

    ax.clear()
    
    detected_k = Afk['k'][idx]
    print(f'k  = {detected_k:.2f}')
    label_k = f'{detected_k:.2f}'
    
    # select peaks detected at this frequency
    mask = np.where(filtered_properties['idx'] == idx)[0]
    peaks = filtered_properties['peaks'][mask]
    
    section = Afk['E'][:,idx]
    norm = section/section.max()
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    
    ax.plot(Afk['f'],norm,'-o',label = label_k)
    ax.plot(Afk['f'],norm_correlated,'-o')
    
    for p in peaks:

        points2fit = indices2fit(norm_correlated,Afk['f'],p,rel_height)
        # plot points
        ax.plot(Afk['f'][p],norm_correlated[p],'s')
        ax.plot(Afk['f'][points2fit],norm[points2fit],'d')  
        if section[p] > section[points2fit].min():
            
            x = Afk['f'][points2fit]
            y_exp = (section[points2fit] - section[points2fit].min())/(section[points2fit].max() - section[points2fit].min())
            
            # gaussian fit 
            popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                                 bounds = ([x[0],1e-5],[x[-1],10]))
            xfit = np.linspace(x[0],x[-1],len(x))
            yth = lorentzian(xfit, popt[0], popt[1])
            yth = yth*(section[points2fit].max() - section[points2fit].min()) + section[points2fit].min()
            err_coeff = np.sqrt(np.diag(pcov))
            print(f'f0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')

            ax.plot(xfit,yth,'r-')  
            
            # compute distance to fit
            dist2fit = np.sum((y_exp - yth)**2)/np.sum(y_exp**2)
            m['d'].append(dist2fit)
            m['A'].append(Afk['E'][p,idx])
            m['f'].append(popt[0])
            m['k'].append(detected_k)
            m['lambda'].append(popt[1])
            m['err_f'].append(err_coeff[0])
            m['err_lambda'].append(err_coeff[1])
        
    plt.pause(1.0)
        
for key in m.keys():
    m[key] = np.array(m[key])    
    
#%% Compute water depth hw

popt,pcov = scipy.optimize.curve_fit(lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi,m['k'],m['f'])
y_exp = bound_harmonicN(m['k'], 1, popt[0])/2/np.pi
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(m['k'],m['f'],'.w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(0.09,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

hw = popt[0]
err_hw = np.sqrt(np.diag(pcov)[0])
title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
ax.plot(m['k'],y_exp,'r',label = title)
ax.legend()

print(f'H = {hw:.2f} ± {err_hw:.2f}')

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}Efk_hw_{hw_txt}_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute spatial attenuation coefficient 

time_att = np.column_stack((m['lambda'],m['err_lambda']))
water_depth = [hw , err_hw]
alpha,err_alpha = time2space_attenuation(time_att, m['k'], water_depth)

# fit by a powerlaw
log_alpha = np.log(alpha)
log_f = np.log(m['f'])
p,cov = np.polyfit(log_f,log_alpha,1,cov = True)
err_coeff = np.sqrt(np.diag(cov))
yth = np.exp(np.polyval(p,log_f))
label_th = r'$y = ' + f'{np.exp(p[1]):.2f}'+ 'f^{' + f'{p[0]:.3f}' + '}$'

fig, ax = plt.subplots()
ax.errorbar(m['f'],alpha,yerr = err_alpha,fmt = 'o')
ax.plot(m['f'],yth,linewidth = 2,label = label_th)
ax.set_title('Spatial attenuation')
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.1,1])
ax.set_ylim([1e-3,1e-1])
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_from_temporal_attenuation_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  

# Save data 
m['alpha'] = alpha
m['err_alpha'] = err_alpha
m['power_law'] = {}
m['power_law']['B'] = np.exp(p[1])
m['power_law']['err_B'] = np.exp(p[1])*err_coeff[1]
m['power_law']['beta'] = p[0]
m['power_law']['err_beta'] = err_coeff[0]

file2save = f'{save_folder}attenuation_from_Efk_respect2f_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)
#%% Create main_results structure and save it 

main_results = {}
main_results['date'] = date
main_results['drone_ID'] = drone_ID
main_results['exp_ID'] = exp_ID
main_results['DRONE'] = data['DRONE']
main_results['GPS'] = data['GPS']
main_results['t0_UTC'] = data['t0_UTC']

main_results['FFT_spectrum'] = FFT_spectrum
main_results['attenuation'] = m

file2save = f'{save_folder}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)






















