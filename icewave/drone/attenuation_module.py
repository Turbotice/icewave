# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 20:24:03 2026

@author: sebas

Python module gathering all functions useful to determine waves attenuation
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 

import scipy
import math

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


def matrix_peak_correlation_detection(M,x,y,x_range,y_range,model_signal,detec_param):
    """ Perform peak detection using gaussian correlation 
    Inputs : - M, numpy.ndarray [ny,nx]: matrix over which we perform peak detection
             - x, array type (nx,) : detection is performed at fixed x
             - y, array type (ny,): correlation is performed with respect to y coordinate
             - x_range,y_range, array (2,) : range of x and y over which we look for maxima
             - model_signal, array type (ny,): reference array used to perform correlation
             - detec_param, dict : parameters dictionnary to be used for find_peaks """
    
    delta_y = np.diff(y)[0]
    
    print(M.shape)
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

def attenuation_from_FK_spectrum(E,peaks_properties,x_range,y_range,signal_model,rel_height,sub_folder):
    """ Compute peak width using peaks detected from an FK spectrum
    Inputs: - E, dictionnary with following keys: 
                    - E: array like, [nx,ny], FK spectrum
                    - x: array like, [nx,], x coordinates values
                    - y: array like, [ny,], y coordinates values
            - peaks_properties, dictionnary with following keys:
                    - peaks: array like, indices of detected peaks
                    - x: array like, x coordinate of detected peaks
                    - y: array like, y coordinate of detected peaks 
            - x_range,y_range : arrays like (2,) range of coordinate values over which we try to fit each detected peak
            - signal_model : array like, model signal used for convolution in lorentzian fit
            - rel_height, float, relative height used to compute width, see function indices2fit 
            - sub_folder, string, folder where figures of each fit will be saved 
    Output : - m, dictionnary, contains all results 
    """
    mask_x = np.where(np.logical_and(E['x'] < x_range[1],E['x'] > x_range[0]))[0]
    mask_y = np.where(np.logical_and(E['y'] < y_range[1],E['y'] > y_range[0]))[0]
    
    x = E['x'][mask_x]
    y = E['y'][mask_y]
    M = E['E'][np.ix_(mask_x,mask_y)]
    
    m = {'A':[],'x':[],'y':[],'alpha':[],'err_x':[],'err_alpha':[],'d':[]}
        
    for idx, current_y in enumerate(y):
        
        # select peaks detected at current wavevector
        mask = np.where(peaks_properties['y'] == current_y)[0]
        peaks = peaks_properties['peaks'][mask]
        
        if len(peaks) > 0:
        
            y_txt = f'{current_y:.2f}'.replace('.','p')
            figname = f'{sub_folder}Lorentzian_fit_y_{y_txt}'
            
            # compute attenuation 
            current_results = fit_peaks_by_lorentzian(M[:,idx], x, peaks, signal_model, rel_height, figname)
            m['A'].append(current_results['A'])
            m['x'].append(current_results['x0'])
            # m['x'].append(x[peaks])
            m['err_x'].append(current_results['err_x0'])
            m['y'].append(current_y)
            m['alpha'].append(current_results['alpha'])
            m['err_alpha'].append(current_results['err_alpha'])
            m['d'].append(current_results['d'])
        
    for key in m.keys():
        m[key] = np.array(m[key])
        m[key] = np.squeeze(m[key])
        
    return m

#---------------------------------------------------------------------------------------------------------------------


def temporal_attenuation_from_FK_spectrum(Efk,peaks_properties,f_range,k_range,signal_model,rel_height,sub_folder):
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
    mask_f = np.where(np.logical_and(Efk['f'] < f_range[1],Efk['f'] > f_range[0]))[0]
    mask_k = np.where(np.logical_and(Efk['k'] > k_range[0],Efk['k'] < k_range[1]))[0]

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
            current_results = fit_peaks_by_lorentzian(M[:,idx], f*2*np.pi, peaks, signal_model, rel_height, figname)
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

