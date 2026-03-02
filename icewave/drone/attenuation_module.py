# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 20:24:03 2026

@author: sebas

Python module gathering all functions useful to determine waves attenuation
from an FK spectrum 
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
        print(m[key])
        m[key] = np.array(m[key])
        m[key] = np.squeeze(m[key])
        
    return m

#---------------------------------------------------------------------------------------------------------------------

def plot_FK_spectrum(Efk):
    """ Plot FK spectrum """
    
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
    
    return fig, ax, c, cbar

#-------------------------------------------------------------------------------------------

def extract_peaks_fixed_k(Efk,gaussian,wavevector_range,frequency_range,file2save,detec_param):
    """ Extract peaks of the wave FK spectrum using fixed values of wavevector k, 
    peaks are detected using fit over frequencies. Detected peaks are then filtered using scatter_filter.py module. 
    Inputs : - Efk, dictionary, with keys: + 'E' = array like [nf,nk], FK spectrum
                                           + 'f' = array like, frequencies
                                           + 'k' = array like, wavevectors
             - gaussian, array like, gaussian signal used to perform convolution over frequencies and detect peaks
             - wavevector_range, array like, range of wavevectors over which we will detect peaks
             - frequency_range, array like, range of frequencies over which we will detect peaks 
             - file2save, string, path where detected peaks properties will be saved 
             - detec_param, dictionnary, collects parameters used to detect maximum convolution peaks
             
    Outputs : - filtered properties, dict, contains peaks properties, it is also saved as a .h5 file"""
    
    Speaks = matrix_peak_correlation_detection(Efk['E'], Efk['k'], Efk['f'], wavevector_range, frequency_range, 
                                          gaussian, detec_param)
    # change key name
    Speaks['k'] = Speaks['x']
    Speaks['f'] = Speaks['y']
    del Speaks['x']
    del Speaks['y']
    
    print('Filter undesired points')
    # Filter detected points, plot results and save filtered peaks
    filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(Speaks['k'], Speaks['f'],Speaks)

    # plot detected peaks
    fig, ax, c, cbar = plot_FK_spectrum(Efk)
    ax.plot(filtered_properties['k'],filtered_properties['f'],'r.')
    
    # save filtered_properties
    rw.save_dict_to_h5(filtered_properties, file2save)
    print('DONE.')
    
    return filtered_properties

#-------------------------------------------------------------------------------------------------------

def extract_peaks_fixed_f(Efk,gaussian,wavevector_range,frequency_range,file2save,detec_param):
    """ Extract peaks of the wave FK spectrum using fixed values of frequency f, 
    peaks are detected using fit over wavevectors. Detected peaks are then filtered using scatter_filter.py module. 
    Inputs : - Efk, dictionary, with keys: + 'E' = array like [nf,nk], FK spectrum
                                           + 'f' = array like, frequencies
                                           + 'k' = array like, wavevectors
             - gaussian, array like, gaussian signal used to perform convolution over frequencies and detect peaks
             - wavevector_range, array like, range of wavevectors over which we will detect peaks
             - frequency_range, array like, range of frequencies over which we will detect peaks 
             - file2save, string, path where detected peaks properties will be saved 
             - detec_param, dictionnary, collects parameters used to detect maximum convolution peaks
             
    Outputs : - filtered properties, dict, contains peaks properties, it is also saved as a .h5 file"""
    
    Speaks = matrix_peak_correlation_detection(Efk['E'].T, Efk['f'], Efk['k'], frequency_range, wavevector_range, 
                                               gaussian, detec_param)
    # change key name
    Speaks['f'] = Speaks['x']
    Speaks['k'] = Speaks['y']
    del Speaks['x']
    del Speaks['y']
    
    print('Filter undesired points')
    filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(Speaks['k'], Speaks['f'],Speaks)
    
    fig,ax,c,cbar = plot_FK_spectrum(Efk)
    ax.plot(filtered_properties['k'],filtered_properties['f'],'r.')
    
    # save filtered_properties
    rw.save_dict_to_h5(filtered_properties, file2save)
    print('DONE.')

    return filtered_properties

#----------------------------------------------------------------------------------------

def temporal_attenuation(Efk,filtered_properties,frequency_range,wavevector_range,
                         gaussian,detec_param,fig_folder):
    """ Compute temporal attenuation from FK spectrum, based on a Lorentzian fit over 
    frequencies. 
    
    Inputs : - Efk, dictionary, with keys: + 'E' = array like [nf,nk], FK spectrum
                                           + 'f' = array like, frequencies
                                           + 'k' = array like, wavevectors
             - filtered_properties, dictionary, properties of detected peaks of FK spectrum, with keys:
                 + 'peaks' = indices of detected peaks
                 + 'f' = array like, frequencies
                 + 'k' = array like, wavevectors
             - gaussian, array like, gaussian signal used to perform convolution over frequencies and detect peaks
             - wavevector_range, array like, range of wavevectors over which we will detect peaks
             - frequency_range, array like, range of frequencies over which we will detect peaks 
             - detec_param, dictionnary, collects parameters used to detect maximum convolution peaks
             - fig_folder, string, path where figures of fit are saved 
             
    Outputs : - m, dictionary, contains results of temporal attenuation """
    
    # get time 
    E = {}
    E['E'] = Efk['E']
    E['x'] = 2*np.pi*Efk['f']
    E['y'] = Efk['k']

    peaks_dico = {}
    peaks_dico['peaks'] = filtered_properties['peaks']
    peaks_dico['x'] = 2*np.pi*filtered_properties['f']
    peaks_dico['y'] = filtered_properties['k']

    x_range = 2*np.pi*np.array(frequency_range)
    y_range = wavevector_range

    rel_height = detec_param['rel_height']
    sub_folder = f'{fig_folder}time_lorentzian_fit/'
    if not os.path.isdir(sub_folder):
        os.mkdir(sub_folder)
        
    m = attenuation_from_FK_spectrum(E, peaks_dico, x_range, y_range, gaussian, rel_height, sub_folder)

    del E
    del peaks_dico

    m['f'] = m.pop('x')/2/np.pi
    m['err_f'] = m.pop('err_x')/2/np.pi
    m['k'] = m.pop('y')
    m['lambda'] = m.pop('alpha')
    m['err_lambda'] = m.pop('err_alpha')
    print(m.keys())
    
    return m

#------------------------------------------------------------------------------------------------

def spatial_attenuation(Efk,filtered_properties,gaussian,wavevector_range,frequency_range,
                        detec_param,fig_folder):
    """ Compute space attenuation from FK spectrum, based on a Lorentzian fit over 
    wavevectors. 
    
    Inputs : - Efk, dictionary, with keys: + 'E' = array like [nf,nk], FK spectrum
                                           + 'f' = array like, frequencies
                                           + 'k' = array like, wavevectors
             - filtered_properties, dictionary, properties of detected peaks of FK spectrum, with keys:
                 + 'peaks' = indices of detected peaks
                 + 'f' = array like, frequencies
                 + 'k' = array like, wavevectors
             - gaussian, array like, gaussian signal used to perform convolution over frequencies and detect peaks
             - wavevector_range, array like, range of wavevectors over which we will detect peaks
             - frequency_range, array like, range of frequencies over which we will detect peaks 
             - detec_param, dictionnary, collects parameters used to detect maximum convolution peaks
             - fig_folder, string, path where figures of fit are saved 
             
    Outputs : - m, dictionary, contains results of temporal attenuation """
    
    E = {}
    E['E'] = Efk['E'].T
    E['x'] = Efk['k']
    E['y'] = 2*np.pi*Efk['f']


    peaks_dico = {}
    peaks_dico['peaks'] = filtered_properties['peaks']
    peaks_dico['y'] = 2*np.pi*filtered_properties['f']
    peaks_dico['x'] = filtered_properties['k']

    y_range = 2*np.pi*np.array(frequency_range)
    x_range = wavevector_range

    rel_height = detec_param['rel_height']
    sub_folder = f'{fig_folder}space_lorentzian_fit/'
    if not os.path.isdir(sub_folder):
        os.mkdir(sub_folder)
        
    m = attenuation_from_FK_spectrum(E, peaks_dico, x_range, y_range, 
                                             gaussian, rel_height, sub_folder)

    del E
    del peaks_dico

    m['k'] = m.pop('x')
    m['err_k'] = m.pop('err_x')
    m['f'] = m.pop('y')/2/np.pi
    print(m.keys())
    
    return m

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

#----------------------------------------------------------------------------------------


def structure_space_attenuation(m):
    """ Add space attenuation computed from temporal attenuation to an existing structure m """
    
    time_att = np.column_stack((m['lambda'],m['err_lambda'])) # time attenuation coefficient and standard deviation
    water_depth = [m['hw'] , m['err_hw']] # water depth and standard deviation
    alpha,err_alpha = time2space_attenuation(time_att, m['k'], water_depth) # compute spatial attenuation

    # fit by a powerlaw
    coeffs,err_coeffs = powerlaw_fit(m['f'], alpha, err_alpha)
    
    # Save data 
    m['alpha'] = alpha
    m['err_alpha'] = err_alpha
    m['power_law'] = {}
    m['power_law']['B'] = coeffs[1]
    m['power_law']['err_B'] = err_coeffs[1]
    m['power_law']['beta'] = coeffs[0]
    m['power_law']['err_beta'] = err_coeffs[0]

    return m

#----------------------------------------------------------------------------------------------

def plot_attenuation_power_law(m,xbounds,ybounds,figname):
    """ Plot power law fit of attenuation law alpha VS frequency """
    
    x_fit = np.linspace(xbounds[0],xbounds[1],100)

    # get powerlaw fit
    coeffs = (m['power_law']['beta'],m['power_law']['B'])
    yth = coeffs[1]*x_fit**coeffs[0]

    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    
    list_key = list(m.keys())
    if 'err_f' not in list_key:
        m['err_f'] = None

    fig, ax = plt.subplots()
    ax.errorbar(m['f'],m['alpha'],yerr = m['err_alpha'],xerr = m['err_f'],fmt = 'o')
    ax.plot(x_fit,yth,linewidth = 2,label = label_th)
    ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
    ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xbounds)
    ax.set_ylim(ybounds)
    ax.legend()

    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')  






# def temporal_attenuation_from_FK_spectrum(Efk,peaks_properties,f_range,k_range,signal_model,rel_height,sub_folder):
#     """ Compute temporal attenuation of detected peaks from a FK spectrum. 
#     Inputs : - Efk, dictionnary, contains follwoing keys :+ E, array like [nf,nk], space-time spectrum 
#                                                           + f, array like (nf,), frequency array 
#                                                           + k, array like (nk,), wavevector array 
#              - peaks_properties, dictionnary, containing at least keys : 'f','k' and 'peaks'. Contains peaks indices 
#              and their associate coordinates (f,k). 
#              - rel_height, float, relative height used to compute width, see function indices2fit 
#              - sub_folder, string, folder where figures of each fit will be saved 
#     Output : - m, dictionnary, contains all results """
    
    
#     # select area of FK where we want to compute attenuation
#     mask_f = np.where(np.logical_and(Efk['f'] < f_range[1],Efk['f'] > f_range[0]))[0]
#     mask_k = np.where(np.logical_and(Efk['k'] > k_range[0],Efk['k'] < k_range[1]))[0]

#     f = Efk['f'][mask_f]
#     k = Efk['k'][mask_k]
#     M = Efk['E'][np.ix_(mask_f,mask_k)]
    
#     m = {'A':[],'f':[],'k':[],'lambda':[],'err_f':[],'err_lambda':[],'d':[]}
        
#     for idx, current_k in enumerate(k):
        
#         # select peaks detected at current wavevector
#         mask = np.where(peaks_properties['k'] == current_k)[0]
#         peaks = peaks_properties['peaks'][mask]
        
#         if len(peaks) > 0:
        
#             k_txt = f'{current_k:.2f}'.replace('.','p')
#             figname = f'{sub_folder}Lorentzian_fit_respectto_f_fixed_k_{k_txt}'
            
#             # compute attenuation 
#             current_results = fit_peaks_by_lorentzian(M[:,idx], f*2*np.pi, peaks, signal_model, rel_height, figname)
#             m['A'].append(current_results['A'])
#             m['f'].append(current_results['x0']/2/np.pi)
#             m['err_f'].append(current_results['err_x0']/2/np.pi)
#             m['k'].append(current_k)
#             m['lambda'].append(current_results['alpha'])
#             m['err_lambda'].append(current_results['err_alpha'])
#             m['d'].append(current_results['d'])
        
#     for key in m.keys():
#         m[key] = np.array(m[key])
#         m[key] = np.squeeze(m[key])
        
#     return m

