# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 09:37:55 2025

@author: sebas

Extract flexural modulus D from averaged Continuous Wavelet Transform. 
"""


import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors 
import matplotlib.cm as cm

import h5py 
import glob
import os 
import pickle

from datetime import datetime
import pytz

import scipy
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pywt

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

# plt.rcParams.update({
#     "text.usetex": True}) # use latex

global g
g = 9.81
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 

def plot_spatio_temp(spatio,fiber_length,extents,cmap):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsh, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots()
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = cmap,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#----------------------------------------------------------------------------------------------------

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega

#-------------------------------------------------------------------------------------------------------

def get_water_height(DAS_water_height,UTC_t,selected_x):
    """ Compute water height for a given UTC datetime and position along the fiber
    Input : - DAS_water_height, dictionnary containing recording of water depth evolution during the whole experiment 
            - UTC_t, datetime object, time at which we look for water height
            - selected_x, float, position along optical fiber in meter 
    Output : H, float, water height at the selected time and position 
    """

    closest_idx = np.argmin(abs(UTC_t - DAS_water_height['UTC_t']))
    closest_idx_pos = np.argmin(abs(selected_x - DAS_water_height['s']))
    H = DAS_water_height['water_height'][closest_idx_pos,closest_idx]
    return H

#-------------------------------------------------------------------------------------------------------

def subpixel_polyfit(x,y,p,deg = 2):
    """ Compute subpixel position of a local maximum / minimum 
    Inputs : - x, array of abscisse 
             - y, array of y-axis values
             - p, index at which a peak has been detected 
    Outputs : - local_argmax, subpixel value for which the maximum of degree 2 polynomial fit is reached 
              - local_max, maximum of the degree 2 polynomial """
    
    x2fit = np.array([x[p-1],x[p],x[p+1]])
    y2fit = np.array([y[p-1],y[p],y[p+1]])
    coeffs = np.polyfit(x2fit,y2fit,deg = 2)
    local_argmax = -coeffs[1]/2/coeffs[0]
    local_max = np.polyval(coeffs,local_argmax)

    return local_argmax,local_max

#-------------------------------------------------------------------------------------------------------

def extract_peaks(FK,freq,k,freq_range,min_prominence,min_width):
    """ Extract peak coordinates from a space-time spectrum 
    Inputs : - FK, numpy 2D array [nf,nk], space-time spectrum
             - freq, numpy array (nf,), frequencies
             - k, numpy array (nk,), wavevectors
             - freq_range, tuple [freq_min, freq_max], frequency range within we look for peaks
             - min_prominence, minimum prominence of peaks once a row is normalized (see scipy.signal.findpeaks)
             - min_width, minimal width of peaks once a row is normalized (see scipy.signal.findpeaks)
             
    Outputs : - k_exp, array of peaks wavevector coordinate
              - f_exp, array of peaks frequency coordinate """
    
    freq_min = freq_range[0]
    freq_max = freq_range[1]
    idx_min = np.argmin(abs(freq - freq_min))
    idx_max = np.argmin(abs(freq - freq_max))

    indices = np.arange(idx_min,idx_max + 1, step = 1)

    # create arrays for peak detection 
    k_exp = []
    f_exp = []
    
    for idx_freq in indices:
    # print(f'frequency = {freq[idx_freq]:.2f} Hz')
        current_max = np.max(abs(FK[idx_freq,:]))
        normalized = abs(FK[idx_freq,:])/current_max
        # find peaks
        peaks,properties = scipy.signal.find_peaks(normalized,
                                                   prominence = min_prominence,width = min_width)
        if len(peaks) != 0:
            for peak in peaks:
                local_argmax,local_max = subpixel_polyfit(k,normalized,peak)
                k_exp.append(local_argmax)
                f_exp.append(freq[idx_freq])
    
    k_exp = np.array(k_exp)
    f_exp = np.array(f_exp)

    return k_exp,f_exp

#-------------------------------------------------------------------------------------------------------------------

def extract_D(mean_CWT,selected_x,s,freq,k,peaks_detec,DAS_water_height,UTC_mid):
    """ Compute flexural modulus D from FK plot. First peaks of maximum amplitude are detected
    in the coordinate system (f,k). Then the detected points are fitted by a shallow hydroelastic 
    dispersion relation. Flexural modulus is computed using a least square method fit. 
    Inputs : - mean_CWT, np array [nf,nk,nx], space-time spectrum averaged over all time chunks. COmputed using 
    COntinuous Wavelet Transform. 
             - selected_x, np array 1D, curvilinear coordinates at which D is computed 
             - s, np array 1D (nx,), curvilinear coordinates
             - freq, np array (nf,), frequencies array
             - k, np array (nk,), wavevectors, array
             - peaks_detec, dictionnary, contains properties for peak detection.
             The dictionnary must contain following keys : 
                     + freq_range, tuple, (f_min,f_max), frequency range over which the hydroelastic dispersion relation
                     is fitted 
                     + min_prominence, float, minimal prominence of peaks
                     + min_width, float, minimal width of peaks 
             - DAS_water_height, dictionnary, contains data relative to water depath below the fiber
             - UTC_mid, datetime object, date at which we compute the water height profile below the fiber
    Outputs : - D, np array 1D, values of flexural modulus 
              - err_D, np array 1D, values of standard deviation for each position """


    freq_range = peaks_detec['freq_range']
    min_prominence = peaks_detec['min_prominence']
    min_width = peaks_detec['min_width']
    
    list_idx_x = [np.argmin(abs(s - x)) for x in selected_x]
    
    D = np.zeros(len(selected_x))
    err_D = np.zeros(len(selected_x))
    
    rho_w = 1027

    for i,idx_x in enumerate(list_idx_x):
        FK = abs(mean_CWT[:,:,idx_x])
    
        k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
        
        if len(k_exp) != 0:
            
            H = get_water_height(DAS_water_height,UTC_mid,selected_x[i])
            popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                                 bounds = (1e5,1e8))
            err_coeff = np.sqrt(np.diag(pcov))
            # print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
            D[i] = popt[0]
            err_D[i] = err_coeff[0]
            
        else:
            D[i] = None
            err_D[i] = None

    return D, err_D

#--------------------------------------------------------------------------------------------------------------------

def plot_FK_panel(data,selected_x):
    """ Plot FK for different positions into a 3 x 2 subplot panel 
    Inputs : - data, dictionnary
             - selected_x, selected positions in meter for which we want to show FK plot """
             
    # Plot (f,k) for different positions 
    s = data['x']
    k = data['k_star']
    freq = data['f']
    mean_CWT = data['mean_CWT']
    
    idx_x = [np.argmin(abs(s - x)) for x in selected_x]

    # plot FK without theory
    fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
    for ax, idx in zip(axs.flatten(),idx_x):
        imsh = ax.pcolormesh(k,freq,abs(mean_CWT[:,:,idx]),cmap = parula_map)
        ax.set_ylim([0,1])
        ax.set_xlim([0.01,0.6])
        # ax.set_ylim([0,10])
        # ax.set_xlim([0,1.0])
        ax.set_title(f'x = {s[idx]:.2f} m')
        # imsh.set_clim([0,300])

    for j in range(axs.shape[1]):
        axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

    for i in range(axs.shape[0]):
        axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$')
    
    return fig, axs

#----------------------------------------------------------------------------------------------------------------------

def plot_FK_panel_with_theory(data,selected_x,peaks_detec,date,DAS_water_height):
    """ Plot FK for different positions into a 3 x 2 subplot panel, superposed with hydroelastic theory
    Inputs : - data, dictionnary
             - selected_x, selected positions in meter for which we want to show FK plot """
             
    freq_range = peaks_detec['freq_range']
    min_prominence = peaks_detec['min_prominence']
    min_width = peaks_detec['min_width']
    
    s = data['x']
    k = data['k_star']
    freq = data['f']
    mean_CWT = data['mean_CWT']
    
    format_date = '%Y-%m-%dT%H-%M-%S'
    UTC_0 = datetime.strptime(data['label'],format_date)
    UTC_0 = UTC_0.replace(tzinfo = pytz.timezone('UTC'))
    
    idx_x = [np.argmin(abs(s - x)) for x in selected_x]

    # select range of k for which we look for a peak
    if date in date_downsampling:
        k_min = 0.05
        k_max = 0.8
        idx_kmin = np.argmin(abs(k - k_min))
        idx_kmax = np.argmin(abs(k - k_max))
        mean_CWT = mean_CWT[:,idx_kmax:idx_kmin,:]
        k = k[idx_kmax:idx_kmin]

    rho_w = 1027

    fig, axs = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True,figsize = (12,8))
    for i,(ax,idx) in enumerate(zip(axs.flatten(),idx_x)):
        FK = mean_CWT[:,:,idx]
        imsh = ax.pcolormesh(k,freq,abs(FK),cmap = parula_map)
        ax.set_ylim([0,1])
        ax.set_xlim([0.01,0.6])
        # ax.set_ylim([0,8])
        # ax.set_xlim([0,0.7])
        ax.set_title(f'x = {s[idx]:.2f} m')
        # imsh.set_clim([0,300])
        
        k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
        ax.plot(k_exp,f_exp,'r.')
        print(f'Detected coordinates :{k_exp}')
        
        if len(k_exp) != 0:
            H = get_water_height(DAS_water_height,UTC_0,selected_x[i])
            popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                                 bounds = (1e5,1e8))
            err_coeff = np.sqrt(np.diag(pcov))
            print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')
        
            k_th = np.linspace(0,1,100)
            omega_th = shallow_hydroelastic(k_th, popt[0], rho_w, H)
            label_th = r'$D = ' + f'{popt[0]:.2e}' + r'$'
            ax.plot(k_th,omega_th/2/np.pi,'r--',label = label_th)
            ax.legend(loc = 'lower right',fontsize = 10)

    for j in range(axs.shape[1]):
        axs[1,j].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

    for i in range(axs.shape[0]):
        axs[i,0].set_ylabel(r'$f \; \mathrm{(Hz)}$') 
    
    return fig, axs

#----------------------------------------------------------------------------------------------------------------

def load_DAS_water_height(file2water):
    DAS_water_height = rw.load_dict_from_h5(file2water)
    print(DAS_water_height.keys())
    # convert UTC string to datetime 
    format_date = '%Y-%m-%dT%H-%M-%S.%f'
    # convert UTC string to datetime 
    DAS_water_height['UTC_datetime'] = []
    for date_txt in DAS_water_height['UTC_t']:
        if date_txt != 'None' :
            datetime_obj = datetime.strptime(date_txt,format_date)
            datetime_obj = datetime_obj.replace(tzinfo = pytz.timezone('UTC'))
        # else :
        #     datetime_obj = None
        DAS_water_height['UTC_datetime'].append(datetime_obj)
        
    DAS_water_height['UTC_t'] = np.array(DAS_water_height['UTC_datetime'])
    del DAS_water_height['UTC_datetime']
    
    return DAS_water_height



#%% Plots and D extraction

def process_data(data,fig_folder,date,DAS_water_height,swell_dict):
    
    set_graphs.set_matplotlib_param('single')
    
    # define label 
    label_UTC0 = data['label']
    
    # change k if we want to account for swell orientation 
    if swell_dict['bool']:
        data['k_star'] = data['k_star']/np.cos(swell_dict['theta']*np.pi/180)
        label_swell = 'swell_corrected'
    else:
        label_swell = ''
    
    current_fig_folder = f'{fig_folder}Wavelet_study_{label_UTC0}/'
    if not os.path.isdir(current_fig_folder):
        os.mkdir(current_fig_folder)
    
    # plot FK for different positions
    # selected_x = np.array([100,180,400,450,500,550]) # 0212
    selected_x = np.array([80,150,300,400,500,550]) # 0211
    
    fig, axs = plot_FK_panel(data,selected_x)
    
    figname = f'{current_fig_folder}{date}_FK_subplots_spectrum_spectrum_avg_{label_swell}_{label_UTC0}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    
    peaks_detec = {}
    peaks_detec['freq_range'] = [0.02,0.8]
    peaks_detec['min_prominence'] = 0.7
    peaks_detec['min_width'] = 3
    
    fig,axs = plot_FK_panel_with_theory(data, selected_x, peaks_detec, date, DAS_water_height)
    
    figname = f'{current_fig_folder}{date}_subplots_D_extraction_subpix_spectrum_avg_{label_swell}_{label_UTC0}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    
    
    # Extract flexural modulus 
    
    # Compute flexural modulus D 
    zone1 = np.arange(60,210,step = 20)
    zone2 = np.arange(280,570,step = 20)
    selected_x = np.concatenate((zone1,zone2))
    
    peaks_detec = {}
    peaks_detec['freq_range'] = [0.02,0.8]
    peaks_detec['min_prominence'] = 0.7
    peaks_detec['min_width'] = 3
    
    mean_CWT = data['mean_CWT']
    s = data['x']
    freq = data['f']
    k = data['k_star']
    
    if date in date_downsampling:
        k_min = 0.05
        k_max = 0.8
        idx_kmin = np.argmin(abs(k - k_min))
        idx_kmax = np.argmin(abs(k - k_max))
        mean_CWT = mean_CWT[:,idx_kmax:idx_kmin,:]
        k = k[idx_kmax:idx_kmin]
    
    format_date = '%Y-%m-%dT%H-%M-%S'
    UTC_0 = datetime.strptime(data['label'],format_date)
    UTC_0 = UTC_0.replace(tzinfo = pytz.timezone('UTC'))
    
    D,err_D = extract_D(mean_CWT, selected_x, s, freq, k, peaks_detec,DAS_water_height,UTC_0)
    # save data 
    results_dict = {}
    results_dict['D'] = D
    results_dict['err_D'] = err_D
    results_dict['x'] = selected_x
    results_dict['UTC_0'] = label_UTC0
    
    file2save = f'{current_fig_folder}{date}_wavelet_flexural_modulus_subpix_{label_swell}_file_{label_UTC0}.h5'
    rw.save_dict_to_h5(results_dict, file2save)
    print(f'{file2save} successfully saved !')
    

def main():
    
    date = '0212'
    # main_path = f'U:'

    main_path = '/media/turbots/Backup25/'
    path2data = f'{main_path}Data/{date}/DAS/'

    path2CWT = f'{path2data}avg_CWT/'
    filelist = glob.glob(f'{path2CWT}avg_CWT*.h5')
    print(filelist)
    
    # Create folder for saving graphs
    fig_folder = f'{path2data}Figures/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
    
    # load water height data
    file2water = glob.glob(f'{path2data}fiber_water_height_GPS_structure_{date}.h5')[0]
    print(file2water)    
    DAS_water_height = load_DAS_water_height(file2water)
    
    # load swell orientation 
    file_swell_orientation = f'{main_path}Data/swell_orientation_beam_forming.h5'
    swell_orientation = rw.load_dict_from_h5(file_swell_orientation)
    
    for file2load in filelist :
        data = rw.load_dict_from_h5(file2load)
        
        swell_dict = {'bool':0, 'theta':swell_orientation[date]}
        process_data(data,fig_folder,date,DAS_water_height,swell_dict)
        print(f'D extracted without swell correction, file {file2load}')
        
        swell_dict = {'bool':1, 'theta':swell_orientation[date]}
        process_data(data,fig_folder,date,DAS_water_height,swell_dict)
        print(f'D extracted with swell correction, file {file2load}')
        
if __name__ == '__main__':
    main()

    