# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 16:17:46 2025

@author: sebas

Extract local maxima of hydro-elastic dispersion relation curves obtained for different positions.
This can be used to perform inversion of Young modulus and thickness coefficients using both active and passive sources 

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

#%% Function section 

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


def process_data(data,fig_folder,date,swell_dict,selected_x):
    
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

    # Extract dispersion relation points

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
    
    freq_range = peaks_detec['freq_range']
    min_prominence = peaks_detec['min_prominence']
    min_width = peaks_detec['min_width']
    
    list_idx_x = [np.argmin(abs(s - x)) for x in selected_x]

    disp_dict = {}
    disp_dict['freq'] = {}
    disp_dict['k'] = {}
    disp_dict['xpos'] = {}

    for i,idx_x in enumerate(list_idx_x):
        FK = abs(mean_CWT[:,:,idx_x])
    
        k_exp,f_exp = extract_peaks(FK,freq,k,freq_range,min_prominence,min_width)
        
        current_key = str(selected_x[i])
        print(current_key)
        disp_dict['freq'][current_key] = f_exp
        disp_dict['k'][current_key] = k_exp
        disp_dict['xpos'][current_key] = s[idx_x]
        
        
    disp_dict['UTC_0'] = label_UTC0
    
    file2save = f'{current_fig_folder}{date}_hydro_elastic_disp_relation_subpix_{label_swell}_file_{label_UTC0}.h5'
    rw.save_dict_to_h5(disp_dict, file2save)
    print(f'{file2save} successfully saved !')
    

def main():
    
    date = '0211'
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
    
    # load swell orientation 
    file_swell_orientation = f'{main_path}Data/swell_orientation_beam_forming.h5'
    swell_orientation = rw.load_dict_from_h5(file_swell_orientation)
    
    # load Ludovic's detection 
    save_path = os.path.join(path2data, f'new_disp_curv_{date}.h5')
    print(save_path)
    ludo_disp = rw.load_dict_from_h5(save_path)
    
    xpos_dict = ludo_disp['xposition'] # dict
    selected_x = np.array([xpos_dict[key] for key in xpos_dict.keys()])
    print(selected_x)
    
    for file2load in filelist :
        data = rw.load_dict_from_h5(file2load)
        
        swell_dict = {'bool':0, 'theta':swell_orientation[date]}
        process_data(data,fig_folder,date,swell_dict,selected_x)
        print(f'D extracted without swell correction, file {file2load}')
        
        swell_dict = {'bool':1, 'theta':swell_orientation[date]}
        process_data(data,fig_folder,date,swell_dict,selected_x)
        print(f'D extracted with swell correction, file {file2load}')
        
if __name__ == '__main__':
    main()
