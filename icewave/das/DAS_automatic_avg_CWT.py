# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 09:53:22 2025

@author: sebas
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

from datetime import datetime
import pytz

import scipy
import pywt

import sys
# sys.path.append('C:/Users/sebas/git')
# sys.path.append('/media/turbots/DATA/thiou/homes/skuchly/git/')

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

# plt.rcParams.update({
#     "text.usetex": False}) # use latex

global g
g = 9.81
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10


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
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
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

def all_chunks_FFT_t(stack_strain,fs):
    """ Perform Fourier Transform in time for all time chunks in a given h5 file. 
    Inputs : - stack_strain, numpy array, dimensions [nb_chunks,nt,nx], contains spatio temporal signals
             - fs, float, sampling frequency 
    Outputs : - FFT_t, numpy array, [nb_chunks,nf,nx] Time Fourier transform of each spatio-temporal chunk 
              - freq, numpy array, dimension (nf,) contains frequencies """

    N = 2**(FT.nextpow2(stack_strain.shape[1]))
    
    mean_sig = np.mean(stack_strain,axis = 1)
    mean_sig = np.tile(mean_sig,(stack_strain.shape[1],1,1)).transpose((1,0,2))
    
    FFT_t = np.fft.fft(stack_strain - mean_sig,n = N,axis = 1)
    FFT_t = FFT_t/stack_strain.shape[1]
    
    # keep only half of frequencies
    FFT_t = FFT_t[:,:N//2,:]
    FFT_t[:,1:-1,:] = 2*FFT_t[:,1:-1,:]
    
    freq = fs*np.arange(0,(N/2))/N
    
    return FFT_t, freq

#-------------------------------------------------------------------------------------------------------

def plot_time_spectrum_chunk(FFT_t,freq,s,UTC_stack,fig_folder,colormap): 
    """ Plot time spectrum for all time chunks 
    Input: - FFT_t, np array, [nb_chunks,nf,nx], Time Fourier transform array
           - freq, np array, (nf,) array of frequencies
           - s, np array, (nx,) array of curvilinear coordinates
           - fig_folder, string, folder where plots will be saved 
           - colormap, colormap object used for plot """

    for chunk in range(FFT_t.shape[0]):
        
        set_graphs.set_matplotlib_param('single')
        fig, ax = plt.subplots()
        imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
                  extent = [freq[0],freq[-1],s[0],s[-1]])
        imsh.set_clim([1e1,3e2])
        
        ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
        ax.set_ylabel(r'$s \; \mathrm{(m)}$')
        ax.set_xlim([0,10])
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh,cax = cax)
        
        title_label = f'UTC : {UTC_stack[chunk,0]}'
        ax.set_title(title_label)
        
        format_date = '%Y-%m-%dT%H-%M-%S'
        label_UTC0 = UTC_stack[chunk,0].strftime(format_date)
        figname = f'{fig_folder}Time_spectrum_UTC0_{label_UTC0}'
        
        plt.savefig(f'{figname}.png',bbox_inches = 'tight')
        plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
        plt.close(fig)
       
#-----------------------------------------------------------------------------------------------------------------

def CWT_all_chunks(FFT_t,scales,wavelet,sampling_period):
    """ Compute Continuous Wavelet Transform over space dimension for each time chunk. See pywt.cwt documentation for more
    precision. This function returns the average over all time chunks of the absolute value of each CWT
    Inputs: - FFT_t, np array, [nb_chunks,nf,nx] Time Fourier Transform for each time chunk and each position
            - scales, np array 1D, scales that will be used during CWT computation 
            - wavelet, string, wavelet used for CWT computation 
            - sampling_period, float, sampling period over space 
    Outputs: - mean_CWT, np array [nf,nk,nx], space-time spectrum averaged over all time chunks
             - k, np array (nk,), array of wavevectors """


    CWT = np.zeros((FFT_t.shape[0],FFT_t.shape[1],len(scales),FFT_t.shape[2]),dtype = 'complex')
    # 0 : chunk_idx, # 1 : freq_idx, # 2 : scales, # 3 : space (curvilinear coordinate)
    
    # loop over all chunks
    for current_chunk in range(FFT_t.shape[0]):
        print(current_chunk)
        current_FFT = FFT_t[current_chunk,:,:]
        
        for i in range(current_FFT.shape[0]): # loop over time frequencies 
            cwtmatr, sigmas = pywt.cwt(current_FFT[i,:], scales, wavelet, sampling_period = sampling_period)
            CWT[current_chunk,i,:,:] = cwtmatr
    
    # Average norm of CWT over all chunks 
    mean_CWT = np.mean(abs(CWT),axis = 0)
    
    # compute wavevector 
    k = 2*np.pi*sigmas

    print('CWT computed successfully !')
    return mean_CWT,k


#%% Load DAS parameter and Data


def process_file(file2load,fig_folder,save_path,path2DAS_param,date): 
    """ Process a file, computation of Time Fourier Transform, Continuous Wavelet Transform, 
    average of FK spectrum and extraction of flexural modulus 
    Inputs : - file2load, string, path to .h5 file 
             - fig_folder, string, folder where plots are saved 
             - save_path, string, folder where avg_CWT is saved
             - date, string, date of acquisition """
             
    fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)
    
    Nb_minutes = 1 # duration of each stack
    stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    
    # decimation for 0210 and 0212
    if date in date_downsampling:
        fs = fs/down_sampling_factor # new value of sampling frequency
        stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
                                                                            stack_time,UTC_stack,down_sampling_factor)
        print(f'New value of sampling frequency, fs = {fs:.2f}')
    
    format_date = '%Y-%m-%dT%H-%M-%S'
    label_UTC0 = UTC_stack[0,0].strftime(format_date)
    current_fig_folder = f'{fig_folder}Wavelet_study_{label_UTC0}/'
    if not os.path.isdir(current_fig_folder):
        os.mkdir(current_fig_folder)
        
    # Perform Fourier transform in time for all time chunks
    
    FFT_t, freq = all_chunks_FFT_t(stack_strain,fs)
    plot_time_spectrum_chunk(FFT_t,freq,s,UTC_stack,current_fig_folder,parula_map)
    
    # Define wavelet 
    
    wavelet = 'cmor1.5-1.5' # mother wavelet 
    [psi, x] = pywt.ContinuousWavelet(wavelet).wavefun(10)
    
    fig, axs = plt.subplots(figsize = (6,4))
    axs.plot(x,np.real(psi))
    axs.plot(x,np.imag(psi))
    axs.set_xlim([-5,5])
    
    figname = f'{current_fig_folder}wavelet_{wavelet}'
    plt.savefig(f'{figname}.png',bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
    plt.close(fig)
    
    # Perform CWT for all time chunks 
    
    # parameters for CWT
    sampling_period = 1/facq_x #spatial sampling wavelength
    #logarithmic scales
    scales = np.geomspace(2,1024,num = 100)
    
    # compute normalized wavenumber (sigma)
    norm_sigma = pywt.scale2frequency(wavelet,scales)
    wavelength_sampled = sampling_period/norm_sigma # typical wavelength that will be sampled 
    print(f'Wavelengths sampled : min = {wavelength_sampled[0]:.2f}, max = {wavelength_sampled[-1]:.2f}')
    
    mean_CWT,k = CWT_all_chunks(FFT_t,scales,wavelet,sampling_period)
    
    # save mean_CWT, frequencies and apparent_wavevectors
    CWT_results = {}
    CWT_results['mean_CWT'] = mean_CWT
    CWT_results['k_star'] = k
    CWT_results['f'] = freq
    CWT_results['x'] = s
    CWT_results['wavelet'] = wavelet
    CWT_results['scales'] = scales
    CWT_results['label'] = label_UTC0
    CWT_results['DAS_param'] = {'fs':fs,'fiber_length':fiber_length,'fx':facq_x}

    file2save = f'{save_path}avg_CWT_{label_UTC0}.h5'
    rw.save_dict_to_h5(CWT_results, file2save)
   
    print(f'{file2save} successfully saved !')
    
    
def process_file_from_args(args):
    
    file2load,fig_folder,save_path,path2DAS_param,date = args
    process_file(file2load, fig_folder, save_path, path2DAS_param, date)
    
    

def main():
    date = '0211'
    print(np.__version__)

    # Load parameters for DAS
    main_path = '/media/turbots/Backup25/'
    path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'
    
    path2data = f'{main_path}Data/{date}/DAS/'
    filelist = glob.glob(path2data + '*UTC.h5')

    # Create folder for saving graphs
    fig_folder = f'{path2data}Figures/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
        
    save_path = f'{path2data}avg_CWT/'
    if not os.path.isdir(save_path):
        os.mkdir(save_path)
    
    # perform parfor loop
    # args_list = [(filelist[idx_file],fig_folder,save_path,path2DAS_param,date) for 
    #              idx_file in range(len(filelist) - 1)]
    
    # print(args_list)
    
    # with ProcessPoolExecutor(max_workers = 20) as executor:
    #     list(tqdm(executor.map(process_file_from_args, args_list), total=len(args_list)))
    
    for idx_file in range(len(filelist)-1):
        file2load = filelist[idx_file]
        print(file2load)
        process_file(file2load, fig_folder, save_path, path2DAS_param, date)
        
if __name__ == '__main__':
    main()
