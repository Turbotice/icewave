# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 15:34:29 2025

@author: sebas
"""


import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
from time import strftime, localtime
from datetime import datetime
import pytz
import pickle

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

font_size_medium = 20
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

import sys
sys.path.append('C:/Users/sebas/git')
# '/home/turbots/Codes/git/' path to icewave package on adour 

#%% FUNCTION SECTION 

def epoch2datetime(t,timezone = None,format_date = '%Y-%m-%d %H:%M:%S.f'):
    """ Convert time since epoch to UTC time 
    Inputs: - t, array or list of times since epoch
            - timezone, pytz timezone in which we want datetime 
            - format_date, string, format in which datetime will be written """  
    
    t_datetime = []
    for t_epoch in t:
        current_string = strftime(format_date, localtime(t_epoch))
        current_time = datetime.strptime(current_string, format_date)
        if timezone is None:
            adjusted_time = current_time
        else:
            adjusted_time = current_time.astimezone(timezone)
            
        t_datetime.append(adjusted_time)
    
    t_datetime = np.array(t_datetime)
    t_datetime = np.reshape(t_datetime,(t.shape))
    return t_datetime

#-------------------------------------------------------------------------------------------------------


def generate_datetime_txt(UTC_0):
    """ Generate a string from a datetime UTC_0"""
    
    txt = f'{UTC_0.year}_{UTC_0.month}_{UTC_0.day}_{UTC_0.hour}_{UTC_0.minute}_{UTC_0.second}'
    return txt


#-------------------------------------------------------------------------------------------------------

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#--------------------------------------------------------------------------------------------------------

def time_stacking(strain_rate,nb_seconds,fiber_length):
    """ Stack strain rate measurements over time
        Inputs : - strain_rate, 3D numpy array, #dim 0 : each second of acquisition, #dim 1 : time sampled at fs within 
        the corresponding second, #dim 2 : space
                 - nb_seconds , int type, number of seconds to stack
                 - fiber_length, float type, fiber length
                 
        Outputs : - spatio_long, 2D numpy array, #dim 0 : time in seconds, #dim 1 : space
                  - s, 1D numpy array, curvilinear coordinate
                  - t_sec, 1D numpy array, time in seconds """
                  
    time_length = np.shape(strain_rate)[0]
    fs = np.shape(strain_rate)[1] # time frequency sampling 
    fx = np.shape(strain_rate)[2]/fiber_length # space frequency sampling
    
    # curvilinear coordinate 
    s = np.arange(0,fiber_length,1/fx)
    
    for i in range(nb_seconds):
        if i == 0 :    
            spatio_long = np.vstack((strain_rate[i,:,:],strain_rate[i + 1,:,:]))
        elif i > 1:
            spatio_long = np.vstack((spatio_long,strain_rate[i,:,:]))
            
    print(f'Spatio-temp computed for a total number of seconds : {nb_seconds}')
    t_sec = np.arange(0,nb_seconds,1/fs)  # time in seconds
    
    return spatio_long,s,t_sec

#---------------------------------------------------------------------------------------------------

def new_stack_strain_rate(strain_rate,t,fiber_length,Nb_minutes = 1):
    """ Structure strain rate array so that several minutes of spatio-temporal recording are stacked in a 3D array
    Inputs: - strain_rate, 3D numpy array, #dim 0 : each second of acquisition, #dim 1 : time sampled at fs within 
            the corresponding second, #dim 2 : space
            - t, numpy array, time since epoch for each seconds of recording
            - fiber_length, float, length of the fiber
            - Nb_minutes, float, duration of each created stack
    Outputs: - stack_strain, 3D numpy array, #dim 0 : each stack of duration Nb_minutes, #dim 1 : time sampled at fs within 
            the corresponding stack, #dim 2 : space
            - stack_time, 2D numpy array # dim 0: each stack of duration Nb_minutes, # dim 1: time in seconds, taking 
            beginning of recording as a reference
            - s, numpy array, curvilinear coordinate """

    fs = np.shape(strain_rate)[1] # time frequency sampling 
    fx = np.shape(strain_rate)[2]/fiber_length # space frequency sampling
    
    # create an array # dim 0: stacks of duration Nb_minutes, 
    # dim 1: time samples within Nb_minutes, dim 2: space samples
    Nb_stacks = int(np.shape(strain_rate)[0]/60/Nb_minutes) 
    print(f'Created {Nb_stacks} stacks of duration {Nb_minutes} minutes')  
    Nb_seconds = int(Nb_minutes*60) # nb_seconds for a single section 
    
    stack_strain = np.zeros((Nb_stacks,Nb_seconds*fs,np.shape(strain_rate)[2]))
    stack_time = np.zeros((Nb_stacks,Nb_seconds*fs))
    
    for n in range(Nb_stacks):
        idx_min = int(n*Nb_seconds)
        idx_max = int((n + 1)*Nb_seconds)
        
        current_section = strain_rate[idx_min : idx_max, : ,:]
        spatio_long,s,t_sec = time_stacking(current_section,Nb_seconds,
                                            fiber_length) # stack for a given duration Nb_minutes
        stack_strain[n,:,:] = spatio_long
    
        # create array of time
        current_time = np.arange(n*Nb_seconds, (n + 1)*Nb_seconds, 1/fs)
        stack_time[n,:] = current_time 
        

    return stack_strain,stack_time,s 
    
#----------------------------------------------------------------------------------------------------

# function for plotting spatio temporal
def plot_spatio_temp(spatio,t,s,fiber_length):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - t, numpy array or list, time array 
            - s, numpy array or list, curvilinear coordinate array
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsho, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
              interpolation = 'gaussian', extent = extents(t) + extents(s))
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#----------------------------------------------------------------------------------------------------


#%% Load data   

def main():
    
    date = '0211'
    
    # Load parameters for DAS
    path2DAS_param = 'U:/Data/parameters_Febus_2025.pkl'
    with open(path2DAS_param,'rb') as pf:
        param = pickle.load(pf)
    print('Parameters file loaded')
    
    # Set parameters
    # fs = np.shape(strain_rate)[1] # time sampling frequency 
    # facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
    fs = param[date]['fs']
    fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
    facq_x = param[date]['facq_x'] 
    
    format_date = '%Y-%m-%d %H:%M:%S.f'
    local_timearea = pytz.timezone('America/Montreal')
    UTC_timearea = pytz.timezone('UTC')
    
    # path2data = f'E:/Data/{date}/DAS_h5/'
    # path2data = 'F:/20250211/'
    path2data = f'U:/Data/{date}/DAS/'
    
    # Create folder for saving graphs
    fig_folder = f'{path2data}Figures/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
    
    filelist = glob.glob(path2data + '*.h5')
    Nb_minutes = 1 # duration of each stack
    
    for idx_file in range(len(filelist)):
        file2load = filelist[idx_file]
        with h5py.File(file2load,'r') as f:
            print(list(f.keys()))
            
            a_group_key = list(f.keys())[0]
            data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])
            
        # shape strain rate : 
        # dim 0 : time in second / dim 1 : time sampled at fs / dim 2 : space 
        strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
        t = data['Source1']['time'] #time since epoch
    
        # Stack strain rate into sections of duration Nb_minutes
        stack_strain,stack_time,s = new_stack_strain_rate(strain_rate, t, fiber_length,Nb_minutes)
        
        stack_epoch = t[0] + stack_time
        UTC_stack = np.empty(stack_time.shape,dtype = 'object')
        for i in range(stack_epoch.shape[0]):
            current_UTC = epoch2datetime(stack_epoch[i,:],timezone = UTC_timearea)
            UTC_stack[i,:] = current_UTC
        
        # Show stacks 
        UTC_0 = UTC_stack[0,0]
        txt_UTC = generate_datetime_txt(UTC_0)
        subfig_folder = f'{fig_folder}idx_file_{idx_file}_UTC_{txt_UTC}/'
        if not os.path.isdir(subfig_folder):
            os.mkdir(subfig_folder)
            
        # loop over all stacks
        for idx in range(stack_strain.shape[0]):
            
            spatio = stack_strain[idx,:,:]
            current_t = stack_time[idx,:]
            current_UTC = UTC_stack[idx,:]
            
            fig,ax,imsh,cbar = plot_spatio_temp(spatio, current_UTC, s, fiber_length)
            imsh.set(clim = [1, 0.6e4])
            ax.set_title(f'{current_UTC[0]} - {current_UTC[-1]}')
            ax.set_xlabel(r'UTC',labelpad = 5)
            
            initial_UTC = generate_datetime_txt(current_UTC[0])
            figname = f'{subfig_folder}spatiotemporal_duration_{Nb_minutes}_UTC_{initial_UTC}'
            
            plt.savefig(figname + '.pdf', bbox_inches='tight')
            plt.savefig(figname + '.png', bbox_inches='tight')
            
            plt.close('all')



if __name__ == '__main__':
    main()
