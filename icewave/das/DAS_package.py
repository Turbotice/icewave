# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:03:00 2025

@author: sebas

Package for DAS analysis. Gathers all functions useful to study DAS signals

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
import scipy

import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.matlab2python as mat2py

# PARULA COLORMAP 
parula_map = matcmaps.parula()


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
        spatio_long,s,_ = time_stacking(current_section,Nb_seconds,
                                        fiber_length) # stack for a given duration Nb_minutes
        stack_strain[n,:,:] = spatio_long
    
        # create array of time
        current_time = np.arange(n*Nb_seconds, (n + 1)*Nb_seconds, 1/fs)
        stack_time[n,:] = current_time 
        

    return stack_strain,stack_time,s 
    
#----------------------------------------------------------------------------------------------------

def stack_data_fromfile(file2load,fiber_length,Nb_minutes, timezone = pytz.timezone('UTC')):
    """ Collect and structure data stored in a .h5 file. Strain rate measurements are stacked into chunks of 
    duration Nb_minutes. The UTC time of recording is also computed for each chunk.
    
    Input: - file2load, string, path to .h5 file to be loaded. 
           - fiber_length, float, length of optical fiber in meter
           - Nb_minutes, int or float, duration in minutes of each chunk
           - timezone, pytz.timezone object, optional argument used to define in which timezone the
           array of time should be refered to. Default is UTC. 
    Output: - stack_strain, 3D numpy array, #dim 0 : each stack of duration Nb_minutes, #dim 1 : time sampled at fs within 
            the corresponding stack, #dim 2 : space
            - stack_time, 2D numpy array # dim 0: each stack of duration Nb_minutes, # dim 1: time in seconds, taking 
            beginning of recording as a reference
            - UTC_stack, 2D numpy array, # dim 0 : each stack / chunk of duration Nb_minutes, # dim 1 : datetime object, 
            UTC based
            - s, numpy array, curvilinear coordinate """
            
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
        current_UTC = epoch2datetime(stack_epoch[i,:],timezone = timezone)
        UTC_stack[i,:] = current_UTC

    return stack_strain,stack_time,UTC_stack,s

#-----------------------------------------------------------------------------------------------------------------

def time_decimation_stack_strain(stack_strain,stack_time,UTC_stack,down_sampling_factor,axis = 1):
    """ Decimate arrays of strain rate, experimental time and UTC datetime (returned by function new_stack_strain_rate)
    Inputs : - stack_strain, 3D numpy array, #dim 0 : each stack of duration Nb_minutes, #dim 1 : time sampled at fs within 
            the corresponding stack, #dim 2 : space
            - stack_time, 2D numpy array # dim 0: each stack of duration Nb_minutes, # dim 1: time in seconds, taking 
            beginning of recording as a reference
            - UTC_stack, 2D numpy array, # dim 0 : each stack / chunk of duration Nb_minutes, # dim 1 : datetime object, 
            UTC based 
            - down_sampling_factor, float or int, down sampling factor 
            - axis, axis over which decimation is performed, default is 1
    Outputs: decimated versions of stakc_strain, stack_time and UTC_stack """
    
    # convert datetime objects into time since epoch
    epoch_matrix = np.zeros(UTC_stack.shape)
    for i in range(UTC_stack.shape[0]):
        t_epoch = [UTC_datetime.timestamp() for UTC_datetime in UTC_stack[i,:]]
        epoch_matrix[i,:] = t_epoch
        
    # decimate
    decimated_epoch = epoch_matrix[:,::down_sampling_factor]
    decimated_time = stack_time[:,::down_sampling_factor]
    decimated_stack_strain = scipy.signal.decimate(stack_strain,down_sampling_factor,axis = 1) # anti-aliasing filter

    # convert back to datetime objects
    UTC_decimated = np.zeros(decimated_epoch.shape,dtype = object)
    for i in range(UTC_decimated.shape[0]):
        UTC_decimated[i,:] = epoch2datetime(decimated_epoch[i,:],timezone = pytz.timezone('UTC'))
    
    print(f'Decimation performed with a factor {down_sampling_factor:.2f}')

    return decimated_stack_strain, decimated_time, UTC_decimated

#--------------------------------------------------------------------------------------------------------------


def get_DAS_parameters(path2DAS_param,date):
    """ Collect DAS acquisition parameters 
    Inputs : - path2DAS_param, string, path to table of DAS acquisition parameters
             - date, string, date for wich parameters are collected, format 'mmdd'
    Outputs : - fs, float, sampling frequency
              - fiber_length, float, fiber length set in the interrogator (in meters)
              - facq_x, float, spatial acquisition frequency (points/meter) """
    
    with open(path2DAS_param,'rb') as pf:
        param = pickle.load(pf)
    print('Parameters file loaded')
    
    # Set parameters
    # fs = np.shape(strain_rate)[1] # time sampling frequency 
    # facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
    fs = param[date]['fs']
    fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
    facq_x = param[date]['facq_x'] 
    
    return fs,fiber_length,facq_x