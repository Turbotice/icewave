# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:41:26 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
import h5py 
import glob
from time import strftime, localtime
from datetime import datetime
import pytz
import time 


import icewave.tools.matlab2python as mat2py


plt.rcParams.update({
    "text.usetex": True}) # use latex
#%%
def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]
#%%

date = '0212'
path2data = f'E:/Data/{date}/DAS_h5/'

filelist = glob.glob(path2data + '*.h5')
i = 1
file2load = filelist[i]
with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
    
    a_group_key = list(f.keys())[0]
    data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])
    

#%% 

# shape strain rate : 
# dim 1 : time in second / dim 2 : time sampled at fs / dim 3 : space 
strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
t = data['Source1']['time'] #time since epoch

#%% Convert time to UTC time

fs = 840 # time sampling frequency 
fx = 1 # spacial sampling frequency
fiber_length = np.shape(strain_rate)[2]

t_1sec = np.arange(0,1,1/fs)
s = np.arange(0,fiber_length,1/fx)

format_date = '%Y-%m-%d %H:%M:%S.f'
local_timearea = pytz.timezone('America/Montreal')
UTC_timearea = pytz.timezone('UTC')

# Convert time since epoch to UTC time 
UTC_t = []
local_t = []
for i,t_epoch in enumerate(t):
    current_string = strftime(format_date, localtime(t_epoch))
    current_time = datetime.strptime(current_string, format_date)
    local_t.append(current_time)
    UTC_time = current_time.astimezone(UTC_timearea)
    UTC_t.append(UTC_time)
    
#%% Show spatio-temp for a given time 

normalization = 'linear'
fig, ax = plt.subplots()
pause = 3

for i0 in range(np.shape(strain_rate)[0]):

    spatio = strain_rate[i0,:,:]
    imsh = ax.imshow(spatio.T,origin = 'lower', aspect = 'auto', norm = normalization,
              extent = extents(t_1sec) + extents(s))
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$s \; \mathrm{(m)}$')
    
    ax.set_title(f'local time : {local_t[i0]}')
    plt.pause(1)
    ax.clear()

    # cbar = plt.colorbar(imsh)

#%% Try to concatenate several seconds 


    
    
    
    