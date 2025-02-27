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
import scipy.signal as signal


import icewave.tools.matlab2python as mat2py


plt.rcParams.update({
    "text.usetex": True}) # use latex
#%%
def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

def time_stacking(strain_rate,nb_seconds,fiber_length):
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
    t_sec = np.arange(0,nb_seconds,1/fs)  
    
#%%

date = '0212'
# path2data = f'E:/Data/{date}/DAS_h5/'
path2data = 'F:/20250216_h5/'

filelist = glob.glob(path2data + '*.h5')
i = 0
file2load = filelist[i]
with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
    
    a_group_key = list(f.keys())[0]
    data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])
    
# shape strain rate : 
# dim 0 : time in second / dim 1 : time sampled at fs / dim 2 : space 
strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
t = data['Source1']['time'] #time since epoch

#%% Convert time to UTC time

fs = np.shape(strain_rate)[1] # time sampling frequency 
fiber_length = 700 # fiber length in meters (set on DAS)
fx = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency

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
    
    ax.set_ylim([450,700])
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$s \; \mathrm{(m)}$')
    
    ax.set_title(f'local time : {local_t[i0]}')
    plt.pause(3)
    ax.clear()

    # cbar = plt.colorbar(imsh)

#%% Try to concatenate several seconds 

i = 0
nb_seconds = 60 # total number of seconds to stack 

for i in range(nb_seconds):
    if i == 0 :    
        spatio_long = np.vstack((strain_rate[i,:,:],strain_rate[i + 1,:,:]))
    elif i > 1:
        spatio_long = np.vstack((spatio_long,strain_rate[i,:,:]))
        
print(f'Spatio-temp computed for a total number of seconds : {nb_seconds}')
t_sec = np.arange(0,nb_seconds,1/fs)    

# create an array of a given set of nb_seconds monitoring 



#%% 

normalization = 'log'
fig,ax = plt.subplots()
imsh = ax.imshow(spatio_long.T,origin = 'lower',aspect = 'auto',norm = normalization,
          extent = extents(t_sec) + extents(s),vmin = 1e-5, vmax = 1e5)
ax.set_ylim([450,700])
cbar = plt.colorbar(imsh)


#%% Build low pass filter 

wn = 4 # critical frequency of low pass filter 
order_filter = 10
b,a = signal.butter(order_filter,wn,'low',fs = fs)
w,h = signal.freqs(b,a)

fig,ax = plt.subplots()
ax.semilogx(w,20*np.log10(abs(h)))
ax.set_xlabel(r'$f \; \mathrm{(rad.s^{-1})}$')
ax.set_ylabel(r'Amplitude (dB)')
ax.grid(which = 'both',axis = 'both')

# Apply filter to signal 

filtered_spatio = signal.filtfilt(b,a,spatio_long,axis = 0)

fig,ax = plt.subplots()
imsh = ax.imshow(spatio_long.T,origin = 'lower',aspect = 'auto',norm = 'linear',
          extent = extents(t_sec) + extents(s),vmin = 1e-5, vmax = 1e5)
ax.set_ylim([450,700])
cbar = plt.colorbar(imsh)

fig,ax = plt.subplots()
imsh = ax.imshow(filtered_spatio.T,origin = 'lower',aspect = 'auto',norm = 'linear',
          extent = extents(t_sec) + extents(s),vmin = 1e-5, vmax = 1e5)
ax.set_ylim([450,700])
cbar = plt.colorbar(imsh)

#%% Look at a given position on the fiber 

s_check = 590 # position (in meter) at which we look the fiber signal 
idx = np.argmin(abs(s - s_check))

profile = spatio_long[:,idx]
profile_filtered = filtered_spatio[:,idx]


fig,ax = plt.subplots()
ax.plot(t_sec,profile)
ax.plot(t_sec,profile_filtered)