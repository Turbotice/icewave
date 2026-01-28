# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 10:56:17 2026

@author: sebas

This script enables the quick computation of sea ice mechanical properties from a line of geophones : 
    - Young's modulus E
    - Poisson coefficient nu
    - ice thickness h

First velocities of QS0 modes (longitudinal/acoustic) and SH0 modes (shear) must have been computed 
We also need to extract coordinates (f_mode,k_mode) of flexural waves dispersion relation. 

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pickle 
import scipy

import icewave.geophone.package_geophone as geopack
#%% Import data 

year = '2025'
date = '0210' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 

path2data = os.path.join('U:/Data/',date,'Geophones/')
signal_length = 1

# load c_QS0 and c_SH0
file2load = f'{path2data}Phase_velocity_dictionnary_acqu_{acqu_numb}_sig_length_' + str(signal_length).replace('.','p') + '.pkl'
with open(file2load,'rb') as pfile:
    s = pickle.load(pfile)

# load (f,k) flexural_mode
for direction in ([1,2]):
    key_dir = f'dir{str(direction)}'
    file2load =  f'{path2data}{year}_{date}_acq{acqu_numb}disp_QS_dir{str(direction)}.pkl'
    with open(file2load,'rb') as pfile:
        f_mode,k_mode = pickle.load(pfile)
        s[key_dir]['f'] = f_mode
        s[key_dir]['k'] = k_mode

#%% Combine velocities from different directions 

results = {}
results['C_shear'] = 0.5*(s['dir1']['C_shear'] + s['dir2']['C_shear'])
results['C_longi'] = 0.5*(s['dir1']['C_longi'] + s['dir2']['C_longi'])
results['uC_shear'] = 0.5*(s['dir1']['uC_shear'] + s['dir2']['uC_shear'])
results['uC_longi'] = 0.5*(s['dir1']['uC_longi'] + s['dir2']['uC_longi'])


###################################################
#%% Computes Young modulus & Poisson coefficient 
###################################################

results['rho_ice'] = 917


def compute_nu(C_shear,C_longi):
    nu = 1-2*(C_shear/C_longi)**2
    return nu

def compute_E(C_longi,nu,rho_ice):
    E = rho_ice*C_longi**2*(1-nu**2)
    return E

def compute_unu(C_shear,C_longi,uC_shear,uC_longi):
    unu = 4*C_shear/(C_longi**2) * uC_shear + 4*C_shear**2/(C_longi**3) * uC_longi
    return unu

def compute_uE(C_longi,nu,rho_ice,uC_longi,unu):
    uE = 2*rho_ice * (C_longi*(1-nu**2)*uC_longi + (C_longi**2)*nu*unu)
    return uE

results['nu'] = compute_nu(results['C_shear'],results['C_longi'])
results['E'] = compute_E(results['C_longi'],results['nu'],results['rho_ice'])
results['unu'] = compute_unu(results['C_shear'],results['C_longi'],results['uC_shear'],results['uC_longi'])
results['uE'] = compute_uE(results['C_longi'],results['nu'],results['rho_ice'],
                           results['uC_longi'],results['unu'])

###################################
#%% --- Invert ice thickness -- 
###################################

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 

# direction = 1
# key_dir = f'dir{str(direction)}'
# f_mode = s[key_dir]['f']
# k_mode = s[key_dir]['k']

# interpolate (f,k) points from dispersion relation in both directions
mini = max(min(s['dir1']['f']),min(s['dir2']['f']))
maxi = min(max(s['dir1']['f']),max(s['dir2']['f']))

freq = np.linspace(mini,maxi , 15)

# Interpolate kQS1 at the new frequency points, without extrapolation
interpolator1 = scipy.interpolate.interp1d(s['dir1']['f'], s['dir1']['k'])
kQS1_interp = interpolator1(freq)

# Interpolate kQS2 at the new frequency points, without extrapolation
interpolator2 = scipy.interpolate.interp1d(s['dir2']['f'], s['dir2']['k'])
kQS2_interp = interpolator2(freq)

# Calculate the average values, ignoring NaNs
kQS = np.nanmean([kQS1_interp, kQS2_interp], axis=0)

# Plot the results
fig, ax = plt.subplots()
ax.plot(s['dir1']['k'], s['dir1']['f'],  'bo', label='kQS1')
ax.plot(s['dir2']['k'], s['dir2']['f'],  'go', label='kQS2')
ax.plot(kQS1_interp, freq,  'b--', label='kQS1 interpolated')
ax.plot(kQS2_interp, freq, 'g--', label='kQS2 interpolated')
ax.plot(kQS, freq, 'r--', label='Average kQS')
ax.set_xlabel(r'$k_{QS} \: \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \: \mathrm{rad.m^{(-1})}$')
ax.legend()

#%% Find ice thickness

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h[i], 
                                                              results['E'], results['nu'],freq,c_w,rho_w)
    error = np.nansum(np.power((kQS-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
print(h_ice)

#%% computes wavevectors k 
k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h_ice, 
                                                          results['E'], results['nu'],freq,c_w,rho_w)

fig, ax = plt.subplots()
ax.plot(k_mode_synthetic,freq,color = 'g')
ax.plot(kQS,freq,linestyle = '--', color = 'r')

results['h_ice'] = h_ice

#%%

fig, ax = plt.subplots()
ax.plot(h,l2_norm)


    
