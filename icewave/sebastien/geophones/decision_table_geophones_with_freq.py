#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:05:16 2024

@author: moreaul
"""

import numpy as np
import os
import csv

#%%

def wavenumbers_stein_squire(rho_ice, h, H, E, nu, freq, c_w, rho_w, equation):
    g = 9.81
    D = E * pow(h, 3) / (12 * (1 - nu ** 2))

    k = np.linspace(1e-12, 5, 200000)

    idx_zero = np.zeros(len(freq))
    flag = 0
    for kf in range(len(freq)):
        omeg = 2 * np.pi * freq[kf]
        if omeg == 0:
            flag = 1
            idx_flag = kf
        else:
            cph = omeg / k
            if equation == 'stein':
                func = rho_w / D * (g - omeg / np.lib.scimath.sqrt((1 / cph) ** 2 - (1 / c_w) ** 2)) - h * omeg ** 2 * rho_ice / D + pow(omeg / cph, 4)
            elif equation == 'squire':
                coth = 1 / np.tanh(k * H)
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + coth) - D * pow(k, 5) / rho_w - g * k
            else:
                print('inappropriate equation name: choose between stein or squire')

            func[func.imag != 0] = -1
            func = func.real
            zero_crossing = np.where(np.diff(np.signbit(func)))[0]
            if zero_crossing.size > 0:
                idx_zero[kf] = zero_crossing[0]
            else:
                idx_zero[kf] = 0

    idx_zero = idx_zero.astype(int)
    k_QS = k[idx_zero]
    if flag:
        k_QS[idx_flag] = 0

    return k_QS


def QS0_nondispersive(rho_ice,E,nu,freq):
    "Compute wavenumber according to non-dispersive QS0 mode"
    
    omeg = 2 * np.pi*freq
    kQS0 = omeg * np.sqrt(rho_ice * (1 - nu**2)/E)
    return kQS0

#%%
    
    
    
# Parameters
rho_ice = 917  # Density of ice in kg/m³
H = 100  # Depth of water in meters
nu = 0.3  # Poisson's ratio
c_w = 1480  # Speed of sound in water in m/s
rho_w = 1000  # Density of water in kg/m³

freq = [5]  # Frequencies in Hz

# E from 2 to 9 GPa by step of 1 GPa (convert to Pa)
E_values = np.arange(2e9, 12e9, 2e9)

# Ice thickness h from 0.5 to 3 meters by step of 0.5 meters
h_values = np.arange(0.5, 4.5 , 0.5)

# Calculate wavenumbers for each combination of E and h
wavenumbers = np.zeros((np.size(h_values),np.size(E_values),np.size(freq)))
wavelength = wavenumbers

for i,h in enumerate(h_values):
    for j,E in enumerate(E_values):
        k_values = wavenumbers_stein_squire(rho_ice, h, H, E, nu, freq, c_w, rho_w, 'stein')
        lambda_values = 2 * np.pi / k_values
        wavenumbers[i,j,:] = k_values
        wavelength[i,j,:] = lambda_values

# Print results
# for (E, h), lambda_values in wavelength.items():
#     print(f"E = {E / 1e9:.1f} GPa, h = {h:.1f} m -> wavelength (m): {lambda_values}")


#%% Write data in a csv file. A single csv for each frequency 

path = 'D:/PC Seb/These PMMH/Arctic_refuge_2024/Expedition_plan/Geophones/Decision_tables/'
for k,freq_value in enumerate(freq) :
    csv_file = path + 'Decision_table_geophones_f'+ str(freq_value) +'.csv'
    
    if not os.path.isdir(path):
        os.mkdir(path)
    
    with open(csv_file, 'w', newline = '') as csvfile: 
        writer = csv.writer(csvfile) 
    
        row = ['h_ice\E']  + list(np.round(E_values*1e-9,2))
        writer.writerow(row)
    
        for i, h_ice_value in enumerate(h_values):
            row = [round(h_ice_value,3)] + list(np.round(wavelength[i,:,k],1))
            writer.writerow(row)

#%% Computation of minimal wavelength of validity of Stein equation (flexural waves)

limit_stein = 50 # validity limit of Stein equation (in Hz.m)
E_values = np.arange(2e9, 12e9, 2e9)
h_values = np.arange(0.5, 4.5 , 0.5)

wavelength = np.zeros((np.size(h_values),np.size(E_values)))

for i,h in enumerate(h_values):
    for j,E in enumerate (E_values):
        freq = [limit_stein/h]
        print(freq)
        k_value = wavenumbers_stein_squire(rho_ice, h, H, E, nu, freq, c_w, rho_w, 'stein')
        wavelength[i,j] = 2 * np.pi / k_value
        
#%% Write data in a csv file 
path = 'D:/PC Seb/These PMMH/Arctic_refuge_2024/Expedition_plan/Geophones/Decision_tables/'

csv_file = path + 'Decision_table_geophones_limit_stein.csv'

if not os.path.isdir(path):
    os.mkdir(path)

with open(csv_file, 'w', newline = '') as csvfile: 
    writer = csv.writer(csvfile) 

    row = ['h_ice\E']  + list(np.round(E_values*1e-9,2))
    writer.writerow(row)

    for i, h_ice_value in enumerate(h_values):
        row = [round(h_ice_value,3)] + list(np.round(wavelength[i,:],1))
        writer.writerow(row)
        
        
#%% Computation of validity limit of QS0 mode (longitudinal waves)

limit_QS0 = 500 # maximal value thickness*frequency to get non-dispersive mode QS0

rho_ice = 917
nu = 0.3
E_values = np.arange(2e9, 12e9, 2e9)
h_values = np.arange(0.5, 4.5 , 0.5)

QS0_wavelength = np.zeros((np.size(h_values),np.size(E_values)))

for i,h in enumerate(h_values):
    for j,E in enumerate (E_values):
        freq = limit_QS0/h
        print(freq)
        k_value = QS0_nondispersive(rho_ice, E, nu, freq)
        QS0_wavelength[i,j] = 2 * np.pi / k_value
        

#%% Write data in a csv file 
path = 'D:/PC Seb/These PMMH/Arctic_refuge_2024/Expedition_plan/Geophones/Decision_tables/'

csv_file = path + 'Decision_table_geophones_QS0_minimal_wavelength.csv'

if not os.path.isdir(path):
    os.mkdir(path)

with open(csv_file, 'w', newline = '') as csvfile: 
    writer = csv.writer(csvfile) 

    row = ['h_ice\E']  + list(np.round(E_values*1e-9,2))
    writer.writerow(row)

    for i, h_ice_value in enumerate(h_values):
        row = [round(h_ice_value,3)] + list(np.round(QS0_wavelength[i,:],1))
        writer.writerow(row)