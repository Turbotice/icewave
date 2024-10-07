# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 09:55:40 2024

@author: sebas

This script aims at understanding difference between different models of flexural waves 

"""

import numpy as np 
import matplotlib.pyplot as plt
#%% ----------------------------- FUNCTION SECTION ---------------------

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

#%%

# Parameters
rho_ice = 917  # Density of ice in kg/m³
H = 100  # Depth of water in meters
nu = 0.3  # Poisson's ratio
c_w = 1480  # Speed of sound in water in m/s
rho_w = 1000  # Density of water in kg/m³

E_values = 4.0e9
h_values = 3.0

Nb_pts = 100
freq = np.linspace(2,300,Nb_pts)

wavenumbers = np.zeros((2,np.size(freq)))

wavenumbers[0,:] = wavenumbers_stein_squire(rho_ice, h_values, H, E_values, nu, freq, c_w, rho_w, 'stein')
wavenumbers[1,:] = wavenumbers_stein_squire(rho_ice, h_values, H, E_values, nu, freq, c_w, rho_w, 'squire')

fig, ax = plt.subplots(figsize = (10,10))
ax.plot(wavenumbers[0,:],freq,'ro',label = 'Stein')
ax.plot(wavenumbers[1,:],freq,'bo',label = 'Squire')
ax.grid()
ax.set_xlabel(r'$k \: \rm (rad.m^{-1})$')
ax.set_ylabel(r'$f \: \rm (Hz)$')
ax.legend()