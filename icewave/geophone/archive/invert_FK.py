#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:23:22 2024

@author: moreaul
"""


from scipy import fft as fft
import numpy as np
import matplotlib.pyplot as plt
import pickle

def wavenumbers_stein( rho_ice, h, E, nu,freq,c_w,rho_w):
    
    g = 9.81
    G = E/(2*(1+nu))
    cS0 = np.sqrt(E/(rho_ice*(1-nu**2)))
    cSH0 = np.sqrt(G/rho_ice)
    D = E*pow(h,3)/(12*(1-nu**2))

    k = np.linspace(1e-12,5,100000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for kf in range(len(freq)):
        omeg = 2*np.pi*freq[kf]
        if omeg == 0:
            flag = 1;
            idx_flag = kf;
        else:
            cph = omeg/k
            
            func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            func[func.imag != 0] = -1
            func = func.real
            idx_zero[kf] = (np.where(np.diff(np.signbit(func)))[0])
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero]       
    if flag:
        k_QS[idx_flag] = 0
        
    k_QS0 = freq/cS0*2*np.pi
    k_SH0 = freq/cSH0*2*np.pi  
    cphQS = freq/k_QS*2*np.pi

    
    return k_QS, k_QS0, k_SH0, cphQS

#%%
rho_ice = 917
E = 7.31e9# 2.43e9
nu = 0.4
c_w = 1450 # waves celerity in water 
rho_w = 1027 # density of water 
h = np.arange(0.1,0.6,0.005)
path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones'
acqu_numb = '0001'

file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir1.pkl'
with open(file2load, "rb") as f:
    data = pickle.load(f)
f_mode1 = data[0] 
k_mode1 = data[1] 
plt.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points')    
    
f_mode1 = f_mode1[1:]
k_mode1 = k_mode1[1:]

# file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir2.pkl'
# with open(file2load, "rb") as f:
#     data = pickle.load(f)
# f_mode2 = data[0] 
# k_mode2 = data[1] 


plt.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points') 
#plt.plot(f_mode2, k_mode2, linestyle='--', color='r', label='Line between points') 

f_mode = f_mode1
k_mode = k_mode1


l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h[i], E, nu,f_mode,c_w,rho_w)
    error = np.sum(np.power((k_mode-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)]
   

k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h_ice, E, nu,f_mode,c_w,rho_w)

plt.plot(f_mode, k_mode_synthetic, color='g', label='Line between points')    
plt.plot(f_mode, k_mode, linestyle='--', color='r', label='Line between points') 

 
plt.plot(h,l2_norm)

print(h_ice)



