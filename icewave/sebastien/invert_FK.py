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
import csv 

import seb
from seb.pickle_m import read, write

def wavenumbers_stein( rho_ice, h, E, nu,freq,c_w,rho_w):
    """ This function computes the wave vectors associated to a given array of frequencies
    It takes as arguments : 
        - rho_ice : ice density 
        - h : a given thickness of ice 
        - E : Young modulus of ice
        - nu : Poisson coefficient of ice
        - freq : an array of frequencies, to which will correspond wave vectors 
        - c_w : waves phase velocity
        - rho_w : water density 
        
    The function returns : 
        - k_QS : wave vectors of the flexural mode
        - k_QS0 : wave vecotrs of the acoustic mode
        - k_SH0 : wave vectors of the shear mode
        - cphQS : phase velocity of the flexural mode"""
    
    
    g = 9.81
    G = E/(2*(1+nu))
    cS0 = np.sqrt(E/(rho_ice*(1-nu**2))) # celerity of longitudinal wave
    cSH0 = np.sqrt(G/rho_ice) # celerity of shear wave 
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    k = np.linspace(1e-12,10,100000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for kf in range(len(freq)):
        omeg = 2*np.pi*freq[kf]
        if omeg == 0:
            flag = 1;
            idx_flag = kf;
        else:
            cph = omeg/k # phase velocity
            # Ludovic version
            func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            # Sebastien version (Stein 1998)
            # func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_ice/D + pow(omeg/cph,4)
            
            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[kf] = (np.where(np.diff(np.signbit(func)))[0]) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero]       
    if flag:
        k_QS[idx_flag] = 0
        
    k_QS0 = freq/cS0*2*np.pi # wave vector associated to longitudinal wave
    k_SH0 = freq/cSH0*2*np.pi   # wave vector associated to shear wave
    cphQS = freq/k_QS*2*np.pi # phase velocity of the flexural wave

    
    return k_QS, k_QS0, k_SH0, cphQS


def hydro_elastic(f,h,rho_ice,E,nu):
    """ Computes the hydro-elastic wave vector associated to an array frequency. 
    The function takes the following arguments : 
        - f : frequency arrray 
        - h = thickness of ice 
        - rho_ice = density of ice
        - E = Young modulus of ice
        - nu = Poisson coefficient 
    """
    
    omega = 2*np.pi*f
    k = omega**2 * (12*rho_ice*(1 - nu**2))/(E*h**3)
    k = np.power(k,1/5)
    
    return k

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#%%
rho_ice = 917
E = 2.6062e9# 2.43e9
nu = 0.326
c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 
h = np.arange(0.1,0.6,0.005) # array of height tested 
#path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones' # path 2 data points saved from FK flexural
acqu_numb = '0001' # acquisition number 

file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir1.pkl'
with open(file2load, "rb") as filename:
    data = pickle.load(filename)
f_mode1 = data[1] 
k_mode1 = data[0] 
plt.plot(f_mode1, k_mode1, linestyle='--', color='r', label='Line between points')    
    
f_mode1 = f_mode1[1:]
k_mode1 = k_mode1[1:]

# file2load = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir2.pkl'
# with open(file2load, "rb") as f:
#     data = pickle.load(f)
# f_mode2 = data[0] 
# k_mode2 = data[1] 


plt.plot(k_mode1, f_mode1, linestyle='--', color='r', label='Line between points') 
#plt.plot(f_mode2, k_mode2, linestyle='--', color='r', label='Line between points') 

f_mode = f_mode1
k_mode = k_mode1


l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h[i], E, nu,f_mode,c_w,rho_w)
    error = np.sum(np.power((k_mode-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
   

k_mode_synthetic, k_QS0, k_SH0, cphQS = wavenumbers_stein( rho_ice, h_ice, E, nu,f_mode,c_w,rho_w)

plt.plot(k_mode_synthetic, f_mode, color='g', label='Line between points')    
plt.plot(k_mode, f_mode, linestyle='--', color='r', label='Line between points') 

 
plt.plot(h,l2_norm)

print(h_ice)

#%% Save all relevant data into a file pickle
pickle_file = path2data + '/' + acqu_numb + '/' + 'Physical_data_direction_1.pkl'

s = dict()
s['C_shear']= C_shear
s['C_longi'] = C_longi
s['E'] = E
s['nu'] = nu
s['h_ice'] = h_ice
s['rho_ice'] = rho_ice
s['c_w'] = c_w
s['rho_w'] = rho_w

write(s,pickle_file)

#%% Save relevant data into a csv file 

csv_file = path2data + '/' + acqu_numb  + '/' + 'Physical_data_direction_1_fk.csv'

with open(csv_file, 'w') as csvfile: 
    writer = csv.DictWriter(csvfile, fieldnames = s.keys()) 
    writer.writeheader() 
    for key in s.keys(): 
        csvfile.write("%.3f"% s[key] + ',')
        
#%% Create an other plot with data points and fitted curve 

frequency_list = np.linspace(1,150,100)
k_flex_theory, k_acoust_theory, k_shear_theory, cphQS = wavenumbers_stein( rho_ice, h_ice, E, nu,frequency_list,c_w,rho_w)

fig, ax = plt.subplots()
FK = FK/np.max(FK)
c1 = ax.contourf(F, K, FK, cmap='gnuplot2', vmin=vmin, vmax=vmax)
ax.scatter(data[0],data[1], label = 'Detected points')
ax.plot(frequency_list,k_flex_theory,color = 'k', label = 'Flexural mode')
plt.colorbar(c1, ax=ax, label= r'Spectrum Magnitude')

ax.set_xlim([0, 100])
ax.grid()
ax.set_xlabel(r'$f \quad \mathrm{(Hz)}$')
ax.set_ylabel(r'$k \quad \mathrm{(m^{-1})}$')
# ax.scatter(f_mode,k_QS0)
# ax.scatter(f_mode,k_SH0)


#%% Superposition of the Unwrapped FK spectrum and the fitted flexural law 


frequency_list = np.linspace(1,200,100)
k_flex_theory, k_acoust_theory, k_shear_theory, cphQS = wavenumbers_stein( rho_ice, h_ice, E, nu,frequency_list,c_w,rho_w)
k_hydro_elast = hydro_elastic(frequency_list, h_ice, rho_ice, E, nu)
fig, ax1 = plt.subplots(1, 1, figsize=(18, 9))
vmin = 0 
vmax = 1 

c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)

#ax1.set_ylim([-1.5, 1.5])
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$', labelpad = 5)
# ax1.set_title(r'Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label = r'$\frac{|\hat{s}|}{|\hat{s}_{max}|}(f,k)$')
ax1.plot(k_flex_theory,frequency_list,linestyle = '--', 
         linewidth = 5 ,color = 'white', label = 'Flexural mode')
ax1.tick_params(axis='both', which='major', pad=7)
# ax1.plot(k_hydro_elast,frequency_list,linestyle = '--', 
#          linewidth = 3, color = 'r')
ax1.set_ylim([0 , 200])

#%%
figname = fig_folder + 'Unwrapped_flexural_with_theory'

plt.savefig(figname + '.pdf', dpi = 1200, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = 1200, bbox_inches = 'tight')


