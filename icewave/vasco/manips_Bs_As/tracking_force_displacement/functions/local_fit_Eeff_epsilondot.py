#%%
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

from load_force_vs_frame import *
from load_matdata_JZmethod import load_displacement_px
from load_displacement_vs_frame import load_disp_vs_frame

from force_disp import load_force_displacement_curve

#%% fonctions moyennes glissantes

def moving_nanmean(arr, window):
    arr = arr.astype(float)
    half = window // 2
    
    result = np.array([
        np.nanmean(arr[max(0, i-half):min(len(arr), i+half+1)])
        for i in range(len(arr))
    ])
    
    return result

def moving_nan_gaussian(arr, sigma, window=None):
    arr = arr.astype(float)
    
    if window is None:
        window = int(6 * sigma + 1)  # standard
    
    half = window // 2
    
    # noyau gaussien
    x = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    
    result = np.empty(len(arr))
    
    for i in range(len(arr)):
        start = max(0, i - half)
        end = min(len(arr), i + half + 1)
        
        slice_arr = arr[start:end]
        
        # adapter le noyau à la taille (bords)
        k = kernel[(start - (i - half)):(end - (i - half))]
        
        mask = ~np.isnan(slice_arr)
        
        if np.any(mask):
            weights = k[mask]
            values = slice_arr[mask]
            result[i] = np.sum(values * weights) / np.sum(weights)
        else:
            result[i] = np.nan
    
    return result

def gradient_k(arr, k):
    arr = arr.astype(float)
    
    grad = np.full_like(arr, np.nan)
    
    grad[k:-k] = (arr[2*k:] - arr[:-2*k]) / (2*k)
    
    return grad
#%% load dict results (for one sample)
#for idx in range(10):

idx = 12

dict_results, _, _, _, _, _, _ = load_force_displacement_curve(idx=idx,savepkl=False,plot=False)

force = dict_results['force_newtons_comidx']
disp = dict_results['displacement_meters_comidx']
idx_frac_common = dict_results['idx_frac_common']
L = dict_results['L']
w = dict_results['w']
h_avg = dict_results['h_avg']
h_std = dict_results['h_std']
facq = dict_results['freq_acq_Hz']

mask = np.arange(len(force)) <= idx_frac_common

force_withnans = np.where(mask, force, np.nan)
disp_withnans = np.where(mask, disp, np.nan)

#disp_smoothed = moving_nanmean(disp_withnans, window=10)
disp_smoothed = moving_nan_gaussian(disp_withnans, sigma=3)


#flexural_rigidity = np.gradient(force_withnans)/np.gradient(disp_smoothed)
flexural_rigidity = gradient_k(force_withnans,k=5)/gradient_k(disp_smoothed,k=5)

#flexural_rigidity_smoothed = moving_nan_gaussian(flexural_rigidity, sigma=3)

plt.figure()
plt.plot(force_withnans, disp_withnans)
plt.plot(force_withnans, disp_smoothed)
plt.show()

plt.figure()
plt.plot(force_withnans, flexural_rigidity)
#plt.plot(force_withnans, flexural_rigidity_smoothed)
plt.ylabel('Rigidité en flexion (Fz/deplacement)')
plt.show()


plt.figure()
plt.plot(np.gradient(disp_smoothed) * facq, flexural_rigidity)
#plt.plot(force_withnans, flexural_rigidity_smoothed)
plt.ylabel('Rigidité en flexion (Fz/deplacement)')
plt.xlabel('vz [m/s]')

#%% test pour extraire seulement quelques points de la figure

sys.path.append('../../../tools')
import clickonfigures
%matplotlib qt
clickonfigures.get_n_points(force_withnans, flexural_rigidity,n_points=4)




# %%
#plt.plot(np.gradient(disp_smoothed) * facq * 6 * h_avg/(L**2), (1/(4*w))*(L/h_avg)**3 * flexural_rigidity,'.', alpha=0.5)
#plt.xlim(1e-4,2e-3)
#plt.ylim(1e9,10e9)
#plt.loglog()
#plt.show()
