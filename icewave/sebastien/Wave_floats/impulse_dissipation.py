# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:58:17 2025

@author: sebas
"""

import numpy as np
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
import pickle
import scipy.signal
from scipy.signal import find_peaks


# import modules 
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs

parula_map = matcmaps.parula()

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/Referee_answer/'

#%% FUNCTION SECTION 

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

#------------------------------------------------------------------------------------------------------------------------

def indices2fit(y,x,p,rel_height):
    """ Compute indices of points to be fitted, based on index of a peak p and the y-values to fit
    Inputs : y, numpy.ndarray, values to fit
             x, numpy.ndarray, abscisse coordinates
             p, index of the detected peak, y[p] is therefore a local maximum
             rel_height, float between 0 and 1, relative prominence at which peak width should be computed if left or right limits
             are too large
    Outputs : points2fit, indices of point to be fitted"""

    # compute left/right limits
    decreasing_curve = 1
    left_limit = p
    right_limit = p
    while decreasing_curve :
        next_left = left_limit-1
        next_right = right_limit+1
        if np.logical_or(y[next_left] >= y[left_limit],
                         y[next_right] >= y[right_limit]) : 
            decreasing_curve = 0
        else : 
            left_limit = next_left
            right_limit = next_right

        if np.logical_or(left_limit == 0,right_limit == len(x)-1):
            decreasing_curve = 0
            widths = scipy.signal.peak_widths(y,[p],rel_height = rel_height)
            left_limit = math.floor(widths[2][0])-1
            right_limit = math.ceil(widths[3][0])+1
            
    # compute left/right distance to peak
    left_ips = x[p] - x[left_limit]
    right_ips = x[right_limit] - x[p] 
    # min_ips = min(left_ips,right_ips)
    points2fit = np.where(np.logical_and(x > x[p] - left_ips,x < x[p] + right_ips))[0]
    # points2fit = np.where(np.logical_and(Afk['k'] > Afk['k'][p] - min_ips,Afk['k'] < Afk['k'][p] + min_ips))[0]
    
    return points2fit



#%%  Load data 

base = 'F:/Waves_reconstruction_wilson/8_hammering_2024_01_02/'
exp_ID = 'D6cm_h15mm'
test_ID = 'test2'

path2data = f'{base}{exp_ID}/center/{test_ID}/'
filelist = glob.glob(f'{path2data}*Height_in_cm_surf.mat')
file2load = filelist[0]

data = {}
with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    for i in range(1,len(list_keys) - 1):
        key = list_keys[i]
        data[key] = mat2py.mat_to_dict(fmat[key],fmat[key])
        
    data['table'] = mat2py.mat_to_dict(fmat['table'],fmat['#refs#'])
    
    data['H_ww'] = np.flip(data['H_ww'],axis = 2)

#%% Show a frame 

frame = 250
field = data['H_ww'][frame,:,:]
extents = np.array([data['X'][0] , data['X'][-1], data['Y'][0], data['Y'][-1]])/data['fx']*1e2

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(field.T,cmap = parula_map, aspect = 'equal', 
                 origin = 'lower', extent = extents, vmin = -0.05, vmax = 0.05) 

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$')


ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)


#%% Compute FFT in time averaged over space 

fps = 200 # frame rate 
idx_start = 140
idx_end = 350
TF_spectrum,freq,FFT_t = FT.temporal_FFT(np.transpose(data['H_ww'],(1,2,0))[:,:,idx_start:idx_end],fps,padding_bool = 1,
                                         add_pow2 = 1,output_FFT = True)

print('Temporal FFT computed !')
#%% Plot spectrum 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([3e-5,1e-2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')


#%% Find peak and gets its width at half maximum

rel_height = 0.5
delta_f = freq[1] - freq[0]
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
prominences = scipy.signal.peak_prominences(TF_spectrum,peaks)
widths = scipy.signal.peak_widths(TF_spectrum,peaks,rel_height = rel_height)


fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')
ax.plot(freq[peaks],TF_spectrum[peaks],'rd')

ax.vlines(x=freq[peaks], ymin=TF_spectrum[peaks] - prominences[0],
        ymax = TF_spectrum[peaks],color = 'g')
ax.hlines(y=widths[1], xmin=widths[2]*delta_f,
        xmax=widths[3]*delta_f,color = 'g')

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim([3e-5,1e-2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

alpha = widths[0]*2*np.pi*delta_f
print(alpha)

#%% Find peak and fit it by a lorentzian 

# find peaks 
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
first_peak = peaks[0] # get first peak 

# normalize spectrum
norm_spectrum =(TF_spectrum - TF_spectrum.min())/(TF_spectrum.max() - TF_spectrum.min())

# find points to fit 
points2fit = indices2fit(norm_spectrum,freq,first_peak,0.85)

# fit by a lorentzian
x = freq[points2fit]*2*np.pi # pulsation
y_exp = norm_spectrum[points2fit]
popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                     bounds = ([x[0],1e-4],[x[-1],10]))
err_coeff = np.sqrt(np.diag(pcov))
print(f'Omega = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')

xfit = np.linspace(x[0],x[-1]) 
yth = lorentzian(xfit, popt[0], popt[1])
scaled_theory = yth*(TF_spectrum.max() - TF_spectrum.min()) + TF_spectrum.min()
label_th = r'$\frac{1}{\sqrt{1 + (\omega - \Omega)^2/\lambda^2}}$'

set_graphs.set_matplotlib_param('single')
fix, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')
ax.plot(freq[points2fit],TF_spectrum[points2fit],'ro')
ax.plot(xfit/2/np.pi,scaled_theory,'r',label = label_th)

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim([-1e-3,1e-2])
ax.set_xlim([-1,40])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

ax.legend(loc = 'upper right')

figname = f'TF_spectrum_impulse_lorentzian_linear_{exp_ID}_{test_ID}_Omega_{popt[0]:.2f}_lambda_{popt[1]:.2f}'
figname = figname.replace('.','p')
figname = f'{fig_folder}{figname}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')