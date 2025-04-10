# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 13:36:48 2025

@author: sebas
"""

import numpy as np
import os
import glob 
import pandas as pd 
import h5py
import pickle
import cv2 as cv
import scipy
import math
# from skimage.filters import meijering, sato, frangi, hessian


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
from matplotlib.path import Path

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 
def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN


def unnormed_gaussian(x,x0,alpha):
    y = np.exp(-0.5*((x - x0)/alpha)**2)
    return y


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

#%% 
date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

path2data = f'K:/Share_hublot/Data/{date}/Drones/{drone_ID}/matData/{exp_ID}/'

file2load = glob.glob(f'{path2data}*scaled.mat')[0]
with h5py.File(file2load,'r') as fmat:
    print('Top-level keys : ', list(fmat.keys()))
    data = mat2py.mat_to_dict(fmat['m'],fmat['m'])

data = mat2py.transpose_PIVmat_fields(data)

#%% Compute Efk using python 
# reduce drone noise
Vs = FT.supress_quadratic_noise(np.transpose(data['Vx'],(1,0,2)),data['x'],data['y'])
Vs = np.transpose(Vs,(1,0,2))
Efk = FT.space_time_spectrum(Vs,1/data['SCALE']['fx'],data['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Define fig_folder

fig_folder = f'{path2data}Plots/attenuation_from_disp_relation/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Plot Afk 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(0.09,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)




#%% Try to detect branches using image processing 
Afk = Efk
# sobel filter
sobel_f = scipy.ndimage.sobel(Afk['E'],axis = 0)
sobel_k = scipy.ndimage.sobel(Afk['E'],axis = 1)
sobel_magnitude = np.hypot(sobel_k,sobel_f)

# laplacian
laplacian = scipy.ndimage.laplace(Afk['E'])

# gaussian filter 
sigma = 1
blurred_1 = scipy.ndimage.gaussian_filter(Afk['E'],sigma)
blurred_2 = scipy.ndimage.gaussian_filter(Afk['E'],sigma + 0.5)
dog = blurred_1 - blurred_2

fig, ax = plt.subplots()
c = ax.imshow(sobel_magnitude, cmap = 'gray', aspect = 'auto',norm = 'log',origin = 'lower',interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))
ax.set_xlim(0.1,3.3)
ax.set_ylim(0.09,1.3)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)
# c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
#               origin = 'lower', interpolation = 'gaussian',
#               extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))


figname = f'{fig_folder}Image_filtered_sobel_magnitude'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')




#%% Select a profile for a given wave vector


wavevectors = [0.65,0.85,1.0,1.5,1.9]
fig, ax = plt.subplots()
for select_k in wavevectors:
    idx = np.argmin(np.abs(select_k - Efk['k']))
    detected_k = Efk['k'][idx]
    print(detected_k)
    section = Efk['E'][:,idx]
    label_freq = f'{detected_k:.2f}'
    ax.plot(Efk['f'],section/section.max(),'-',label = label_freq)
    # detected_idx = np.argmin(np.abs(select_freq - Harm['f']))
    # detected_points = np.where(Harm['f'] == Harm['f'][detected_idx])[0]
    # print(Harm['k'][detected_points])
    # ax.plot(Harm['k'][detected_points],Harm['A'][detected_points],'s',mec = 'k')

# ax.set_yscale('log')
ax.set_xlim([1e-4,2])
ax.legend()




#%% Select a profile for a given frequency

frequencies = [0.25,0.35,0.45,0.6,0.7]
fig, ax = plt.subplots()
for select_freq in frequencies:
    idx = np.argmin(np.abs(select_freq - Efk['f']))
    detected_f = Efk['f'][idx]
    print(detected_f)
    section = Efk['E'][idx,:]
    label_freq = f'{detected_f:.2f}'
    ax.plot(Efk['k'],section,'-o',label = label_freq)
    ax.plot(Efk['k'],blurred_1[idx,:])
    # detected_idx = np.argmin(np.abs(select_freq - Harm['f']))
    # detected_points = np.where(Harm['f'] == Harm['f'][detected_idx])[0]
    # print(Harm['k'][detected_points])
    # ax.plot(Harm['k'][detected_points],Harm['A'][detected_points],'s',mec = 'k')

ax.set_yscale('log')
ax.set_ylim([1e-4,1e-1])
ax.legend()

figname = f'{fig_folder}Observation_several_freq_gaussian_filter_1'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Method for peaks detection
# step 1 - find peaks using gaussian convolution
# step 2 - filter and keep only relevant peaks
# step 3 - for each relevant peak, compute points on which gaussian will be fitted
# step 4 - save alpha coefficient for each detected peak 



#%% Select a profile for a given frequency and perform peak detection with gaussian fit

Afk = Efk
delta_k = np.diff(Efk['k'])[0]

frequencies = [0.2,0.9]
idx_min = np.argmin(np.abs(frequencies[0] - Afk['f']))
idx_max = np.argmin(np.abs(frequencies[1] - Afk['f']))
indices = np.arange(idx_min,idx_max,step = 1)

# Build Gaussian to be convoluted
gaussian_width = 2
N = 6
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)

# parameters to find peaks on correlation curve
min_prominence = 1e-2
S = {'idx':[],'peaks':[],'A':[],'f':[],'k':[]}
# create dictionnary in which results are saved 

fig, ax = plt.subplots()
for i in range(len(indices)):
    ax.clear()
    
    idx = indices[i]
    
    detected_f = Afk['f'][idx]
    print(detected_f)
    section = Afk['E'][idx,:]
    label_freq = f'{detected_f:.2f}'
    norm = section/section.max()
    
    ax.plot(Afk['k'],norm,'-o',label = label_freq)
    # x_gauss = np.arange(0,gaussian_width*N,step = 1)*delta_k
    # ax.plot(x_gauss,gaussian,'-k')
    
    # Apply gaussian filter 
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    ax.plot(Afk['k'],norm_correlated)
    # find peaks of correlation
    peaks,properties = scipy.signal.find_peaks(norm_correlated,prominence = min_prominence)
    ax.plot(Afk['k'][peaks],norm_correlated[peaks],'d')
    
    ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
    ax.set_ylabel(r'$A_{norm}$')
    ax.legend()
    
    # save parameters in dictionnary
    for i,p in enumerate(peaks) :
        S['idx'].append(idx)
        S['peaks'].append(p)
        S['A'].append(norm[p]*section.max())
        S['f'].append(detected_f)
        S['k'].append(Afk['k'][p])
    
    # f_txt = f'{detected_f:.2f}'.replace('.','p')
    # figname = f'{fig_folder}Section_f{f_txt}_corr_gaussian_sigma{gaussian_width}'
    # plt.savefig(figname + '.pdf', bbox_inches='tight')
    # plt.savefig(figname + '.png', bbox_inches='tight')

for key in S.keys():
    S[key] = np.array(S[key])
    
    
#%% Show detected peaks superposed with space-time spectrum

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(S['k'], S['f'],S)


#%%
fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_respectto_k_fixed_f.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(filtered_properties,pf)

#%% Select a peak and cascade down to compute width over which we perform fitting

rel_height = 0.4

m = {'A':[],'f':[],'k':[],'sigma':[],'err_k':[],'err_sigma':[]}

fig, ax = plt.subplots()
for idx in indices:

    ax.clear()
    
    detected_f = Afk['f'][idx]
    print(f'f  = {detected_f:.3f}')
    label_freq = f'{detected_f:.2f}'
    
    # select peaks detected at this frequency
    mask = np.where(filtered_properties['idx'] == idx)[0]
    peaks = filtered_properties['peaks'][mask]
    
    section = Afk['E'][idx,:]
    norm = section/section.max()
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    
    ax.plot(Afk['k'],norm,'-o',label = label_freq)
    ax.plot(Afk['k'],norm_correlated,'-o')
    
    for p in peaks:
        ax.plot(Afk['k'][p],norm_correlated[p],'s')
        
        points2fit = indices2fit(norm_correlated,Afk['k'],p,rel_height)
    
        # plot points
        ax.plot(Afk['k'][points2fit],norm[points2fit],'d')  
        if section[p] > section[points2fit].min():
            x = Afk['k'][points2fit]
            y_exp = (section[points2fit] - section[points2fit].min())/(section[p] - section[points2fit].min())
            
            # gaussian fit 
            popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : unnormed_gaussian(x, x0, sigma),x,y_exp,
                                                 bounds = ([x[0],1e-5],[x[-1],1]))
            xfit = np.linspace(x[0],x[-1],100)
            yth = unnormed_gaussian(xfit, popt[0], popt[1])
            err_coeff = np.sqrt(np.diag(pcov))
            print(f'k0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and sigma = {popt[1]:.2e} ± {err_coeff[1]:.2e}')
            # fig, ax1 = plt.subplots()
            # ax1.plot(x,y_exp,'-o')
            # ax1.plot(xfit,yth)  
            
            m['A'].append(Afk['E'][idx,p])
            m['f'].append(detected_f)
            m['k'].append(popt[0])
            m['sigma'].append(popt[1])
            m['err_k'].append(err_coeff[0])
            m['err_sigma'].append(err_coeff[1])
        
    plt.pause(0.5)
        
for key in m.keys():
    m[key] = np.array(m[key])
        
        
#%% DISCRIMINATE DIFFERENT HARMONICS
N = np.arange(1,4)
hw = 4.5
nb_points = len(m['k'])
k_list = np.linspace(m['k'].min(),m['k'].max(),nb_points)
harmonics = np.zeros((k_list.shape[0],N.shape[0]))
for i in range(len(N)):
    harmonics[:,i] = bound_harmonicN(k_list, N[i], hw)

fig, ax = plt.subplots()
ax.plot(m['k'],m['f'],'o')
for i in range(len(N)):
    ax.plot(k_list,harmonics[:,i]/2/np.pi)

# build a distance to each harmonic
D = np.zeros((nb_points,len(N)))
for i in range(len(N)):
    D[:,i] = np.abs(m['f'] - bound_harmonicN(m['k'], N[i], hw)/2/np.pi)
    
closest_harmonic = np.argmin(D,axis = 1)
# correct closestharmonic
closest_harmonic[np.where(m['k']<0.4)] = 0

mask = np.where(np.logical_and(closest_harmonic == 2,m['k'] < 0.67))
closest_harmonic[mask] = 1

m['closest_harmonic'] = closest_harmonic

#%% 
fig, ax = plt.subplots()
for i in range(len(N)):
    mask = m['closest_harmonic'] == i
    ax.scatter(m['k'][mask],m['f'][mask])

fig, ax = plt.subplots()
for i in range(3):
    mask = m['closest_harmonic'] == i
    alpha = np.sqrt(2)*m['sigma'][mask]
    ax.scatter(m['f'][mask],alpha) 
  
ax.set_xlabel(r'$f$')
ax.set_ylabel(r'$\alpha$')
ax.set_xscale('log')
ax.set_yscale('log')




#%% TRY PEAK DETECTION WITH RESPECT TO F? AT FIXED K

  


# f_txt = f'{detected_f:.2f}'.replace('.','p')
# figname = f'{fig_folder}Section_f{f_txt}_corr_gaussian_cascade'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')

###########################################################################################################################
###########################################################################################################################


#%%
# find width over which we perform the fit
rel_height = 0.5
widths = scipy.signal.peak_widths(norm_correlated,peaks,rel_height = rel_height)

max_width = 25
mask = np.where(widths[0] > max_width)[0]
new_widths = scipy.signal.peak_widths(norm_correlated,peaks[mask],rel_height = 0.32)
for i in range(len(widths)):
    widths[i][mask] = new_widths[i]

ax.hlines(y=widths[1], xmin=widths[2]*delta_k,
        xmax=widths[3]*delta_k,color = 'g')

    
#%% order peaks and right_ips / left_ips to compute best rel_height such that ips are correctly ordered for each peak






#%%


print('Select space-time spectrum from Python')



gaussian_width = 3
N = 6
min_prominence = 2e-3
width_max = 30 # in samples
rel_height = 0.3
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)

plt.ion()
fig, ax = plt.subplots()

S = {'idx':[],'peaks':[],'A':[],'f':[],'k':[],'w':[],'right_ips':[],'left_ips':[]}
for idx in range(idx_min,idx_max + 1):

    ax.clear()
    
    detected_f = Afk['f'][idx]
    print(detected_f)
    section = Afk['E'][idx,:]
    label_freq = f'{detected_f:.2f}'
    norm = section/section.max()
    ax.plot(norm,'-o',label = label_freq)

    # Perform correlation
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    peaks,properties = scipy.signal.find_peaks(norm_correlated,prominence = min_prominence)
    widths,width_heights,left_ips,right_ips = scipy.signal.peak_widths(norm_correlated,peaks,rel_height = rel_height)
    w_properties = {'widths':widths,'width_heights':width_heights,'left_ips':left_ips,'right_ips':right_ips}
    
    # keep only peaks with widths below a threshold
    mask = np.where(w_properties['widths'] < width_max)[0]
    peaks = peaks[mask]
    for key in properties.keys():
        properties[key] = properties[key][mask]
    
    for key in w_properties.keys():
        w_properties[key] = w_properties[key][mask] 
        
    # save data in main dictionnary
    for i,p in enumerate(peaks) :
        S['idx'].append(idx)
        S['peaks'].append(p)
        S['A'].append(norm[p]*section.max())
        S['f'].append(detected_f)
        S['k'].append(Afk['k'][p])
        S['w'].append(w_properties['widths'][i])
        S['right_ips'].append(w_properties['right_ips'][i])
        S['left_ips'].append(w_properties['left_ips'][i])
    
    
    ax.plot(norm_correlated,'-')
    ax.vlines(x=peaks, ymin=norm_correlated[peaks] - properties["prominences"],
           ymax = norm_correlated[peaks],color = 'g')
    ax.hlines(y=w_properties['width_heights'], xmin=w_properties['left_ips'],
           xmax=w_properties['right_ips'],color = 'g')
    ax.plot(peaks,norm[peaks],'s')
    ax.legend()

    plt.draw()
    plt.pause(0.3)
    
for key in S.keys():
    S[key] = np.array(S[key])


#%% Superpose detected peaks with space-time spectrum
plt.ion()
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['k'].min(),Afk['k'].max(),Afk['f'].min(),Afk['f'].max()))
ax.scatter(S['k'],S['f'],color = 'w')
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(0.09,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

#%% Filter detected peaks

outliers = np.logical_and(S['f'] > 0.32, S['k'] < 0.35)
outliers = np.logical_or(outliers,S['k']> 2.33)

# I need to create a function to filter points by hand (clicking in 4 different points)

fig, ax = plt.subplots()
ax.plot(S['k'][~outliers],S['f'][~outliers],'o')
            
for key in S.keys():
    S[key] = S[key][~outliers]

#%% Detect peaks and width over which gaussian can be fitted 

frequencies = [0.4,0.9]
idx_min = np.argmin(np.abs(frequencies[0] - Afk['f']))
idx_max = np.argmin(np.abs(frequencies[1] - Afk['f']))

delta_k = Efk['k'][1] - Efk['k'][0] 
width_max = 30 # in samples
rel_height = 0.5

plt.ion()
fig, ax = plt.subplots()

for idx in range(idx_min,idx_max + 1):

    ax.clear()
    
    detected_f = Afk['f'][idx]
    print(detected_f)
    section = Afk['E'][idx,:]
    label_freq = f'{detected_f:.2f}'
    norm = section/section.max()
    ax.plot(Afk['k'],norm,'-o',label = label_freq)
    
    # plot peaks
    mask = np.where(S['idx'] == idx)[0]
    # create a smaller structure for this frequency
    s = {}
    for key in S.keys():
        s[key] = S[key][mask]

    ax.plot(s['k'],s['A']/section.max(),'s')
    
    prominences,left_bases,right_bases = scipy.signal.peak_prominences(norm,s['peaks'])
    widths,width_heights,left_ips,right_ips = scipy.signal.peak_widths(norm,s['peaks'],rel_height = rel_height,
                                                                       prominence_data = (prominences,left_bases,right_bases))
    
    for peak_idx in range(len(s['peaks'])):

        # select points to fit using width of the correlation curve
        low_bound = left_ips[peak_idx]*delta_k
        up_bound = right_ips[peak_idx]*delta_k
        
        points2fit = np.where(np.logical_and(Afk['k'] > low_bound,Afk['k'] < up_bound))
        ax.plot(Afk['k'][points2fit],norm[points2fit],'gd')
         
    ax.legend()

    plt.draw()
    plt.pause(1.0)
    

#%% Perform gaussian fit over a given peak 

select_freq = 0.6
delta_k = Efk['k'][1] - Efk['k'][0] 


# compute associated idx
idx = np.argmin(np.abs(select_freq - Afk['f']))
# select peaks 
mask = np.where(S['idx'] == idx)[0]
# create a smaller structure for this frequency
s = {}
for key in S.keys():
    s[key] = S[key][mask]

section = Afk['E'][idx,:]
norm = section/section.max()
detect_f = s['f'][0]
label_freq = f'{detect_f:.2f}'

for peak_idx in range(len(s['A'])):

    # select points to fit using width of the correlation curve
    low_bound = s['left_ips'][peak_idx]*delta_k
    up_bound = s['right_ips'][peak_idx]*delta_k
    
    points2fit = np.where(np.logical_and(Afk['k'] > low_bound,Afk['k'] < up_bound))
    fig, ax = plt.subplots()
    ax.plot(Afk['k'],norm,'-o',label = label_freq)
    ax.plot(s['k'],s['A']/section.max(),'s')
    ax.plot(Afk['k'][points2fit],norm[points2fit],'gd')
    
    x = Afk['k'][points2fit]
    y_exp = (section[points2fit] - section[points2fit].min())/(s['A'][peak_idx] - section[points2fit].min())
    
    fig, ax = plt.subplots()
    ax.plot(x,y_exp,'-o')
    
    # gaussian fit 
    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,alpha : unnormed_gaussian(x, x0, alpha),x,y_exp,
                                         bounds = ([low_bound,1e-3],[up_bound,1]))
    xfit = np.linspace(x[0],x[-1],100)
    yth = unnormed_gaussian(xfit, popt[0], popt[1])
    err_coeff = np.sqrt(np.diag(pcov))
    print(f'k0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and sigma = {popt[1]:.2e} ± {err_coeff[1]:.2e}')
    ax.plot(xfit,yth)  
    

    
    
    
    
    
    
#%% perform a peak detection directly 
frequencies = [0.2,0.9]
idx_min = np.argmin(np.abs(frequencies[0] - Afk['f']))
idx_max = np.argmin(np.abs(frequencies[1] - Afk['f']))

prominence_limits = (5e-2,1.0)
width_limits = (3,100) # in samples number 
width_max = 15
rel_height = 0.4

plt.ion()
fig, ax = plt.subplots()

S = {'idx':[],'A':[],'f':[],'k':[],'w':[],'right_ips':[],'left_ips':[]}
for idx in range(idx_min,idx_max + 1):

    ax.clear()
    
    detected_f = Afk['f'][idx]
    print(detected_f)
    section = Afk['E'][idx,:]
    label_freq = f'{detected_f:.2f}'
    norm = section/section.max()
    ax.plot(norm,'-o',label = label_freq)

    peaks,properties = scipy.signal.find_peaks(norm,prominence = prominence_limits,width = width_limits)
    widths,width_heights,left_ips,right_ips = scipy.signal.peak_widths(norm,peaks,rel_height = rel_height)
    w_properties = {'widths':widths,'width_heights':width_heights,'left_ips':left_ips,'right_ips':right_ips}
    
    # keep only peaks with widths below a threshold
    mask = np.where(w_properties['widths'] < width_max)[0]
    peaks = peaks[mask]
    for key in properties.keys():
        properties[key] = properties[key][mask]
    
    for key in w_properties.keys():
        w_properties[key] = w_properties[key][mask] 
        
    # save data in main dictionnary
    for i,p in enumerate(peaks) :
        S['idx'].append(idx)
        S['A'].append(norm[p]*section.max())
        S['f'].append(detected_f)
        S['k'].append(Afk['k'][p])
        S['w'].append(w_properties['widths'][i])
        S['right_ips'].append(w_properties['right_ips'][i])
        S['left_ips'].append(w_properties['left_ips'][i])
      
    ax.vlines(x=peaks, ymin=norm[peaks] - properties["prominences"],
           ymax = norm[peaks],color = 'g')
    ax.hlines(y=w_properties['width_heights'], xmin=w_properties['left_ips'],
           xmax=w_properties['right_ips'],color = 'g')
    ax.plot(peaks,norm[peaks],'s')
    ax.legend()

    plt.draw()
    plt.pause(0.3)
    

#%% Show detected peaks

fig, ax = plt.subplots()

c = ax.imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['k'].min(),Afk['k'].max(),Afk['f'].min(),Afk['f'].max()))
ax.scatter(S['k'],S['f'],color = 'w')
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(0,3.3)
ax.set_ylim(0.09,1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

#%% Filter detected peaks

mask = np



#%% NEXT STEPS

# Keep track of correlation peak width 
# select a given width of the detected peak
# Fit gaussian on this detected peak (without prefactor) 







# Plot vertical line in the space-time spectrum 

# set_graphs.set_matplotlib_param('single')
# fig, ax = plt.subplots()
# c = ax.imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
#               origin = 'lower', interpolation = 'gaussian',
#               extent = (Afk['k'].min(),Afk['k'].max(),Afk['f'].min(),Afk['f'].max()))

# ax.axhline(detected_f,color = 'w')
# # ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlim(0.1,3.3)
# ax.set_ylim(0.09,1.3)

# ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
# ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

# cbar = plt.colorbar(c,ax = ax)
# cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

 



#%% Create a gaussian function 

































