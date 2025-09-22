# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 14:28:51 2024

@author: sebas
"""

import numpy as np
import cmath
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize

import cv2 as cv
import glob
import os
import pickle 
import h5py

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

# import seb plotting package 
import icewave.sebastien.set_graphs as set_graphs

#%% Parameters for plots
# picture on whole page width : fig_size = (8,6), font_size_medium = 22, font_size_small = 0.75*font_size_medium

# picture on half page width : fig_size = (8,6), font_size_medium = 30, font_size_small  =  0.75*font_size_medium

# quadrant of 4 figures : fig_size = (12,9), font_size_medium = 20, font_size_small = 0.75*font_size_medium 

 
font_size_medium = 20
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 100 # dpi to save images 

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
# PARULA COLORMAP 
parula_map = matcmaps.parula()

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Figures_article_BicWin2024/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% Import data 
base = 'F:/Rimouski_2024/Data/'

date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = path2data + 'PIV_processed_i00_N0_Dt6_b1_W32_xROI1_width3839_yROI1_height2159_scaled.mat'

with h5py.File(file2load, 'r') as fmat:
    data_wave = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])

#%% Save data as pkl file 

file2save = file2load.replace('.mat','.pkl')

with open(file2save,'wb') as pfile :
    pickle.dump(S, pfile)
    
print('Data saved as pkl file')

#%% FUNCTION SECTION 

def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN


    
#%% Show a frame of DIC and fit it by a quadratic function 

V = np.transpose(S['Vx'],(1,2,0))
Vs = FT.supress_quadratic_noise(V, S['x'], S['y'])

fig, ax = plt.subplots()
c = ax.imshow(Vs[:,:,1], origin = 'lower', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]))
plt.colorbar(c)

#%% Compute Fourier spectrum of the velocity field 

fps = S['SCALE']['facq_t']
TF_spectrum,freq_time,FFT_t = FT.temporal_FFT(V,fps,padding_bool = 1,add_pow2 = 0,output_FFT = True)

fig, ax = plt.subplots()
ax.loglog(freq_time,abs(TF_spectrum),'-',zorder = 2)
ax.set_ylim(1e-4,10)
ax.grid(zorder = 1)

#%% Show demodulated field for a given frequency 

f_demod = 0.212 # frequency chosen for demodulation 
idx = np.argmin(abs(freq_time - f_demod))
field = np.real(FFT_t[:,:,idx])

fig, ax = plt.subplots()
c = ax.imshow(field, cmap = parula_map, aspect = 'equal', origin = 'lower', 
              interpolation = 'gaussian', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]))
cbar = plt.colorbar(c, ax = ax, shrink = 0.6)

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)
cbar.set_label(r'$\hat{V}_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)

fig, ax = plt.subplots()
ax.plot(S['x'],field[20,:])

#%% Compute space-time spectrum 

facq_x = 1/S['SCALE']['fx']
fps = S['SCALE']['facq_t']

result = FT.space_time_spectrum(Vs, facq_x, fps) 

#%% Show space FFT for a given frequency

f_demod = 0.25
idx = np.argmin(abs(result['freq'] - f_demod))

fig, ax = plt.subplots()
c = ax.imshow(np.abs(result['shift'][:,:,idx]), cmap = parula_map, aspect = 'equal', origin  = 'lower',
          extent = (result['kx'][0],result['kx'][-1],result['ky'][0],result['ky'][-1]))
# ax.plot(kx[round(x0)],ky[round(y0)],'o',mfc = 'red')
cbar = plt.colorbar(c, ax = ax)

#%% Load E(f,k) from Matlab 
path2Efk = path2data + 'Plots/'
file2load = path2Efk + 'Data_A_fk_0226_mesange_23-waves_012.mat'

with h5py.File(file2load, 'r') as fmat:
    Afk = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    Afk['E'] = mat2py.mat_to_dict(fmat['E'],fmat['E'])
    Afk['f'] = mat2py.mat_to_dict(fmat['f'],fmat['f'])
    Afk['k'] = mat2py.mat_to_dict(fmat['k'],fmat['k'])
    Afk['omega'] = mat2py.mat_to_dict(fmat['omega'],fmat['omega'])
    
Afk['E'] = np.transpose(Afk['E'],axes = (1,0))

#%%

k_list = np.linspace(0.08,5,100)
h_w = 5.0 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, h_w)

fig, ax = plt.subplots()
c = ax.imshow(result['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (result['k'].min(),result['k'].max(),2*np.pi*result['freq'].min(),
                        2*np.pi*result['freq'].max()))
ax.plot(k_list,harmonic1,'--',color = 'white',lw = 2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(2*np.pi*0.09,2*np.pi*1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)


#%% Plot Data from Matlab 

k_list = np.linspace(0.08,5,100)
h_w = 5.0 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, h_w)

fig, ax = plt.subplots()
c = ax.imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['k'].min(),Afk['k'].max(),Afk['omega'].min(),Afk['omega'].max()))
ax.plot(k_list,harmonic1,'--',color = 'white',lw = 2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(2*np.pi*0.09,2*np.pi*1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

#%% Plot Data using Seb plotting package

k_list = np.linspace(0.08,5,100)
h_w = 5.0 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, h_w)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['k'].min(),Afk['k'].max(),Afk['omega'].min(),Afk['omega'].max()))
ax.plot(k_list,harmonic1,'--',color = 'white',lw = 2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,3.3)
ax.set_ylim(2*np.pi*0.09,2*np.pi*1.3)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)



#%% Load initial image 

# path2img = path2data + 'DJI_20240226213559_0938_D_exemple.tiff'
path2img = 'K:/Share_hublot/Data/0226/Drones/mesange/23-waves_012/DJI_20240226213559_0938_D_exemple.tiff'

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
img = np.flip(img,axis = 0)

fig, ax = plt.subplots(figsize = fig_size)
ax.imshow(img,aspect = 'equal',origin = 'lower', extent = [S['x'][0],S['x'][-1],S['y'][0],S['y'][-1]])
ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)

#%% Load data of detected harmonics 

path2harmonics = path2data + 'Plots/' 
file2load = path2harmonics + 'Data_plot_selected_harmonics_0226_mesange_23-waves_012.mat'

with h5py.File(file2load) as fmat : 
    
    print('Top-level keys : ', list(fmat.keys()))
    Harm = mat2py.mat_to_dict(fmat['S_disp'],fmat['S_disp'])
    
#%% Load water height associated to drone GPS position 


#%% Plot Harmonics with colorbar 

k_th = np.linspace(0,2,100)
omega1 = bound_harmonicN(k_th, N = 1, h_w = 5.0)


fig, ax = plt.subplots()
scatter = ax.scatter(Harm['k'][Harm['closest_harmonic'] == 1], Harm['omega'][Harm['closest_harmonic'] == 1],
                     c = Harm['A'][Harm['closest_harmonic'] == 1], s = 70, cmap = new_blues,
                     edgecolors = 'k',zorder = 2)

ax.plot(k_th,omega1,'--', color = 'blue',zorder = 1)

cbar = plt.colorbar(scatter,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5.0)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)

N = 1
hw_min = scipy.optimize.curve_fit(lambda x,h_w : bound_harmonicN(x,N,h_w),Harm['k'][Harm['closest_harmonic'] == 1], 
                                   Harm['omega'][Harm['closest_harmonic'] == 1])

print(f'Best fit : hw = {hw_min[0]} +- {2*np.sqrt(hw_min[1])}')

#%% Create quadrant figure for article 

fig, ax = plt.subplots(2,2,figsize = (12,9),height_ratios= [1,0.6], width_ratios= [1,1],
                       layout = 'constrained') #,gridspec_kw = {'width_ratios':[1,1],'height_ratios':[1.5,1]})
# Image
img_window = (0,0)
ax[img_window].imshow(img,aspect = 'equal',origin = 'lower', extent = [S['x'][0],S['x'][-1],S['y'][0],S['y'][-1]])
ax[img_window].set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax[img_window].set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)


# Demodulated field 
demod_window = (0,1)
f_demod = 0.212 # frequency chosen for demodulation 
idx = np.argmin(abs(freq_time - f_demod))
field = np.real(FFT_t[:,:,idx])

c = ax[demod_window].imshow(field, cmap = parula_map, aspect = 'equal', origin = 'lower', vmin = -2.2, vmax = 2.2,
              interpolation = 'gaussian', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]))

# Use make_axes_locatable to adjust the colorbar
# divider = make_axes_locatable(ax[demod_window])
# cax = divider.append_axes("right", size="4%", pad=0.2)
# cbar = fig.colorbar(c, cax=cax)

cbar = plt.colorbar(c,ax = ax[demod_window],shrink = 0.51)
ax[demod_window].set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax[demod_window].set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)
cbar.set_label(r'$\hat{V}_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)

# E(k,f)
Efk_window = (1,0)
k_list = np.linspace(0.08,5,100)
h_w = 5.0 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, hw_min[0])


c = ax[Efk_window].imshow(Afk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 3e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['k'].min(),Afk['k'].max(),Afk['omega'].min(),Afk['omega'].max()))
ax[Efk_window].plot(k_list,harmonic1,'--',color = 'white',lw = 2)
ax[Efk_window].set_xscale('log')
ax[Efk_window].set_yscale('log')
ax[Efk_window].set_xlim(0.1,3.1)
ax[Efk_window].set_ylim(2*np.pi*0.09,2*np.pi*1.3)

ax[Efk_window].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax[Efk_window].set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax[Efk_window])
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# Detected branch 
harmonic_window = (1,1)
k_th = np.linspace(0,2,100)
omega1 = bound_harmonicN(k_th, N = 1,h_w = hw_min[0])

scatter = ax[harmonic_window].scatter(Harm['k'][Harm['closest_harmonic'] == 1], Harm['omega'][Harm['closest_harmonic'] == 1],
                     c = Harm['A'][Harm['closest_harmonic'] == 1], s = 70, cmap = new_blues,
                     edgecolors = 'k',zorder = 2)

ax[harmonic_window].plot(k_th,omega1,'--', color = 'royalblue',zorder = 1)

cbar = plt.colorbar(scatter,ax = ax[harmonic_window])
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)
ax[harmonic_window].set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax[harmonic_window].set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)
ax[harmonic_window].set_xlim(0,1.7)
ax[harmonic_window].set_ylim(0,4.4)


#%% Create a figure using gridspec
use_omega = 0
if use_omega :
    factor_2pi = 1
else:
    factor_2pi = 2*np.pi 

fig = plt.figure(figsize = (12,9),layout = 'constrained')

# Define a GridSpec layout
gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios = [1,0.6], width_ratios = [1,1], wspace=0, hspace=0)

# add image 
ax1 = fig.add_subplot(gs[0,0])
c = ax1.imshow(img,aspect = 'equal',origin = 'lower', extent = [S['x'][0],S['x'][-1],S['y'][0],S['y'][-1]])
ax1.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax1.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)


# add demodulated field 
ax2 = fig.add_subplot(gs[0,1])

f_demod = 0.212 # frequency chosen for demodulation 
idx = np.argmin(abs(freq_time - f_demod))
field = np.real(FFT_t[:,:,idx])

c = ax2.imshow(field, cmap = parula_map, aspect = 'equal', origin = 'lower', vmin = -2.2, vmax = 2.2,
              interpolation = 'gaussian', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]))

# Use make_axes_locatable to adjust the colorbar

# divider = make_axes_locatable(ax2)
# cax = divider.append_axes("right", size="4%", pad=0.2)
# cbar = fig.colorbar(c, cax=cax)
cbar = plt.colorbar(c,ax = ax2,shrink = 0.435)
cbar.set_label(r'$\hat{V}_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)

ax2.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax2.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)

# add E(k,f)

k_list = np.linspace(0,5,100)
h_w = 5.0 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, hw_min[0])/factor_2pi

ax3 = fig.add_subplot(gs[1,0])
c = ax3.imshow(Afk['E'].T, cmap = parula_map , aspect = 'auto', norm = 'log', vmin = 6e-4, vmax = 1e-1,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Afk['omega'].min()/factor_2pi,Afk['omega'].max()/factor_2pi,Afk['k'].min(),Afk['k'].max()))
ax3.plot(harmonic1,k_list,'--',color = 'white',lw = 2)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_ylim(0.07,2.2)
ax3.set_xlim(0.5/factor_2pi,6/factor_2pi)


ax3.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax3.set_xlabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax3)
cbar.set_label(r'$|\hat{V}_x| (k,f) \; \mathrm{(u.a.)}$',labelpad = 5)

# get axis position 
# pos3 = ax3.get_position().bounds
# pos3[0] = pos3[0] + 0.05

# add detected harmonics

k_th = np.linspace(0,5,100)
omega1 = bound_harmonicN(k_th, N = 1, h_w = hw_min[0])/factor_2pi

ax4 = fig.add_subplot(gs[1,1])
scatter = ax4.scatter(Harm['omega'][Harm['closest_harmonic'] == 1]/factor_2pi, Harm['k'][Harm['closest_harmonic'] == 1],
                     c = Harm['A'][Harm['closest_harmonic'] == 1], s = 70, cmap = new_blues,
                     edgecolors = 'k',zorder = 1)

ax4.plot(omega1,k_th,'--', color = 'red',lw = 2,zorder = 2)

cbar = plt.colorbar(scatter,ax = ax4)
cbar.set_label(r'$|\hat{V}_x| (k,f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax4.set_ylabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax4.set_xlabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

ax4.set_ylim(0.07,2.3)
ax4.set_xlim(0.5/factor_2pi,6/factor_2pi)
ax4.set_xscale('log')
ax4.set_yscale('log')

# figname = fig_folder + 'Subplot_Drone_wave_field_kf'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')


#%% Fit data with shallow water dispersion relation 









