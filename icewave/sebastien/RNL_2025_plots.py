# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:30:59 2025

@author: sebas
"""
import numpy as np
import cmath
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize

import cv2 as cv
import glob
import os
import pickle 
import h5py
import re

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

# import seb plotting package 
import icewave.sebastien.set_graphs as set_graphs

#%% Set plotting parameters 
 
font_size_medium = 40
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
marker_size_plot = 12

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
# PARULA COLORMAP 
parula_map = matcmaps.parula()

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Conferences/GDR_Geofluides_2025/Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 

def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN

def power_law(x,A,alpha):
    return A*x**(alpha)

def exponential(x,A,alpha):
    return A*np.exp(alpha*x)

#%% Figures associated to Aerial observations 

# Load data 
# base = 'E:/sauvegarde_partielle_hublot_20250206/Share_hublot/Data/'
base = 'F:/Rimouski_2024/Data/'

date = '0226'
drone_ID = 'mesange'
exp_ID = '10-waves_005'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = path2data + 'PIV_processed_i00_Dt5_b1_W32_xROI600_width3240_yROI1_height2159_scaled.mat'

with h5py.File(file2load, 'r') as fmat:
    data_wave = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])

#%% Load image 

path2img = path2data + 'DJI_20240226191555_0648_D_exemple.tiff'

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
img = np.flip(img,axis = 0)

fig, ax = plt.subplots(figsize = fig_size)
ax.imshow(img,aspect = 'equal',origin = 'lower', extent = [S['x'][0],S['x'][-1],S['y'][0],S['y'][-1]])
ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)

figname = f'{fig_folder}example_MIZ_scaled_{date}_{drone_ID}_{exp_ID}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


#%% Show a frame of DIC and fit it by a quadratic function 

V = np.transpose(S['Vx'],(1,2,0))
Vs = FT.supress_quadratic_noise(V, S['x'], S['y'])

fig, ax = plt.subplots()
c = ax.imshow(Vs[:,:,1], origin = 'lower', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]))
plt.colorbar(c)

#%% Compute Fourier spectrum of the velocity field 

fps = S['SCALE']['facq_t']
TF_spectrum,freq_time,FFT_t = FT.temporal_FFT(Vs,fps,padding_bool = 1,add_pow2 = 0,output_FFT = True)
#%%
fig, ax = plt.subplots(figsize = fig_size)
ax.loglog(freq_time,abs(TF_spectrum),'-',zorder = 2)
ax.set_ylim(1e-3,1)
ax.set_xlim(1e-3,100)
ax.grid(zorder = 1)

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V_x}| \rangle _{x,y} \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}TF_spectrum_{date}_{drone_ID}_{exp_ID}'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%% Show demodulated field for a given frequency 

f_demod = 0.31 # frequency chosen for demodulation 
idx = np.argmin(abs(freq_time - f_demod))
print(f'{freq_time[idx]:.3f}')
field = np.real(FFT_t[:,:,idx])

fig, ax = plt.subplots(figsize = fig_size)
c = ax.imshow(field, cmap = parula_map, aspect = 'equal', origin = 'lower', 
              interpolation = 'gaussian', extent = (S['x'][0], S['x'][-1],S['y'][0], S['y'][-1]),vmin = -1.2,vmax = 1.2)
cbar = plt.colorbar(c, ax = ax, shrink = 0.72)

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(m)}$',labelpad = 5)
cbar.set_label(r'$\hat{V}_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)

f_txt = f'fedmod_{f_demod:.3f}'
f_txt = f_txt.replace('.','p')
figname = f'{fig_folder}demodulated_field_{date}_{drone_ID}_{exp_ID}_{f_txt}'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%% Fit exponential to the given demodulated field 

profile = np.mean(FFT_t[:,:,idx],axis = 0)
profile_abs = np.mean(abs(FFT_t[:,:,idx]),axis = 0)
fig, ax = plt.subplots(figsize = fig_size)
# add a phase to demodulated field
phase = np.pi*0
profile_phase = profile*np.exp(1j*phase)

ax.plot(S['x'],np.real(profile_phase))
# ax.plot(S['x'],abs(profile))

# ax.plot(S['x'],np.mean(abs(FFT_t[:,:,idx]),0))
logA = np.log(abs(profile))
popt,pcov = np.polyfit(S['x'],logA,deg = 1,cov = True)
yth = np.exp(np.polyval(popt,S['x']))
label_expfit = r'$y(x) = ' f'{np.exp(popt[1]):.2f}' + r'e^{' + f'{popt[0]:.4f}' + 'x}$'
ax.plot(S['x'],yth,label = label_expfit)

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V_x}| \rangle _{y} (x,f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax.legend()

freq_txt = f'{freq_time[idx]:.4f}'
freq_txt = freq_txt.replace('.','p')
figname = f'{fig_folder}Amplitude_profile_fdemod_{freq_txt}_python'

plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')


#%% Compute E(f,k)


N = Vs.shape[2]
Efk = FT.space_time_spectrum(Vs,1/S['SCALE']['fx'],S['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%%

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
Amin = 1e-5 # change to adjust colormap
Amax = 1e-2 # change to adjust colormap
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),2*np.pi*Efk['f'].min(),2*np.pi*Efk['f'].max()))

ax.set_xscale('log')
ax.set_yscale('log')
kbounds = [0.1,3.3] # bounds for k axis on Efk plot
fbounds = [0.55,8] # bounds for f axis on Efk plot
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

















#%% Load data of detected harmonics 

path2harmonics = path2data + 'Plots/' 
file2load = path2harmonics + f'Data_plot_selected_harmonics_{date}_{drone_ID}_{exp_ID}.mat'

with h5py.File(file2load) as fmat : 
    
    print('Top-level keys : ', list(fmat.keys()))
    Harm = mat2py.mat_to_dict(fmat['S_disp'],fmat['S_disp'])
    
    
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


#%% Load E(f,k) from Matlab 
path2Efk = path2data + 'Plots/'
file2load = path2Efk + f'Data_A_fk_{date}_{drone_ID}_{exp_ID}.mat'

with h5py.File(file2load, 'r') as fmat:
    Afk = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    Afk['E'] = mat2py.mat_to_dict(fmat['E'],fmat['E'])
    Afk['f'] = mat2py.mat_to_dict(fmat['f'],fmat['f'])
    Afk['k'] = mat2py.mat_to_dict(fmat['k'],fmat['k'])
    Afk['omega'] = mat2py.mat_to_dict(fmat['omega'],fmat['omega'])
    
Afk['E'] = np.transpose(Afk['E'],axes = (1,0))


#%% Plot E(f,k) using Matlab 

k_list = np.linspace(0.08,5,100)
h_w = 4.07 # water depth on field, in meter
harmonic1 = bound_harmonicN(k_list, 1, h_w)

fig, ax = plt.subplots(figsize = fig_size)
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


figname = f'{fig_folder}Efk_{date}_{drone_ID}_{exp_ID}'

plt.savefig(f'{figname}.png',bbox_inches = 'tight',dpi = 300)
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight',dpi = 300)

#%% Plot section of spatio-temporal for a given frequency

f_demod = 0.3
idx = np.argmin(abs(f_demod - Afk['f']))
print(f'f = {f_demod:.3f}')


section = Afk['E'][idx,:]
fig,ax = plt.subplots(figsize = fig_size)
ax.loglog(Afk['k'],section,'o')





#%% Load Amplitude vs x profile for each frequencies 
base = 'E:/sauvegarde_partielle_hublot_20250206/Share_hublot/Data/'

date = '0226'
drone_ID = 'mesange'
exp_ID = '23-waves_012'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/Plots/'
file2load = f'{path2data}Data_Ax_0226_mesange_23-waves_012_0p15Hz_to_0p6Hz_harmonic_1.mat'

with h5py.File(file2load, 'r') as fmat:
    data_wave = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    S_Ax = mat2py.mat_to_dict(fmat['S_Ax'],fmat['S_Ax'])


#%% Select a frequency and plot the associated profile 

selected_freq = 0.35853
idx = np.argmin(abs(selected_freq - S_Ax['f']))
freq = S_Ax['f'][idx]
print(f'Selected frequency : {freq:.4f}')

A = S_Ax['S_A']['A'][idx][0]
x = S_Ax['S_A']['x'][idx][0]

# fit by exponential
Alog = np.log(A)
popt,pcov = np.polyfit(x,Alog,1,cov = True)

yth = np.exp(np.polyval(popt,x))
# popt,pcov = scipy.optimize.curve_fit(lambda x,A,alpha : exponential(x,A,alpha),x,A)
# yth = exponential(x,popt[0],popt[1])

fig, ax = plt.subplots(figsize = fig_size)
ax.plot(x,A,'o',markersize = marker_size_plot,markeredgecolor = 'k')
label_theory = r'$y(x) = ' + f'{np.exp(popt[1]):.2f}' + 'e^{' + f'{popt[0]:.3f}x' + '}$'
ax.plot(x,yth,'r',label = label_theory)

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V_x}| \rangle _{y} (x,f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax.legend()

freq_txt = f'{freq:.4f}'
freq_txt = freq_txt.replace('.','p')
figname = f'{fig_folder}Amplitude_profile_fdemod_{freq_txt}'

plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%% Load data from one experiment Aurore set up 

fx = 0.337 #scaling in mm/pix
facq_x = 1/fx
facq_t = 67 
# facq_t = 37.215 # frame rate 
ft = 1/facq_t # temporal scaling in s/frame

root = 'Y:/Banquise/Aurore/amplitude_vagues_laser/'

e = 5.0 # frasil thickness in mm
key_thickness = str(e)
thck_string = f'{e}mm/'
path2data = f'{root}{thck_string}'

#Load dictionnary of data 
units_dict = {'alpha':'mm-1','f_ex':'Hz','err_alpha':'mm-1','A0':'mm-1','err_A0':'mm-1'}

filename_data = f'{root}data_attenuation_coeff_e_f_ex.pkl'
if not os.path.isfile(filename_data):
    print('Data dictionnary does not already exist')
    data = {}
    data[key_thickness] = {'f_ex':[],'f0': [],'alpha':[],'err_alpha':[],'A0':[],'err_A0':[],
                           'units':units_dict}
else :
    print('Data dictionnary already exists')
    with open(filename_data,'rb') as pf :
        data = pickle.load(pf)
        
    if key_thickness not in data.keys():
        print(f'Thichkness e = {e} mm does not exist in data dictionnary')
        data[key_thickness] = {'f_ex':[],'f0':[],'alpha':[],'err_alpha':[],'A0':[],'err_A0':[],
                               'units':units_dict}
    else :
        print(f'Thichkness e = {e} mm already exists in data dictionnary')

#%% Loop over all frequencies

xstart = 30 # in mm, start of relevant part of the signal

explist = glob.glob(f'{path2data}*')
# get data for a given frequency
for exp_idx in range(0,len(explist)):
    exp_name = explist[exp_idx]
    
    filelist = glob.glob(f'{exp_name}/*.pkl')
    
    f_ex = float(re.findall(r'(\d+\.\d+)Hz',exp_name)[0]) # excitation_frequency 
    print(f'Excitation frequency = {f_ex} Hz')
    
    if f_ex >= 4.2 :
        facq_t = 37.215
        
    pickle_file = filelist[0]
    
    with open(pickle_file,'rb') as pf :
        Q = pickle.load(pf)
        
    x = np.arange(0,np.shape(Q)[0],step = 1)*fx # distance x in mm 
    t = np.arange(0,np.shape(Q)[1],step = 1)*ft
    
    # Keep only relevant part of the profile 
    
    istart = np.argmin(abs(x - xstart))
    Q_new = Q[istart:,:]
    x_new = x[istart:]
    # compute mean value for each position 
    mean_profile = np.mean(Q_new,axis = -1)
    pmean = np.polyfit(x_new,mean_profile,1)
    interp_mean = np.polyval(pmean,x_new) # interpolate mean value along x
    
    # Perform demodulation 
    
    demod_profile,f_demod = FT.get_demodulated_main_freq(Q_new,t,facq_t)
    
    fig, ax = plt.subplots(figsize = fig_size)
    ax.plot(x_new,np.real(demod_profile))
    # ax.plot(x_new,abs(demod_profile))
    
    # fit demodulated profile by an exponential 
    logA = np.log(abs(demod_profile))
    p,V = np.polyfit(x_new,logA,deg = 1,cov = True)
    yfit = np.exp(np.polyval(p, x_new))
    
    A_0 = np.exp(p[-1])
    alpha = abs(p[0])
    err_alpha = np.sqrt(V[0,0])
    err_A0 = np.sqrt(np.exp(V[1,1]))
    
    label_expfit = r'$y(x) = {:.2f}'.format(A_0) + r'e^{{{:.4f}x}}$'.format(alpha)
    ax.plot(x_new,yfit,label = label_expfit)
    
    # set correctly the plot and save it
    ax.set_xlabel(r'$x \; \mathrm{(mm)}$')
    ax.set_ylabel(r'$\hat{\xi} \; \mathrm{(mm)}$')
    ax.legend()
    
    file2save = f'demod_profile_e{e}_f_demod{f_demod}'
    file2save = file2save.replace('.','p')
    file2save = f'{fig_folder}{file2save}'
    plt.savefig(file2save + '.pdf', bbox_inches='tight')
    plt.savefig(file2save + '.png', bbox_inches='tight')
    
    print(f'Main frequency is {f_demod} Hz')
    plt.close('all')





