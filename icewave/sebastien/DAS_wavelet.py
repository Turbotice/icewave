# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 14:51:42 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
from time import strftime, localtime
from datetime import datetime
import pytz
import time 
import scipy
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pickle

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

global g
g = 9.81



#%% FUNCTION SECTION 

# function for plotting spatio temporal
def plot_spatio_temp(spatio,t,s,fiber_length,extents):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - t, numpy array or list, time array 
            - s, numpy array or list, curvilinear coordinate array
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsho, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#---------------------------------------------------------------------------------------------------

def wavenumbers_hydro(freq,rho_w,rho_ice,E,h,nu,H,c_w,equation):
    """ Compute wavenumbers associated to hydroelastic waves """ 
    
    g = 9.81
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    k = np.linspace(1e-4,10,200000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for i in range(len(freq)):
        omeg = 2*np.pi*freq[i]
        if omeg == 0:
            flag = 1
            idx_flag = i
        else:
            
            if equation == 'Squire_deep':
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + 1) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Squire_shallow':
                coth = 1/np.tanh(H*k)
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + coth) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Stein':
                cph = omeg/k
                func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            else : 
                print('inappropriate equation name, choose between : Squire_deep / Squire_shallow / Stein')

            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[i] = (np.where(np.diff(np.signbit(func)))[0]) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero] # wave vector associated to flexural mode 
    if flag:
        k_QS[idx_flag] = 0

    return k_QS  

#-----------------------------------------------------------------------------------------------------------------------------

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega





#%% load DAS parameters and data 

date = '0211'
# fig_folder = f'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2025/DAS/{date}/Figures/'


# Load parameters for DAS
path2DAS_param = 'U:/Data/parameters_Febus_2025.pkl'
with open(path2DAS_param,'rb') as pf:
    param = pickle.load(pf)
print('Parameters file loaded')

# Set parameters
# fs = np.shape(strain_rate)[1] # time sampling frequency 
# facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
fs = param[date]['fs']
fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
facq_x = param[date]['facq_x'] 

format_date = '%Y-%m-%d %H:%M:%S.f'
local_timearea = pytz.timezone('America/Montreal')
UTC_timearea = pytz.timezone('UTC')

# path2data = f'E:/Data/{date}/DAS_h5/'
# path2data = f'F:/Rimouski_2025/Data/{date}/DAS/'
path2data = f'U:/Data/{date}/DAS/'

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(path2data + '*.h5')
Nb_minutes = 1 # duration of each stack


idx_file = 1
file2load = filelist[idx_file]
with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
        
    a_group_key = list(f.keys())[0]
    data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])
            
# shape strain rate : 
# dim 0 : time in second / dim 1 : time sampled at fs / dim 2 : space 
strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
t = data['Source1']['time'] #time since epoch

# Stack strain rate into sections of duration Nb_minutes
stack_strain,stack_time,s = DS.new_stack_strain_rate(strain_rate, t, fiber_length,Nb_minutes)

stack_epoch = t[0] + stack_time
UTC_stack = np.empty(stack_time.shape,dtype = 'object')
for i in range(stack_epoch.shape[0]):
    current_UTC = DS.epoch2datetime(stack_epoch[i,:],timezone = UTC_timearea)
    UTC_stack[i,:] = current_UTC
        
#%% Plot spatio 

set_graphs.set_matplotlib_param('single')
chunk = 6
extents = [UTC_stack[chunk,0], UTC_stack[chunk,-1],s[0],s[-1]]
fig, ax , imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:],UTC_stack[chunk,:],s,fiber_length,extents = extents)

imsh.set(clim = [-1e4, 1e4])

# figname = f'spatio_file_{idx_file}_chunk{chunk}_{Nb_minutes}min'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute FFT in time 

N = 2**(FT.nextpow2(stack_strain.shape[1]))

mean_sig = np.mean(stack_strain,axis = 1)
mean_sig = np.tile(mean_sig,(stack_strain.shape[1],1,1)).transpose((1,0,2))

FFT_t = np.fft.fft(stack_strain - mean_sig,n = N,axis = 1)
FFT_t = FFT_t/stack_strain.shape[1]

# keep only half of frequencies
FFT_t = FFT_t[:,:N//2,:]
FFT_t[:,1:-1,:] = 2*FFT_t[:,1:-1,:]

freq = fs*np.arange(0,(N/2))/N

fig, ax = plt.subplots()
imsh = ax.imshow(abs(FFT_t[chunk,:,:]).T,origin = 'lower',cmap = parula_map, aspect = 'auto',
          extent = [freq[0],freq[-1],s[0],s[-1]])
imsh.set_clim([1e1,2e3])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$s \; \mathrm{(m)}$')
ax.set_xlim([0,10])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

# figname = f'FFT_t_file_{idx_file}_chunk{chunk}_{Nb_minutes}min'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Show a profile for a frequency 

selected_freq = 0.2
idx_freq = np.argmin(abs(freq - selected_freq))

profile = FFT_t[0,idx_freq,:]

fig, ax = plt.subplots()
ax.plot(s,np.real(profile))
ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$\hat{\dot{\epsilon}}(s)$')

#%% Create a spatial Short-Time Fourier Transform 

# create hann window
M = 100 # number of points of window (in samples)
hann_win = scipy.signal.windows.hann(M,sym = True)

# create SFT transfrom
hop = 25 # hops between two computations (in samples)
padding = 2**(FT.nextpow2(M) + 1)
SFT = scipy.signal.ShortTimeFFT(hann_win, hop,facq_x, fft_mode = 'onesided',mfft = padding,scale_to = 'magnitude')


#%% Plot spatial spectrogram for a given frequency 

Sx = SFT.spectrogram(np.real(profile),detr = 'constant',axis = -1)

# Plot spectrogram 
spectro_log = np.log10(abs(Sx))
fig, ax = plt.subplots()
imsh = ax.imshow(spectro_log,origin = 'lower',aspect = 'auto',extent = SFT.extent(stack_strain.shape[2]),
          cmap = parula_map)

imsh.set_clim([1,6])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(|\dot{\epsilon}|)$')

ax.set_xlabel(r'$s \; \mathrm{(m)}$')
ax.set_ylabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')

# figname = f'Spatial_spectro_file_{idx_file}_chunk{chunk}_{Nb_minutes}min_fdemod_{freq[idx_freq]}_M{M}_hop{hop}'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compute Spectrogram for all frequencies 

chunk = 0
Sfx = SFT.spectrogram(np.real(FFT_t[chunk,:,:]),detr = 'constant',axis = -1)

x = SFT.t(stack_strain.shape[2])

selected_x = 80 # distance from interogator
idx_x = np.argmin(abs(selected_x - SFT.t(stack_strain.shape[2])))
print(idx_x)

SFT_extents = SFT.extent(stack_strain.shape[2])

fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(Sfx[:,:,idx_x]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [SFT_extents[2],SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([1,6])

# ax.set_xlim([0,0.1])
# ax.set_ylim([0,1])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'FK for $s = ' +  f'{x[idx_x]:.2f}' + r'\; \mathrm{m}$')

# figname = f'FK_file_{idx_file}_s_{x[idx_x]}m_chunk{chunk}_{Nb_minutes}min_M{M}_hop{hop}_long_range'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Average FK over all chunks 

mean_FFT = np.mean(FFT_t,axis = 0)
mean_Sfx = SFT.spectrogram(np.real(mean_FFT),detr = 'constant',axis = -1)

#%% Plot FK averaged over time for a specific abscisse

selected_x = 80 # distance from interogator
idx_x = np.argmin(abs(selected_x - SFT.t(stack_strain.shape[2])))
print(idx_x)

SFT_extents = SFT.extent(stack_strain.shape[2])

fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx_x]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [SFT_extents[2],SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([1,5])

# ax.set_xlim([0,0.1])
# ax.set_ylim([0,1.0])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'Average FK for $s = ' +  f'{x[idx_x]:.2f}' + r'\; \mathrm{m}$')


# Select points on the graph
print('Select points on FK plots')
Nb_points = 10
points = plt.ginput(Nb_points)

sigma,f_exp = zip(*points)
k_exp = 2*np.pi*np.array(sigma)
ax.plot(sigma,f_exp,'ro')

# fit experimental points by shallow hydroelastic 
rho_w = 1027
H = 4
popt,pcov = scipy.optimize.curve_fit(lambda x,D : shallow_hydroelastic(x, D, rho_w, H)/2/np.pi,k_exp,f_exp,
                                     bounds = (1e5,1e8))

err_coeff = np.sqrt(np.diag(pcov))
print(f'D = {popt[0]:.2e} ± {err_coeff[0]:.2e}')

k_th = np.linspace(0,1,100)
omega_th = shallow_hydroelastic(k_th, popt, rho_w, H)
label_th = r'$D = ' + f'{popt[0]:.2e}' + r'$ , $H = ' + f'{H:.1f}' + r'$'
ax.plot(k_th/2/np.pi,omega_th/2/np.pi,'w--',label = label_th)

ax.legend()
# figname = f'FK_file_{idx_file}_s_{x[idx_x]}m_{Nb_minutes}min_M{M}_hop{hop}_short_range'
# figname = f'{fig_folder}{figname}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')







#%% Fit FK curve by Squire theory 

# compute theory
E = 3e9
rho_w = 1027
nu = 0.3
h = 0.4
D = E*h**3/12/(1-nu**2)

H = 2.0

k = np.linspace(0,1,100)
omega_QS = shallow_hydroelastic(k, D, rho_w, H)

# plot experiments + theory
fig, ax = plt.subplots()
imsh = ax.imshow(np.log10(mean_Sfx[:,:,idx_x]),origin = 'lower',aspect = 'auto',cmap = parula_map, 
          extent = [SFT_extents[2],SFT_extents[3],freq[0],freq[-1]])
imsh.set_clim([1,5])

ax.plot(k/2/np.pi,omega_QS/2/np.pi,'r--')

ax.set_xlim([0,0.1])
ax.set_ylim([0,1])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\log10(S_{fx})$')

ax.set_xlabel(r'$1/\lambda \; \mathrm{(m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_title(r'Average FK for $s = ' +  f'{x[idx_x]:.2f}' + r'\; \mathrm{m}$')














        