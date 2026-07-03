# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 17:45:35 2026

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy

import os 
import glob  
import pickle
import cmocean

import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% FUNCTION SECTION 

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

def exponential_decay(x,A,alpha):
    y = A*np.exp(-alpha*x)
    return y

def lorentzian_fit(cut,k,k0_idx,bounds_alpha):
    """ Perform Lorentzian fit along a profile. 
    Inputs : - cut, 1D numpy array, profile along which lorentzian fit must be performed
             - k, 1D numpy array or list, list of wavevectors
             - k0_idx, index of initial guess for center of lorentzian
             - bounds_alpha, boundaries between which lorentzian width must be checked """
    
    y_exp = (cut - cut.min())/(cut.max() - cut.min())

    bounds_kx0 = (k[k0_idx - 10],k[k0_idx + 10])
    bounds_curvefit = ([bounds_kx0[0],bounds_alpha[0]],[bounds_kx0[1],bounds_alpha[1]])
    # fit by a lorentzian
    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),k,y_exp,
                                         bounds = bounds_curvefit)
    err_coeff = np.sqrt(np.diag(pcov))
    
    return popt,err_coeff

def exponential_fit(y,x,bounds):
    """ Perform exponential fit of data """
    popt,pcov = scipy.optimize.curve_fit(lambda x,A,alpha : exponential_decay(x,A,alpha),x,y,bounds = bounds)
    err_coeff = np.sqrt(np.diag(pcov))

    return popt,err_coeff

def FFT2D_positive_freq(spatio,facq):
    """ Compute FFT2 of spatio-temporal signal, and keep only part with positive frequency
    """
    shift,k,omega = FT.fft_2D(spatio,facq)
    freq = omega/2/np.pi
    
    positive_shift = shift[:,shift.shape[1]//2:]
    positive_freq = freq[len(freq)//2:]
    
    return positive_shift,k,positive_freq


#%% Set fig_folder

fig_folder = 'F:/PhD_Manuscript/ch3/Attenuation/illustration/'
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)

#%% Load data 

base =  'F:/Stagiaires/Mariya/'
date = '2026_06_17'
h = 7.0 # layer height 
a = 15 # motor amplitude

path2data = f'{base}{date}_e_{h}mm_laser_a{a}mm/Laser_extraction/'
filelist = glob.glob(f'{path2data}scaled*.pkl')
print(filelist)

i = 6
file2load = filelist[i]
with open(file2load,'rb') as pf:
    data = pickle.load(pf)
    
    
h = data['h']
f_ex = data['f_ex']
a = data['amplitude']

suffixe = f'h_{h}_fex_{f_ex}_a_{a}'

#%% Plot spatio-temporal 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(data['spatio'].T*1e3, origin = 'lower', 
          extent = [data['x'][0], data['x'][-1], data['t'][0], data['t'][-1]],
          cmap = cmocean.cm.balance,aspect = 'auto')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\xi \; \mathrm{(mm)}$')
imsh.set_clim([-4,4])

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$t \; \mathrm{(s)}$')

ax.set_ylim([0,6])

#%% Plot FK-spectrum 

# Compute FFT2D
facq = [data['SCALE']['facq_pix'],data['SCALE']['facq_t']]
shift,k,freq = FFT2D_positive_freq(data['spatio'],facq)

# keep frequencies higher than threshold
threshold_freq = 0.5
mask = freq > threshold_freq
positive_freq = freq[mask]
positive_shift = shift[:,mask]

# find maximum
idx_max = np.argmax(abs(positive_shift).flatten())
unravel_coords = np.unravel_index(idx_max,positive_shift.shape)
f_demod = positive_freq[unravel_coords[1]]
k0 = - k[unravel_coords[0]]

# Fit cut for a given frequency by a lorentzian
cut = abs(positive_shift)[:,unravel_coords[1]]
bounds_alpha = (0,1e2)
popt,err_coeff = lorentzian_fit(cut,k,unravel_coords[0],bounds_alpha)

label_fit = r'$\alpha = ' + f'{popt[1]:.1f}' + '\; \mathrm{(m^{-1})}$'
print(label_fit)
k_fit = np.linspace(k.min(),k.max(),4000)
yth = lorentzian(k_fit,popt[0],popt[1])*(cut.max() - cut.min()) + cut.min()

# Plot 
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(abs(positive_shift).T*1e3,cmap = 'gray_r',aspect = 'auto',norm = 'linear',origin = 'lower',
               extent = (-k[0],-k[-1],positive_freq.min(),positive_freq.max()))

ax.set_ylim([0,5])
ax.set_xlim([-150,150])
# ax.plot(k[unravel_coords[0]],f_demod,'ro')
ax.axvline(k0,ymin = 0, ymax = 1,ls = '--', color = 'r', lw = 2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\hat{\xi} (k,f) \; \mathrm{(mm)}$',labelpad = 1)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)

# Show lorentzian fit 
axins = inset_axes(ax,width = '40%',height = '40%', loc = 'lower left',
                   bbox_to_anchor=(0.1, 0.1, 1, 1),
                   bbox_transform=ax.transAxes)
axins.plot(-k,cut*1e3,'o-')
axins.plot(-k_fit,yth*1e3,'r',label = label_fit)
# axins.legend()

# axins.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
# axins.set_ylabel(r'$|\hat{\xi}|(k) \; \mathrm{(mm)}$')
axins.set_xlim([k0 -150, k0 +150])


figname = f'{fig_folder}FK_insert_lorentzian_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


# =============================================================================
# %% Subplot : spatio + FK with insert
# =============================================================================

set_graphs.set_matplotlib_param('double')
fig, axs = plt.subplots(nrows = 1, ncols = 2,layout = 'constrained',figsize = (14,6))

# plot spatio
ax = axs[0]
imsh = ax.imshow(data['spatio'].T*1e3, origin = 'lower', 
          extent = [data['x'][0], data['x'][-1], data['t'][0], data['t'][-1]],
          cmap = cmocean.cm.balance,aspect = 'auto')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\xi \; \mathrm{(mm)}$')
imsh.set_clim([-4,4])

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$t \; \mathrm{(s)}$')
ax.set_ylim([0,6])

# plot FK with insert 
ax = axs[1]
imsh = ax.imshow(abs(positive_shift).T*1e3,cmap = 'gray_r',aspect = 'auto',norm = 'linear',origin = 'lower',
               extent = (-k[0],-k[-1],positive_freq.min(),positive_freq.max()))

ax.set_ylim([0,5])
ax.set_xlim([-150,150])
# ax.plot(k[unravel_coords[0]],f_demod,'ro')
ax.axvline(k0,ymin = 0, ymax = 1,ls = '--', color = 'r', lw = 2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\hat{\xi} (k,f) \; \mathrm{(mm)}$',labelpad = 1)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)

# Show lorentzian fit 
axins = inset_axes(ax,width = '40%',height = '40%', loc = 'lower left',
                   bbox_to_anchor=(0.1, 0.1, 1, 1),
                   bbox_transform=ax.transAxes)
axins.plot(-k,cut*1e3,'o-')
axins.plot(-k_fit,yth*1e3,'r',label = label_fit)
# axins.legend()

# axins.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
# axins.set_ylabel(r'$|\hat{\xi}|(k) \; \mathrm{(mm)}$')
axins.set_xlim([k0 -150, k0 +150])


figname = f'{fig_folder}subplot_spatio_FK_insert_lorentzian_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
# %% Demodulate profile and fit by an exponential 
# =============================================================================
demod_profile = np.mean(data['spatio']*np.exp(1j*2*np.pi*f_demod*data['t']),axis = -1)

bounds = ([1e-6,1e-4],[1e-1,1e2])
x = data['x']
mask = np.where(x > 0.015)[0]
popt_exp,err_coeff_exp = exponential_fit(abs(demod_profile)[mask], x[mask], bounds)

xth = np.linspace(x.min(),x.max(),200)
yth = exponential_decay(xth,popt_exp[0],popt_exp[1])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(x,np.real(demod_profile))
# ax.plot(x,abs(demod_profile))
ax.plot(xth,yth,'r')

ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$\hat{\xi}(x) \; \mathrm{(mm)}$',labelpad = 5)
ax.set_ylim([-2e-3,3.5e-3])

axins = inset_axes(ax, width="40%", height="40%")
axins.plot(x,abs(demod_profile),color = 'k')
axins.plot(xth,yth,'r')
axins.set_yscale('log')

axins.set_xlabel(r'$x$')
axins.set_ylabel(r'$|\hat{\xi}|$')

figname = f'{fig_folder}Demod_profile_insert_log_amp_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')



