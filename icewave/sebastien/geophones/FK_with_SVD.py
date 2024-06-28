#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:04:44 2024

@author: moreaul
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:49:15 2024

@author: moreaul
"""

from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib.patches import PathPatch
from matplotlib import ticker

plt.rcParams['text.usetex'] = True
# Parameters for plots
font_size_medium = 30
font_size_small = 26
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (12,9)
img_quality = 1000 # dpi to save images 



def fn_svd(signals, fs, xs, rang,name, issaving, *varargin):
    if varargin:
        if varargin[0] == 'threshold':
            threshold_user = varargin[1]
        elif varargin[0] == 'rang':
            rang = varargin[1]
        else:
            print('varargin(1) unknown')
            return
    else:
        print('varargin empty')

    Nreceiv, Nt, Nemit = signals.shape # number of recepter, time steps and emitters 

    # time domain fft
    Nf = 2048
    f = (np.arange(Nf) / Nf) * (fs if fs else 1)
    f_axename = 'f/fs' if not fs else 'f'
    SIGNALS = fft(signals, Nf, axis=1)
    SIGNALS = SIGNALS[:, :Nf + 1, :]

    # svd
    U = np.zeros((Nreceiv, Nf, Nemit), dtype=complex)
    S = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    V = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    D = np.zeros((Nemit, Nf), dtype=complex)

    for ii in range(Nf):
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False)
        D[:, ii] = np.diag(S[:, ii, :])

    for ne in range(Nemit):
        titi = 20 * np.log10(np.abs(D[ne, :]) / np.max(np.abs(D[0, :])))
        plt.plot(f, titi, label=f'Slice {ne}')

    if threshold_user is None:
        plt.xlabel('frequency (Hz)')
        plt.ylabel('Singular values (in dB of peak value)')
        [fcut, sigmacutsup] = plt.ginput(5)
        fcut = [f[0]] + fcut.tolist() + [f[-1]]
        sigmacutsup = [sigmacutsup[0]] + sigmacutsup.tolist() + [sigmacutsup[-1]]
        sigmacutsup = np.interp(f, fcut, sigmacutsup)
    else:
        sigmacutsup = np.full(int(Nf/2+1), threshold_user)
        sigmacutsup = threshold_user
        # print(sigmacutsup)              
    # plt.plot(f, sigmacutsup)

    for ne in range(Nemit):
        titi = 20 * np.log10(D[ne, :] / np.max(D[0, :]))
        idx = np.where(titi <= sigmacutsup)[0]
        U[:, idx, ne] = 0

    if issaving:
        plt.savefig(name + '_sv')

    # projection onto each singular vector
    Nk = 2048
    k = (np.arange(Nk) / Nk) * (2 * np.pi / xs)  if xs else np.arange(Nk + 1)
    k_axename = 'k/ks' if not xs else 'k'

    projections = ifft(U, Nk, axis=0)#np.fft.fftshift(ifft(U, Nk, axis=0), axes=0)
    # projections = projections[:Nk + 1, :, :]
    # projections = projections[:, :-1, :]
    projections_sum = np.zeros((Nf, Nk, Nemit))

    # for kemit in range(Nemit):
    #     for ii in range(Nf + 1):
    #         projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]) ** 2


    for kemit in rang:
        for ii in range(Nf):
            max_value = 1  # np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[:, ii, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2

    projections_sum = np.abs(np.mean(projections_sum, axis=2))

    return f, k, projections_sum

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#%%

rang = [0,1,2] # rank of singular values 
geophones_spacing = 3 # spacing between two geophones in meter
signals = np.transpose(seismic_matrix, (0, 2, 1))
#signals = np.transpose(seismic_matrices_interp, (0, 2, 1))

f, k, FK = fn_svd(signals, fs, geophones_spacing,rang,'ExampleName', 0, 'threshold',-60)
# dimension 1 : f
# dimension 2 : k
F, K = np.meshgrid(f, k)
#%%
######################################################
# ! NEEDS FLIPUD DEPENDING ON DIRECTION OF PROPAGATION
######################################################

#%%
FK = FK/np.max(FK)
# FK = np.flipud(FK/np.max(FK))

vmin = 0     # Minimum value for colormap
vmax = 1 # Maximum value for colormap (adjust as needed)
fig, ax1 = plt.subplots(1, 1, figsize = fig_size)

# Computes extents for imshow
x = extents(k)
y = extents(f)

c1 = ax1.imshow(np.transpose(FK), aspect = 'auto', cmap='gnuplot2',origin = 'lower',extent = x + y,vmin = vmin, vmax = vmax)

#ax1.set_ylim([-1.5, 1.5])
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$',labelpad = 5)
# ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}|_{max}}(f,k)$')
ax1.tick_params(axis='both', which='major', pad=7)
ax1.set_ylim([0, 400])

#%%
figname = fig_folder + 'FK_plot_dir1_' + channel_correspondance[channel]
plt.savefig(figname + '.pdf',dpi = img_quality, bbox_inches='tight')
plt.savefig(figname + '.png',dpi = img_quality, bbox_inches='tight')



# #%%
# FK = FK/np.max(FK)
# # FK = np.flipud(FK/np.max(FK))

# vmin = 0     # Minimum value for colormap
# vmax = 1 # Maximum value for colormap (adjust as needed)
# fig, axes = plt.subplots(1, 1, figsize = (16,10))
# ax1 = axes  # No need for indexing when there's only one subplot
# c1 = ax1.contourf(K, F, FK, cmap='gnuplot2', vmin=vmin, vmax=vmax)
# #ax1.set_ylim([-1.5, 1.5])
# ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
# ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
# # ax1.set_title('Spectrum with SVD filter')
# plt.colorbar(c1, ax=ax1, label= r'$\frac{|\hat{s}|}{|\hat{s}_{max}|}(f,k)$')

# ax1.set_ylim([0, 400])

# #%%
# figname = fig_folder + 'FK_plot_dir1_' + channel_correspondance[channel]
# plt.savefig(figname + '.pdf',dpi = 1000)
# plt.savefig(figname + '.png',dpi = 1000)
# #points = plt.ginput(2, timeout=-1)


#%%
# #----------------------- UNWRAPPING SPECTRUM - HORIZONTAL STACKING ------------------
nb_stacking = 2 # number of times we want to stack the FK plot horizontally
idx_stacking = 1
FK_uwp = np.vstack((FK, FK))
while idx_stacking < nb_stacking :
    FK_uwp = np.vstack((FK_uwp, FK))
    idx_stacking += 1
    
k_uwp= np.linspace(0,(nb_stacking + 1)*max(k),FK_uwp.shape[0])


vmin = 0     # Minimum value for colormap
vmax = 1# Maximum value for colormap (adjust as needed)
fig, ax1 = plt.subplots(1, 1, figsize=(10, 10))

c1 = ax1.imshow(np.transpose(FK_uwp), aspect = 'auto', cmap='gnuplot2',
                origin = 'lower',extent = extents(k_uwp) + extents(f),vmin = vmin, vmax = vmax)
#ax1.set_ylim([-1.5, 1.5])
ax1.set_xlabel(r'$k \; \mathrm{(m^{-1})}$')
ax1.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')

ax1.set_ylim([0, 250])
#%% UNWRAPPING SPECTRUM - VERTICAL STACKING
nb_stacking = 2 # number of times we want to stack the FK plot vertically
idx_stacking = 1
FK_uwp = np.vstack((FK, FK))
while idx_stacking < nb_stacking :
    FK_uwp = np.vstack((FK_uwp, FK))
    idx_stacking += 1
    
k_uwp= np.linspace(0,(nb_stacking + 1)*max(k),FK_uwp.shape[0])
F, K_uwp = np.meshgrid(f,k_uwp )
vmin = 0     # Minimum value for colormap
vmax = 1# Maximum value for colormap (adjust as needed)
fig, axes = plt.subplots(1, 1, figsize=(10, 10))
ax1 = axes  # No need for indexing when there's only one subplot
c1 = ax1.contourf(F, K_uwp, FK_uwp, cmap='gnuplot2', vmin=vmin, vmax=vmax)
#ax1.set_ylim([-1.5, 1.5])
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Wavenumber (k)')
ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')


# ####################################
#%% Only for Flexural mode (mode Z) Select points on the FK plot 
# ####################################
points = plt.ginput(10, timeout=-1)

# Extract x and y coordinates of the selected points
f_mode, k_mode = zip(*points)
# Plot a line between the two selected points
plt.plot(f_mode, k_mode, linestyle='--', color='r', label='Line between points')

file2save = path2data +'/' +acqu_numb+'/'+'dispersion_QS_dir1.pkl'
with open(file2save,'wb') as f:
    pickle.dump([f_mode, k_mode], f)



    
    
    