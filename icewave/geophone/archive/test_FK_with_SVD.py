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

    Nreceiv, Nt, Nemit = signals.shape

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
    projections.shape
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










rang = [0,1]
geophones_spacing = 3
signals = np.transpose(seismic_matrix, (0, 2, 1))
f, k, FK = fn_svd(signals, fs, geophones_spacing,rang,'ExampleName', 0, 'threshold',-60)
F, K = np.meshgrid(f, k)
# ! NEEDS FLIPUD DEPENDING ON DIRECTION OF PROPAGATION
FK = FK/np.max(FK)
#FK = np.flipud(FK/np.max(FK))

vmin = 0     # Minimum value for colormap
vmax = 1# Maximum value for colormap (adjust as needed)
fig, axes = plt.subplots(1, 1, figsize=(10, 10))
ax1 = axes  # No need for indexing when there's only one subplot
c1 = ax1.contourf(F, K, FK, cmap='gnuplot2', vmin=vmin, vmax=vmax)
#ax1.set_ylim([-1.5, 1.5])
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Wavenumber (k)')
ax1.set_title('Spectrum with SVD filter')
plt.colorbar(c1, ax=ax1, label='Spectrum Magnitude')



points = plt.ginput(10, timeout=-1)
# Extract x and y coordinates of the selected points
f_QS, k_QS = zip(*points)
# Plot a line between the two selected points
plt.plot(f_QS, k_QS, linestyle='--', color='r', label='Line between points')

file2save = path2data +'/' +acqu_numb+'/'+'dispersion_QS.pkl'
with open(,'wb') as f:
    pickle.dump([f_QS, k_QS], f)



    
    
    