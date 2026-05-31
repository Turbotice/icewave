# -*- coding: utf-8 -*-
"""
Created on Sat May 30 16:21:34 2026

@author: sebas

Plot histogramm obtained from Matlab ice floes detection 
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import h5py

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Load data

base ='F:/Amundsen_RA_2024/Data/'
date = '0921'
drone_ID = 'bernache'
exp_ID = '06-waves_001'
suffixe = f'{date}_{drone_ID}_{exp_ID}'

path2data = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/Figures/'
filelist = glob.glob(f'{path2data}*Data*.mat')

fig_folder = path2data
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

list_dict = []
for file2load in filelist:
    print(file2load)
    
    # load file 
    with h5py.File(file2load, 'r') as fmat:
        S = {}
    
        print('Top-level keys : ', list(fmat.keys()))
    
        S = mat2py.mat_to_dict(fmat['results'],fmat['results'])
    list_dict.append(S)

    
#%% Plot Histogramm

# wave field characteristics
wavelength = 125 
R_eff = wavelength/4
A_eff = np.pi*R_eff**2

# Possible flexural length 
E = 4.8e9
nu = 0.3
rho_w = 1000
h_ice = 4.5
D = E*h_ice**3/(1-nu**2)/12
Ld = (D/rho_w)**0.25
print(Ld)

for i,file2load in enumerate(filelist):
    S = list_dict[i]
    
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots()
    hist = ax.hist(S['stats']['area_real'].T, bins=20,range = (300,8000),
            density=True, alpha=0.5,color='skyblue', edgecolor='black', label='Histogramme')
    ax.axvline(A_eff, color='red', linestyle='--', linewidth=2)
    # ax.axvline(np.pi*(Ld/2)**2, color='tab:blue', linestyle='--', linewidth=2)
    ax.set_xlabel(r'Area $\mathrm{(m^2)}$')
    ax.set_ylabel(r'PDF')
    
    string_list = file2load.split('Figures')
    img_name = string_list[-1][1:7]
    print(img_name)
    
    figname = f'{fig_folder}FSD_histogram_{suffixe}_{img_name}'
    # plt.savefig(figname + '.pdf', bbox_inches='tight')
    # plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute ice thickness assuming equality between wavelength and Ld

