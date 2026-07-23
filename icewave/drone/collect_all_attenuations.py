# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 18:08:34 2026

@author: sebas
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

import os
import glob
import pickle
import sys
sys.path.append('C:/Users/sebas/git/')

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df
import icewave.sebastien.set_graphs as set_graphs

#%% Function section 

def get_date_drone_exp(file2load):
    
    chain = file2load.split('\\')
    date = chain[1]
    drone = chain[3]
    exp = chain[5]
    
    return date, drone, exp

#%% Collect data

year = '2025'
if year == '2025':
    disk = 'Backup25'
    
base = df.find_path(disk,year)

component = 'uz'
dimension = 'space'
folder_attenuation = f'Figures_attenuation_{component}/'
model_path = f'{base}**/Drones/**/matData/**/{folder_attenuation}/attenuation_data_{dimension}*.pkl'

filelist = glob.glob(model_path,recursive = True)

#%% Load data 

data = {}

for i, file2load in enumerate(filelist):
    print(file2load)
    with open(file2load,'rb') as pf: 
        m = pickle.load(pf)
    
    date,drone,exp = get_date_drone_exp(file2load)
    suffixe = f'{date}_{drone}_{exp}'
    
    data[suffixe] = m

#%% Plot data 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for key,m in data.items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1e0])
ax.set_ylim([1e-3,5e-1])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
    



