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

def collect_yearly_attenuation(year,component = 'ux',dimension = 'space',base = None):
    """ Collect all attenuation results for a given year, return a dictionnary. 
    Inputs: 
        - year, str type
        - component, 'ux' or 'uz', corresponds to the wave velocity component used to compute attenuation
        - dimension, 'space' or 'time', corresponds to the dimension over which
        attenuation is fitted using a Lorentzian fit 
        - base, default None, or can be the path to a specific Data base
    Output: 
        - data, dictionnary containing attenuation results """
    
    if base == None:
        if year == '2025':
            disk = 'Backup25'
            base = df.find_path(disk,year)
        elif year == '2024':
            base = 'F:/Rimouski_2024/Data/'
        
    folder_attenuation = f'Figures_attenuation_{component}/'
    model_path = f'{base}**/Drones/**/matData/**/{folder_attenuation}/attenuation_data_{dimension}*.pkl'

    filelist = glob.glob(model_path,recursive = True)
    
    data = {}

    for i, file2load in enumerate(filelist):
        print(file2load)
        with open(file2load,'rb') as pf: 
            m = pickle.load(pf)
        
        date,drone,exp = get_date_drone_exp(file2load)
        suffixe = f'{year}_{date}_{drone}_{exp}'
        m['year'] = year
        m['date'] = date
        m['drone_ID'] = drone    
        m['exp_ID'] = exp
        
        data[suffixe] = m
    
    return data

def collect_all_years(component = 'ux', dimension = 'space'):
    
    years = ['2024','2025']
    data = {}
    for year in years:
        current_data = collect_yearly_attenuation(year,component = component,dimension = dimension)
        data.update(current_data)
        
    return data

#%% Collect data for a given year

year = '2025'
component = 'ux'
dimension = 'time'

data = collect_yearly_attenuation(year,component = component,dimension = dimension)
data_1 = collect_yearly_attenuation('2024',component = component,dimension = dimension)
data.update(data_1)

#%% Plot data 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for key,m in data.items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)

ax.legend(fontsize = 10)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1e0])
ax.set_ylim([1e-3,5e-1])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
    


#%% Collect attenuation from all components and years 

data = {}
components = ['uz','ux']
dims = ['space','time']

for component in components:
    for dim in dims:
        
        key_process = f'{component}_{dim}'
        data[key_process] = collect_all_years(component = component,dimension = dim)
        

# save data 
path2save = 'F:/PhD_Manuscript/ch3/Attenuation/'
file2save = f'{path2save}attenuation_field_main_data.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(data,pf)
    
#%% load data

path = 'F:/PhD_Manuscript/ch3/Attenuation/'
file2load = f'{path}attenuation_field_main_data.pkl'

with open(file2load,'rb') as pf:
    data = pickle.load(pf)
 
#%% Plot data

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for key,m in data[key_process].items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)

ax.legend(fontsize = 10)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1e0])
ax.set_ylim([1e-3,5e-1])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')




















# =============================================================================
# %% Further developments / tests
# =============================================================================



#%%
year = '2025'
if year == '2025':
    disk = 'Backup25'
    base = df.find_path(disk,year)
elif year == '2024':
    base = 'F:/Rimouski_2024/Data/'
    

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
    suffixe = f'{year}_{date}_{drone}_{exp}'
    m['year'] = year
    m['date'] = date
    m['drone_ID'] = drone    
    m['exp_ID'] = exp
    
    data[suffixe] = m
    


