# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 09:28:28 2025

@author: sebas
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os 
import glob 
import pickle
import h5py
import icewave.tools.matlab2python as mat2py

# import seb plotting package 
import icewave.sebastien.set_graphs as set_graphs

#%% 

def powerfun(x,l):
    return l[1] * x**l[0]

#%% Set plotting parameters 
 
font_size_medium = 30
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

#%%

dates = ['0223','0226']

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Attenuation_BicWin2024/'


#%% Load results for a given day 

data = {}

for date in dates :
    path2data = 'K:/Share_hublot/Data/' + date + '/Drones/'
    # search all .mat file 
    filelist = glob.glob(path2data + '**/attenuation*.pkl',recursive = True)

    # build a dictionnary of results 
    data[date] = {}
    
    for i in range(len(filelist)):
    
        current_file = filelist[i]
        
        with open(current_file, 'rb') as pf:
            current_dict = pickle.load(pf)
            
            print('Top-level keys : ', list(current_dict.keys()))
            
        drone_ID = current_dict['drone_ID']
        if drone_ID == 'Bernache':
            drone_ID = 'bernache'
            current_dict['drone_ID'] = 'bernache'
            
        exp_ID = current_dict['exp_ID']
        
        if drone_ID not in data[date].keys(): # create key if drone_ID does not already exists
            data[date][drone_ID] = {}
            
        data[date][drone_ID][exp_ID] = current_dict # store current_dictionnary

#%% Load results from .mat files

data = {}

for date in dates :
    path2data = 'K:/Share_hublot/Data/' + date + '/Drones/'
    # search all .mat file 
    filelist = glob.glob(path2data + '**/Drone_PIV_main_results*.mat',recursive = True)

    # build a dictionnary of results 
    data[date] = {}
    
    for i in range(len(filelist)):
    
        current_file = filelist[i]
        
        with h5py.File(current_file, 'r') as fmat:
            current_dict = {}
            
            print('Top-level keys : ', list(fmat.keys()))
        
            current_dict = mat2py.mat_to_dict(fmat['main_results'],fmat['main_results'])
            
        drone_ID = current_dict['drone_ID']
        if drone_ID == 'Bernache':
            drone_ID = 'bernache'
            current_dict['drone_ID'] = 'bernache'
            
        exp_ID = current_dict['exp_ID']
        
        if drone_ID not in data[date].keys(): # create key if drone_ID does not already exists
            data[date][drone_ID] = {}
            
        data[date][drone_ID][exp_ID] = current_dict # store current_dictionnary
            
    
#%% Create a colormap and a norm base on drone height
maps = ['Blues','Oranges']
maps_date = {}

for i,date in enumerate(dates):
    maps_date[date] = {}
    
    h_drone = []
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
            h = data[date][key_drone][exp_ID]['DRONE']['h_drone']
            h_drone.append(h)

    h_drone = np.array(h_drone)
    
    full_map = mpl.colormaps[maps[i]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    cnorm = mcolors.Normalize(vmin = min(h_drone), vmax = max(h_drone))
    
    maps_date[date]['new_map'] = new_map
    maps_date[date]['cnorm'] = cnorm
    
# Create markers for each drone 

markers = {'bernache':'o','mesange':'^'}
    
#%% Create a colormap and a norm based on slope of each curve 
maps = ['Blues','Oranges']
maps_date = {}

for i,date in enumerate(dates):
    maps_date[date] = {}
    
    slope = []
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
            current_att = data[date][key_drone][exp_ID]['attenuation']
            current_slope = abs(current_att['power_law']['beta'])
            
            slope.append(current_slope)

    slope = np.array(slope)
    
    full_map = mpl.colormaps[maps[i]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    cnorm = mcolors.Normalize(vmin = min(slope), vmax = max(slope))
    
    maps_date[date]['new_map'] = new_map
    maps_date[date]['cnorm'] = cnorm

#%% superpose TF_spectrum 

fig, ax = plt.subplots()

set_graphs.set_matplotlib_param('single')
plt.rc('legend', fontsize= 12)    # legend fontsize

for date in dates :
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
            x = data[date][key_drone][exp_ID]['FFT_spectrum']['f']
            y = data[date][key_drone][exp_ID]['FFT_spectrum']['TF_spectrum']
            
            ax.loglog(x,y,'-',label = exp_ID)

ax.set_ylim([1e-4,1])
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{V} \rangle _{x,y} (f) $ (u.a.)')

ax.legend(loc = 'lower left')

figname = fig_folder + 'Superposition_TF_spectrum'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% superpose dispersion relation 
set_graphs.set_matplotlib_param('single')
plt.rc('legend', fontsize= 12)    # legend fontsize

fig, ax = plt.subplots()
for date in dates:
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
            h = data[date][key_drone][exp_ID]['DRONE']['h_drone']
            x = data[date][key_drone][exp_ID]['disp_relation']['k']
            y = data[date][key_drone][exp_ID]['disp_relation']['omega']
            
            current_color = maps_date[date]['new_map'](maps_date[date]['cnorm'](h))
            
            current_label = r'$h_w = {x:.2f}$ m / {ID}'.format(x = data[date][key_drone][exp_ID]
                                                        ['disp_relation']['param']['h_w'],ID = exp_ID)
            
            ax.loglog(x,y,markers[key_drone],color = current_color,label = current_label)
    
    
ax.set_xlim([3e-2 , 5])
ax.set_ylim([3e-1 , 10])
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')
ax.legend()

figname = fig_folder + 'Superposition_dispersion_relation'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% superpose attenuation laws
set_graphs.set_matplotlib_param('single')
plt.rc('legend', fontsize= 12)    # legend fontsize

fig, ax = plt.subplots()

f_fit = np.linspace(1e-1,1,100)

d_threshold = 0.05
for date in dates :
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
        
            h = data[date][key_drone][exp_ID]['DRONE']['h_drone']
            current_att = data[date][key_drone][exp_ID]['attenuation']
            # mask = np.logical_and(current_att['d'] < d_threshold, abs(current_att['alpha']) > 1e-3)
            mask = np.where(current_att['mask'])
            x = current_att['f'][mask]
            y = abs(current_att['alpha'][mask])
            
            coeff = abs(current_att['power_law']['beta'])
            current_label = r'$\alpha = {x:.2f}$'.format(x = coeff)
            
            current_color = maps_date[date]['new_map'](maps_date[date]['cnorm'](h))
            ax.loglog(x,y,markers[key_drone],color = current_color,
                      label = current_label)
            power_fit = powerfun(f_fit,[current_att['power_law']['beta'] , current_att['power_law']['B']  ])
            ax.loglog(f_fit,power_fit,'--',color = current_color)
        
ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')

ax.set_xlim([1e-1 , 1])
ax.set_ylim([1e-3 , 7e-1])

figname = fig_folder + 'Superposition_attenuation_law'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% superpose attenuation laws using exponant as colormap

# set_graphs.set_matplotlib_param('single')
# plt.rc('legend', fontsize= 12)    # legend fontsize

fig, ax = plt.subplots(figsize = fig_size)

f_fit = np.linspace(1e-1,1,100)
keywords_list = ['0223_bernache_12-waves_010','0226_mesange_10-waves_005']
colors = ['tab:blue','tab:orange']

d_threshold = 0.05
count = 0
for date in dates :
    for key_drone in data[date].keys():
        for exp_ID in data[date][key_drone].keys():
            
            key_word = f'{date}_{key_drone}_{exp_ID}'
            print(key_word)
            
            if key_word in keywords_list:

                current_att = data[date][key_drone][exp_ID]['attenuation']
                # mask = np.logical_and(current_att['d'] < d_threshold, abs(current_att['alpha']) > 1e-3)
                mask = np.where(current_att['mask'])
                x = current_att['f'][mask]
                y = abs(current_att['alpha'][mask])
                
                coeff = abs(current_att['power_law']['beta'])
                prefactor = current_att['power_law']['B'] 
                print(f'{coeff:.2f}')
                # current_label = r'$\beta = {x:.2f}$'.format(x = coeff)
                current_label = r'$\alpha (f) = ' + f'{prefactor:.1f}' + r'f^{' + f'{coeff:.1f}' + r'}$'
                
                current_slope = current_att['power_law']['beta']
                # current_color = maps_date[date]['new_map'](maps_date[date]['cnorm'](current_slope))
                current_color = colors[count]

                ax.loglog(x,y,'o',color = current_color,markeredgecolor = 'k',markersize = marker_size_plot,label = current_label)
                power_fit = powerfun(f_fit,[current_att['power_law']['beta'] , current_att['power_law']['B']  ])
                ax.loglog(f_fit,power_fit,'--',color = current_color)
                
                count = count + 1
        
ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')

ax.set_xlim([1e-1 , 1])
ax.set_ylim([1e-3 , 7e-1])

figname = fig_folder + 'Superposition_attenuation_law_cmap_slope_selection'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Select only a few experiments




