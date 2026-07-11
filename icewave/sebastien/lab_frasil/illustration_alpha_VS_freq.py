# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 10:53:41 2026

@author: sebas
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
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

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 

def get_folder_date(foldername):
    chain = foldername.split('\\')
    date = chain[1][:10]
    return date

def get_boundaries(h_list):
    # define bounds for colormap
    test = [(h_list[i] + h_list[i+1])/2 for i in range(len(h_list)-1)]
    bounds = [h_list[0] - (test[0] - h_list[0])]
    for elem in test:
        bounds.append(elem)
    bounds.append(2*h_list[-1] - test[-1])
    bounds = np.array(bounds)
    return bounds

def collect_all_experiments(filelist):
    data_list = []
    for file in filelist:
        with open(file,'rb') as pf:
            data = pickle.load(pf)
        data_list.append(data)
        
    return data_list

def classify_dataset(data_set,selection_keys,rho,method = 'laser'):
    """ Classify data from a set of experiments
    Inputs : - data_set, list of dictionnaries, each dictionnary for a single experiment
             - selection_keys, list of keys that are kept in the returned dictionnary
             - rho, float, water density used for all experiments
             - method, string, method used to extract all data, default is laser 
    Output : - main_dict, dictionnary"""
             
    main_dict = {}
    for data in data_set : 
        h = data['h']
        f_ex = data['f_ex']
        amp = data['amplitude']
        
        key_exp = f'h{h:.1f}_fex{f_ex:.1f}_amp{amp:.1f}_rho{rho}_{method}'
        print(key_exp)
        
        sub_dict = {}
        sub_dict[key_exp] = {}
        for key in selection_keys:
            current_dict = {key : data[key]}
            sub_dict[key_exp].update(current_dict)

        sub_dict[key_exp]['method'] = method
        sub_dict[key_exp]['rho'] = rho
        
        main_dict.update(sub_dict)
    
    return main_dict

def get_results_folder(folder,format_filename,selection_keys,rho_dict,method = 'laser'):
    """ Get results for all frequencies associated to a given thickness and sort them in a dictionnary """
    
    path2data = folder + format_filename
    filelist = glob.glob(path2data)
    print(filelist)
    date = get_folder_date(folder)
    rho = rho_dict[date] 
    data_set = collect_all_experiments(filelist)  
    
    main_dict = classify_dataset(data_set,selection_keys,rho,method)
    for key in main_dict:
        main_dict[key]['date'] = date
    
    return main_dict

def get_subset_dict(data, selection):
    # On ne garde que les expériences qui valident TOUS les critères de sélection
    return {
        key: exp_data 
        for key, exp_data in data.items()
        if all(exp_data.get(param) == value for param, value in selection.items())
    }

#%% Set fig_folder

fig_folder = 'F:/PhD_Manuscript/ch3/Attenuation/illustration/'
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)

#%% Gather Mariya data

base =  'F:/Stagiaires/Mariya/'
folderlist = glob.glob(f'{base}*e*laser*/')

# create a table of density
rho_dict = {'2026_06_11':1.067,'2026_06_12':1.067,
            '2026_06_17':1.03,'2026_06_18':1.03,'2026_07_06':1.135}

# selected keys
selection_keys = ['h','f_ex','amplitude','alpha','err_alpha','k0','err_k0','f_demod',
                  'lorentz','exponential','Amax']
data = {}
for folder in folderlist:

    main_dict = {}
    format_filename = 'Laser_attenuation/**/laser_attenuation*.pkl'
    main_dict_laser = get_results_folder(folder,format_filename,selection_keys,rho_dict)
    main_dict.update(main_dict_laser)
    
    data.update(main_dict)

file2save = f'{base}attenuation_results.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(data,pf)


# =============================================================================
# %% Load Mariya's data
# =============================================================================

base =  'F:/Stagiaires/Mariya/'
path2data = f'{base}attenuation_results.pkl'

with open(path2data,'rb') as pf:
    data = pickle.load(pf)

# =============================================================================
# %% Test get_subset_dict
# =============================================================================
amp = 25 # in mm
h = 22.0 # in mm

selection = {'amplitude':amp,'h':h}
sub_data = get_subset_dict(data, selection)
print(sub_data.keys())

# =============================================================================
# %% Plot data for a given density and amplitude     
# =============================================================================

amp = 25 # in mm
rho = 1.03 # water density

selection = {'amplitude':amp,'rho':rho}
sub_data = get_subset_dict(data, selection)

h_list = []
for key in sub_data:
    if sub_data[key]['h'] not in h_list:
        h_list.append(sub_data[key]['h'])
h_list = np.array(h_list)

bounds = get_boundaries(h_list)

cmap = new_blues
# cnorm = mcolors.Normalize(h_list.min(),h_list.max())
cnorm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for i,h in enumerate(h_list):
    selection = {'h':h}
    current_dict = get_subset_dict(sub_data, selection)
    
    alpha = np.array([current_dict[key]['alpha'] for key in current_dict.keys()])
    err_alpha = np.array([current_dict[key]['err_alpha'] for key in current_dict.keys()])
    f = np.array([current_dict[key]['f_demod'] for key in current_dict.keys()])
    
    color = cmap(cnorm(h))
    ax.errorbar(f,alpha,yerr = 5*err_alpha,color = color, fmt = 'o',
                ecolor ='k',mec = 'k')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([2e0,7e1])
ax.set_xlim([1.8e0,6e0])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
# midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
# cbar.set_ticks(midpoints)
cbar.set_ticks(h_list)
cbar.set_ticklabels([f'{mid:.1f}' for mid in h_list])
cbar.set_label(r'$h \; \mathrm{(mm)}$')
    
# # create a scalar mappable
# sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
# # sm.set_array([])  # Only needed for the colorbar
# cbar = fig.colorbar(sm, cax=cax)
# cbar.set_label(r'$h \; \mathrm{(mm)}$')

xth = np.linspace(2,3)
yth = 0.8*xth**3.5
ax.plot(xth,yth,'k--')

figname = f'{fig_folder}alpha_VS_freq_rho_{rho}_amp_{amp}mm'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot data for different amplitudes 

selection = {'h':8.0,'rho':1.067}

sub_data = get_subset_dict(data, selection)

amp_list = []
for key in sub_data:
    if sub_data[key]['amplitude'] not in amp_list:
        amp_list.append(sub_data[key]['amplitude'])

A_list = np.array([sub_data[key]['exponential']['A'] for key in sub_data.keys()])
# A_list = np.array([sub_data[key]['Amax'] for key in sub_data.keys()])*1e3

k_list = np.array([sub_data[key]['k0'] for key in sub_data.keys()])
omega_list = np.array([sub_data[key]['f_demod'] for key in sub_data.keys()])*2*np.pi

gamma_list = A_list*k_list
gamma_dot_list = gamma_list*omega_list

cnorm = mcolors.Normalize(vmin = A_list.min()*1e3,vmax = A_list.max()*1e3)
full = mpl.colormaps['Reds'].resampled(256)
new_cmap = mcolors.ListedColormap(full(np.linspace(0.2,0.8,256)))
cmap = new_cmap

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
marker_list = ['o','s','D','^']

for i,amp in enumerate(amp_list):
    
    selection = {'amplitude':amp}
    current_dict = get_subset_dict(sub_data, selection)
    
    alpha = np.array([current_dict[key]['alpha'] for key in current_dict.keys()])
    err_alpha = np.array([current_dict[key]['err_alpha'] for key in current_dict.keys()])
    f = np.array([current_dict[key]['f_demod'] for key in current_dict.keys()])
    A = np.array([current_dict[key]['exponential']['A'] for key in current_dict.keys()])
    # A = np.array([current_dict[key]['Amax'] for key in current_dict.keys()])*1e3
    
    k = np.array([current_dict[key]['k0'] for key in current_dict.keys()])
    gamma = A*k
    gamma_dot = gamma*f*2*np.pi
    ax.scatter(f,alpha,s= 50,c = A*1e3,cmap = cmap,norm = cnorm,
               marker = marker_list[i], edgecolor = 'k')

    
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([3e0,1e2])
ax.set_xlim([2.3e0,6e0])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=cnorm)
# sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label(r'$A \; \mathrm{(mm)}$')

xth = np.linspace(2.8,3.8)
yth = 0.8*xth**3.5
ax.plot(xth,yth,'k--')

figname = f'{fig_folder}alpha_VS_freq_rho_{rho}_h_{h}mm_colorbar_A'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
    
    