# -*- coding: utf-8 -*-
"""
Created on Tue May  6 10:59:19 2025

@author: sebas
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

import pickle
import os
import glob
import re
import cv2 as cv 

from concurrent.futures import ProcessPoolExecutor

import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))
#%% FUNCTION SECTION 
def collect_all_experiments(filelist):
    data_list = []
    for file in filelist:
        with open(file,'rb') as pf:
            data = pickle.load(pf)
        data_list.append(data)
        
    return data_list

def classify_dataset(data_set,selection_keys,method = 'laser'):
    """ Classify data from a set of experiments
    Inputs : - data_set, list of dictionnaries, each dictionnary for a single experiment
             - selection_keys, list of keys that are kept in the returned dictionnary
             - method, string, method used to extract all data, default is laser 
    Output : - main_dict, dictionnary"""
             
    main_dict = {}
    for data in data_set : 
        h = data['h']
        f_ex = data['f_ex']
        amp = data['amplitude']
        
        key_exp = f'h{h:.1f}_fex{f_ex:.1f}_amp{amp:.1f}_{method}'
        print(key_exp)
        
        sub_dict = {}
        sub_dict[key_exp] = {}
        for key in selection_keys:
            current_dict = {key : data[key]}
            sub_dict[key_exp].update(current_dict)

        sub_dict[key_exp]['method'] = method
        
        main_dict.update(sub_dict)
    
    return main_dict

def affine(x,a,b):
    y = a*x + b
    return y

def powerlaw_fit(x,y,err_y = None):
    """ Fit data using a power law, taking into account standard deviation of y """
    log_x = np.log(x)
    log_y = np.log(y)
    
    if err_y is None :
        popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y)
    
    else : 
        
        err_log_y = err_y/y
        popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y,sigma = err_log_y,
                                             absolute_sigma = True)
        
    err_affine = np.sqrt(np.diag(pcov))
    beta = popt[0]
    err_beta = err_affine[0]
    B = np.exp(popt[1])
    err_B = B*err_affine[1]
    
    coeffs = (beta,B)
    err_coeffs = (err_beta,err_B)
    return coeffs,err_coeffs


def get_results_folder(folder,format_filename,selection_keys,method = 'laser'):
    """ Get results for all frequencies associated to a given thickness and sort them in a dictionnary """
    
    path2data = folder + format_filename
    filelist = glob.glob(path2data)
    print(filelist)
    data_set = collect_all_experiments(filelist)  
    
    main_dict = classify_dataset(data_set,selection_keys,method)
    return main_dict

#%% Load data 
h = 5.0
date = '2024_07_10'
main_path = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/'

# data from laser 
path2data_laser = f'{main_path}Laser_attenuation/'
filelist = glob.glob(f'{path2data_laser}**/laser_attenuation*.pkl')
# data from PIV
path2data_PIV = f'{main_path}matData/'
filelist_PIV = glob.glob(f'{path2data_PIV}**/PIV_attenuation_results_h*.pkl')

data_laser = collect_all_experiments(filelist)
data_PIV = collect_all_experiments(filelist_PIV)


#%% Create a single dictionnary

main_dict = {}
main_dict['laser'] = {}
keys = ['alpha','err_alpha','k0','f_demod','f_ex']
for key in keys:
    main_dict['laser'][key] = np.array([data[key] for data in data_laser ])

main_dict['PIV'] = {}
keys = ['alpha','err_alpha','k0','f_demod','f_ex']
for key in keys:
    main_dict['PIV'][key] = np.array([data[key] for data in data_PIV ])

#%% Plot data 

power_laws = {}
coeffs,err_coeffs = powerlaw_fit(main_dict['laser']['f_demod'], main_dict['laser']['alpha'], main_dict['laser']['err_alpha'])
power_laws['laser'] = {}
power_laws['laser']['coeffs'] = coeffs
power_laws['laser']['err_coeffs'] = err_coeffs
coeffs,err_coeffs = powerlaw_fit(main_dict['PIV']['f_demod'], main_dict['PIV']['alpha'], main_dict['PIV']['err_alpha'])
power_laws['PIV'] = {}
power_laws['PIV']['coeffs'] = coeffs
power_laws['PIV']['err_coeffs'] = err_coeffs

yth = {}
xth = np.linspace(1e0,1e1,200)
for key in ['laser','PIV']:
    yth[key] = power_laws[key]['coeffs'][1]*xth**power_laws[key]['coeffs'][0]


fig, ax = plt.subplots(figsize = (12,9))
ax.plot(main_dict['laser']['f_demod'],main_dict['laser']['alpha'],'o',color = 'tab:blue',label = 'laser')
coeffs = power_laws['laser']['coeffs']
label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
ax.plot(xth,yth['laser'],'-',color = 'tab:blue',label = label_th)
ax.plot(main_dict['PIV']['f_demod'],main_dict['PIV']['alpha'],'d',color = 'tab:orange',label = 'PIV')
coeffs = power_laws['PIV']['coeffs']
label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
ax.plot(xth,yth['PIV'],'-',color = 'tab:orange',label = label_th)
ax.legend()

ax.set_xlabel(r'$f_{\mathrm{demod}} \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2e0,8e0])
ax.set_ylim([1e0,1e2])

figname = f'{main_path}Attenuation_law_PIV_laser_h{h:.2f}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')





#############################################
#%% Create a dictionnary more "horizontal"
#############################################

folderlist = glob.glob('U:/Aurore_frasil/*e_*mm_laser/')
selection_keys = ['h','f_ex','amplitude','alpha','err_alpha','k0','err_k0','f_demod',
                  'lorentz','exponential','Amax']
data = {}
for folder in folderlist:
    
    main_dict = {}
    format_filename = 'Laser_attenuation/**/laser_attenuation*.pkl'
    main_dict_laser = get_results_folder(folder,format_filename,selection_keys)
    main_dict.update(main_dict_laser)
    
    format_filename = 'matData//**/PIV_attenuation_results_h*.pkl'
    main_dict_PIV = get_results_folder(folder,format_filename,selection_keys,method = 'PIV')
    main_dict.update(main_dict_PIV)
    
    data.update(main_dict)
    
file2save = 'U:/Aurore_frasil/attenuation_results.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(data,pf)

#%% Create folder where figures are saved 
fig_folder = 'U:/Aurore_frasil/Results_Seb/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)    

#%% Plot attenuation coeff for a single method and all thicknesses

set_graphs.set_matplotlib_param('single')
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots(figsize = (12,9))


for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = np.array([data[key]['alpha'] for key in valid_keys])
    err_alpha = np.array([data[key]['err_alpha'] for key in valid_keys])
    f_demod = np.array([data[key]['f_demod'] for key in valid_keys])
    
    
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(f_demod,alpha,yerr = err_alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
    xth = np.linspace(2.0,6.0,100)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha, err_alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    ax.plot(xth,yth,'-',color = current_color,label = label_th)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([2.7,5.3])
ax.set_ylim([2e0,70])


figname = f'{fig_folder}attenuation_all_thicknesses_logscale_{method}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Compare attenuation coeff 

set_graphs.set_matplotlib_param('single')
method = 'PIV'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots()

for i,h in enumerate(thickness_list):
    valid_keys = {'PIV':[],'laser':[]}
    for key in data.keys():
        if data[key]['h'] == h and data[key]['method'] == 'PIV':
            valid_keys['PIV'].append(key)
        if data[key]['h'] == h and data[key]['method'] == 'laser':
            valid_keys['laser'].append(key)
    
    alpha_PIV = [data[key]['alpha'] for key in valid_keys['PIV']]
    alpha_laser = [data[key]['alpha'] for key in valid_keys['laser']]
        
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(alpha_laser,alpha_PIV,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')


ax.set_ylabel(r'$\alpha_{PIV} \; \mathrm{(m^{-1}})$',labelpad = 5)
ax.set_xlabel(r'$\alpha_{laser} \; \mathrm{(m^{-1}})$',labelpad = 5)
limits = [1e0,60]
yth = np.linspace(limits[0],limits[1],200)
ax.plot(limits,limits,'r--')
ax.set_xlim([1e0,60])
ax.set_ylim([1e0,60])
# ax.set_xscale('log')
# ax.set_yscale('log')


figname = f'{fig_folder}Comparison_attenuation_coeff_PIV_VS_laser'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')


#%% Plot A_max VS frequency 

set_graphs.set_matplotlib_param('single')
method = 'PIV'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots()


for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    Amax = [data[key]['Amax'] for key in valid_keys]
    f_demod = [data[key]['f_demod'] for key in valid_keys]
    
    
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(f_demod,Amax,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')


ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$A_{max} \; \mathrm{(u.a.)}$',labelpad = 5)
scalebase = 'linear'
ax.set_xscale(scalebase)
ax.set_yscale(scalebase)
ax.set_title(f'{method}')
# ax.set_xlim([2.7,5.3])
# ax.set_ylim([2e0,70])


figname = f'{fig_folder}Amax_VS_freq_method_{method}_scale_{scalebase}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Plot A_fit VS frequency

set_graphs.set_matplotlib_param('single')
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots()


for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    Afit = [data[key]['exponential']['A'] for key in valid_keys]
    f_demod = [data[key]['f_demod'] for key in valid_keys]
    
    
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(f_demod,Afit,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')


ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$A_{fit} \; \mathrm{(u.a.)}$',labelpad = 5)
scalebase = 'linear'
ax.set_xscale(scalebase)
ax.set_yscale(scalebase)
ax.set_title(f'{method}')
# ax.set_xlim([2.7,5.3])
# ax.set_ylim([2e0,70])


figname = f'{fig_folder}Afit_VS_freq_method_{method}_scale_{scalebase}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Compare A_max between PIV and laser 

set_graphs.set_matplotlib_param('single')
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots()

for i,h in enumerate(thickness_list):
    valid_keys = {'PIV':[],'laser':[]}
    for key in data.keys():
        if data[key]['h'] == h and data[key]['method'] == 'PIV':
            valid_keys['PIV'].append(key)
        if data[key]['h'] == h and data[key]['method'] == 'laser':
            valid_keys['laser'].append(key)
    
    Amax_PIV = [data[key]['Amax'] for key in valid_keys['PIV']]
    Amax_laser = [data[key]['Amax'] for key in valid_keys['laser']]
        
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(Amax_laser,Amax_PIV,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')


ax.set_ylabel(r'$A_{PIV} \; \mathrm{(u.a.)}$',labelpad = 5)
ax.set_xlabel(r'$A_{laser} \; \mathrm{(m)}$',labelpad = 5)
scalebase = 'linear'
ax.set_xscale(scalebase)
ax.set_yscale(scalebase)
# limits = [1,60]
# yth = np.linspace(limits[0],limits[1],200)
# ax.plot(limits,limits,'r--')
# ax.set_xlim([1e0,60])
# ax.set_ylim([1e0,60])


figname = f'{fig_folder}Comparison_Amax_PIV_VS_laser'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Plot exponential alpha VS frequency 

set_graphs.set_matplotlib_param('single')
method = 'PIV'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots(figsize = (12,9))


for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    alpha = [data[key]['exponential']['alpha'] for key in valid_keys]
    # err_alpha = [data[key]['exponential']['err_alpha'] for key in valid_keys]
    f_demod = [data[key]['f_demod'] for key in valid_keys]
    
    
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(f_demod,alpha,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)
    coeffs,err_coeff = powerlaw_fit(f_demod, alpha)
    yth = coeffs[1]*xth**coeffs[0]
    label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'
    print(label_th)
    ax.plot(xth,yth,'-',color = current_color,label = label_th)
    
# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.legend()

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$',labelpad = 5)
scalebase = 'log'
ax.set_xscale(scalebase)
ax.set_yscale(scalebase)
ax.set_title(f'{method}')
ax.set_xlim([2.7,5.3])
ax.set_ylim([1e-1,1e2])


figname = f'{fig_folder}alpha_VS_freq_exponential_method_{method}_scale_{scalebase}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')

#%% Plot dispersion relation 

set_graphs.set_matplotlib_param('single')
method = 'laser'
thickness_list = [2.5,5.0,7.5,10.0,12.5,15.0]

bounds = np.array([1.25,3.75,6.25,8.75,11.25,13.75,16.25]) # bounds used for colorbar
norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
cmap = new_blues
marker_list = ['o','s','D','^','h','p']
fig, ax = plt.subplots()


for i,h in enumerate(thickness_list):
    valid_keys = []
    for key in data.keys():
        if data[key]['method'] == method and data[key]['h'] == h:
            valid_keys.append(key)
            
    k0 = [data[key]['k0'] for key in valid_keys]
    err_k0 = [data[key]['err_k0'] for key in valid_keys]
    f_demod = [data[key]['f_demod'] for key in valid_keys]
    
    
    current_color = cmap(norm(h))
    marker = marker_list[i]
    ax.errorbar(k0,f_demod,xerr = err_k0,fmt = marker,color = current_color,markeredgecolor = 'k',ecolor = 'k',
                markersize = 8)

g = 9.81
H = 1e-1
kth = np.linspace(10,200,200)
# yth = np.sqrt(g*kth)/2/np.pi
yth = np.sqrt(g*kth*np.tanh(kth*H))/2/np.pi
ax.plot(kth,yth,'r--')

# create a scalar mappable
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$h \; \mathrm{(mm)}$')

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
scalebase = 'linear'
ax.set_xscale(scalebase)
ax.set_yscale(scalebase)
ax.set_title(f'{method}')
ax.set_xlim([25,150])
ax.set_ylim([2.5,5.5])


figname = f'{fig_folder}freq_VS_k_method_{method}_scale_{scalebase}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
plt.savefig(f'{figname}.svg', bbox_inches = 'tight')