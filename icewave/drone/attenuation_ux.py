# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 16:09:51 2026

@author: sebas

This script enables to compute waves attenuation law alpha(f) from the computation of 
horizontal velocity field ux(x,y,t)

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy

import pickle
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.weather as weather
import icewave.drone.drone_projection as dp
import icewave.drone.drone_tools as drone_tools
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter
import icewave.tools.rw_data as rw
import icewave.drone.attenuation_module as att_mod

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% FUNCTION SECTION 


def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN

#----------------------------------------------------------------------------------------------------------
def fit_water_height(f,k,fun,err_f = None):
    """ Compute water height from dispersion relation fit 
    Inputs : - f, array like, frequencies
             - k, array like, wavevectors
             - fun, function of k and other parameters, hw must be the first parameter of this function 
             - err_f, optional, error on f 
    Outputs : - hw and err_hw, water height and its standard deviation computed from fit """
    popt,pcov = scipy.optimize.curve_fit(fun,k,f,sigma = err_f,absolute_sigma = True)
    
    hw = popt[0]
    err_hw = np.sqrt(np.diag(pcov)[0])  
    return hw,err_hw



#%% Import data

base = 'F:/Rimouski_2024/Data/'
date = '0226'
drone_ID = 'mesange'
exp_ID = '10-waves_005'

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
suffixe = f'{date}_{drone_ID}_{exp_ID}'
filelist = glob.glob(f'{path2data}uz_ux*scaled.h5')
print(filelist)

file2load = filelist[0]
S = rw.load_dict_from_h5(file2load)

#%% Set fig_folder 

fig_folder = f'{path2data}Figures_attenuation_ux/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Compute space-time spectrum for ux

N = S['ux'].shape[2]
Efk = FT.space_time_spectrum(S['ux'],1/S['SCALE']['fx'],S['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Plot space time spectrum 

fig,ax,c,cbar = att_mod.plot_FK_spectrum(Efk)

figname = f'{fig_folder}ux_spacetime_spectrum_raw_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Detect peaks from FK spectrum, for fixed k, fit over frequencies

gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.01,2.5] # range of wavevector spanned
frequency_range = [0,1.5] # range of frequency over which we look for peaks

file2save = f'{fig_folder}Filtered_peaks_time_fixed_k_{suffixe}.h5'

filtered_properties = att_mod.extract_peaks_fixed_k(Efk,gaussian,wavevector_range,
                                            frequency_range,file2save,detec_param)

#%% Get temporal attenuation from FK spectrum 
m = att_mod.temporal_attenuation(Efk, filtered_properties, frequency_range, wavevector_range, 
                         gaussian, detec_param, fig_folder)

# Fit water height from dispersion relation curve

fun = lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi
hw,err_hw = fit_water_height(m['f'], m['k'], fun, err_f = m['err_f'])

fig, ax, c, cbar = att_mod.plot_FK_spectrum(Efk)
title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
k_fit = np.linspace(wavevector_range[0],wavevector_range[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
ax.legend()

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}ux_FK_spectrum_time_detection_hw_{hw_txt}_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw
    
m = att_mod.structure_space_attenuation(m)

xbounds = [0.1,1]
ybounds = [1e-3,1e0]
figname = f'{fig_folder}Spatial_attenuation_from_time_detection_{suffixe}'
att_mod.plot_attenuation_power_law(m, xbounds, ybounds, figname)

# save structure
file2save = f'{fig_folder}attenuation_data_time_detection_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)

# =============================================================================
# %% Detect peaks using fixed frequency and determining directly spatial attenuation
# =============================================================================

gaussian_width = 2
N = 6
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.4} # parameters for find peaks function
frequency_range = [0.1,1.0] 
wavevector_range = [0.05,3.0]
file2save = f'{fig_folder}Filtered_peaks_space_fixed_f_{suffixe}.h5'

filtered_properties = att_mod.extract_peaks_fixed_f(Efk, gaussian, wavevector_range, frequency_range, 
                                            file2save, detec_param)

#%% Compute spatial attenuation 

m = att_mod.spatial_attenuation(Efk, filtered_properties, gaussian, wavevector_range,
                        frequency_range, detec_param, fig_folder)

# Fit water height from dispersion relation curve

fun = lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi
hw,err_hw = fit_water_height(m['f'], m['k'], fun)

fig, ax, c, cbar = att_mod.plot_FK_spectrum(Efk)
title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
k_fit = np.linspace(wavevector_range[0],wavevector_range[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
ax.legend()

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}ux_FK_spectrum_space_detection_hw_{hw_txt}_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# fit attenuation law by a power law 
coeffs,err_coeffs = att_mod.powerlaw_fit(m['f'], m['alpha'], m['err_alpha'])

# Save data 
m['power_law'] = {}
m['power_law']['B'] = coeffs[1]
m['power_law']['err_B'] = err_coeffs[1]
m['power_law']['beta'] = coeffs[0]
m['power_law']['err_beta'] = err_coeffs[0]

# plot power law
xbounds = [0.1,1]
ybounds = [1e-3,1e0]
figname = f'{fig_folder}Spatial_attenuation_from_spatial_detection_{suffixe}'
att_mod.plot_attenuation_power_law(m, xbounds, ybounds, figname)

file2save = f'{fig_folder}attenuation_data_space_detection_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)

# =============================================================================
# %% Collect attenuation structure and build a main structure
# =============================================================================

# load both set of coordinates (f,k)
m = {}
file2save = f'{fig_folder}attenuation_data_time_detection_{suffixe}.pkl'
with open(file2save,'rb') as pf:
    m['time'] = pickle.load(pf)

file2save = f'{fig_folder}attenuation_data_space_detection_{suffixe}.pkl'
with open(file2save,'rb') as pf:
    m['space'] = pickle.load(pf)

file2save = f'{path2data}main_results_{suffixe}.pkl'
if os.path.isfile(file2save):
    print(f'{file2save} already exists, loading..')
    with open(file2save,'rb') as pf :
        main_results = pickle.load(pf)

else:
    print(f'{file2save} does not exits, creation in progress..')
    main_results = {}
    main_results['date'] = date
    main_results['drone_ID'] = drone_ID
    main_results['exp_ID'] = exp_ID
    main_results['DRONE'] = S['DRONE']
    main_results['SCALE'] = S['SCALE']
    main_results['GPS'] = S['GPS']

    # get real water height from bathymetry and tides data
    path2SRT = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'
    UTC_t0 = drone_tools.get_UTC0_from_SRT(path2SRT,drone_ID,exp_ID)
    main_results['t0_UTC'] = UTC_t0
    
    # convert string to datetime object 
    GPS_D = (main_results['GPS']['latitude'],main_results['GPS']['longitude'])
    GPS_coords = dp.image_center_gps(GPS_D,main_results['DRONE']['h_drone'],main_results['DRONE']['alpha_0'])
    real_hw = weather.get_water_height(GPS_coords,UTC_t0,disk = 'Elements',year = '2024')
    main_results['real_hw'] = real_hw

# main_results['Efk'] = Efk
main_results['attenuation_ux'] = m

file2save = f'{path2data}main_results_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(main_results,pf)

print(f'{file2save} file saved !')



# =============================================================================
# %% Compare attenuation results using both detection methods
# =============================================================================

#%% Compare dispersion relations 

m = main_results['attenuation_ux']
fig, ax = plt.subplots()
for key in m.keys():
    ax.plot(m[key]['k'],m[key]['f'],'.',label = key)
    
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.legend()

figname = f'{fig_folder}comparison_detection_methods_dispersion_relation_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compare attenuation laws
xbounds = np.array([1e-1,1e0])
ybounds = np.array([1e-3,1e0])
xfit = np.linspace(xbounds[0],xbounds[-1],100)

color_tab = {'time':'tab:blue','space':'tab:orange'}

fig, ax = plt.subplots()
for key in m.keys():
    ax.errorbar(m[key]['f'],m[key]['alpha'],yerr = m[key]['err_alpha'],fmt = '.',
                color = color_tab[key],label = key)
    # yth = m[key]['power_law']['B']*xfit**m[key]['power_law']['beta']
    # ax.plot(xfit,yth,'-',color = color_tab[key])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}comparison_detection_methods_attenuation_law_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')














