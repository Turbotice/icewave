# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 11:47:10 2026

@author: sebas

This script enables to compute waves attenuation law alpha(f) from the computation of 
vertical velocity profile, obtained using one UAV uz(x,t)

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

#----------------------------------------------------------------------------------------------------------
def affine(x,a,b):
    y = a*x + b
    return y

def powerlaw_fit(x,y,err_y):
    """ Fit data using a power law, taking into account standard deviation of y """
    log_x = np.log(x)
    log_y = np.log(y)
    err_log_y = err_y/y
    
    popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y,sigma = err_log_y,absolute_sigma = True)
    err_affine = np.sqrt(np.diag(pcov))
    beta = popt[0]
    err_beta = err_affine[0]
    B = np.exp(popt[1])
    err_B = B*err_affine[1]
    
    coeffs = (beta,B)
    err_coeffs = (err_beta,err_B)
    return coeffs,err_coeffs

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

fig_folder = f'{path2data}Figures_attenuation_uz/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Perform 2D FFT of spatiotemporal signal 

facq = np.array([1/S['SCALE']['fx'],S['SCALE']['facq_t']])
shift, k, omega = FT.fft_2D(S['uz'], facq,add_pow2=[1,0])

shift = np.flip(shift,axis = 0)
# keep only positive values 
FK = abs(shift[shift.shape[0]//2:,shift.shape[1]//2:])
k_pos = k[len(k)//2:]
omega_pos = omega[len(omega)//2:]
f_pos = omega_pos/2/np.pi

Efk = {}
Efk['E'] = FK.T
Efk['k'] = k_pos
Efk['f'] = f_pos


#%% Plot FK spectrum

fig,ax,c,cbar = att_mod.plot_FK_spectrum(Efk)
c.set_clim(vmin = 2e-5,vmax = 2e-2)
cbar.set_label('$|\hat{u}_z| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}uz_spacetime_spectrum_raw_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Detect peaks from FK spectrum, for fixed k, fit over frequencies

gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.01,2.0] # range of wavevector spanned
frequency_range = [0,1.0] # range of frequency over which we look for peaks

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
figname = f'{fig_folder}uz_FK_spectrum_time_detection_hw_{hw_txt}_{suffixe}'
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

gaussian_width = 6
N = 6
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.4} # parameters for find peaks function
frequency_range = [0.1,1.0] 
wavevector_range = [0.05,2.5]
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
figname = f'{fig_folder}uz_FK_spectrum_space_detection_hw_{hw_txt}_{suffixe}'
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
main_results['attenuation_uz'] = m

file2save = f'{path2data}main_results_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(main_results,pf)

print(f'{file2save} file saved !')
    
    



#%% Detect peaks on FK spectrum 

# Build Gaussian to be convoluted (can change gaussian width to detect more or less peaks)
gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.1,2] # range of wavevector spanned
frequency_range = [0,1] # range of frequency over which we look for peaks

Speaks = att_mod.matrix_peak_correlation_detection(FK.T, k_pos, f_pos, wavevector_range, frequency_range, 
                                      gaussian, detec_param)
# change key name
Speaks['k'] = Speaks['x']
Speaks['f'] = Speaks['y']
del Speaks['x']
del Speaks['y']

#%% Filter detected points, plot results and save filtered peaks

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(Speaks['k'], Speaks['f'],Speaks)

fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')
c = ax.imshow(FK.T, cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (k_pos[0],k_pos[-1],f_pos[0],f_pos[-1]))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_time_fixed_k_{suffixe}.h5'
rw.save_dict_to_h5(filtered_properties, file2save)
print('DONE.')

#%% Compute temporal attenuation from Efk using Lorentzian fit 

E = {}
E['E'] = FK.T
E['x'] = 2*np.pi*f_pos
E['y'] = k_pos

peaks_dico = {}
peaks_dico['peaks'] = filtered_properties['peaks']
peaks_dico['x'] = 2*np.pi*filtered_properties['f']
peaks_dico['y'] = filtered_properties['k']

x_range = 2*np.pi*np.array(frequency_range)
y_range = wavevector_range

rel_height = detec_param['rel_height']
sub_folder = f'{fig_folder}time_lorentzian_fit/'
if not os.path.isdir(sub_folder):
    os.mkdir(sub_folder)
    
m = att_mod.attenuation_from_FK_spectrum(E, peaks_dico, x_range, y_range, gaussian, rel_height, sub_folder)

del E
del peaks_dico

m['f'] = m.pop('x')/2/np.pi
m['err_f'] = m.pop('err_x')/2/np.pi
m['k'] = m.pop('y')
m['lambda'] = m.pop('alpha')
m['err_lambda'] = m.pop('err_alpha')
print(m.keys())

#%% Compute water depth hw

fun = lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi
hw,err_hw = fit_water_height(m['f'], m['k'], fun, err_f = m['err_f'])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.errorbar(m['k'],m['f'],yerr = m['err_f'],fmt = '.',color = 'w')
c = ax.imshow(FK.T, cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (k_pos[0],k_pos[-1],f_pos[0],f_pos[-1]))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{u}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'

k_fit = np.linspace(kbounds[0],kbounds[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
ax.legend()

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}uz_FK_spectrum_time_detection_hw_{hw_txt}_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw


#%% Compute spatial attenuation coefficient 

time_att = np.column_stack((m['lambda'],m['err_lambda'])) # time attenuation coefficient and standard deviation
water_depth = [hw , err_hw] # water depth and standard deviation
alpha,err_alpha = att_mod.time2space_attenuation(time_att, m['k'], water_depth) # compute spatial attenuation

xbounds = [0.1,1]
ybounds = [1e-3,1e0]
x_fit = np.linspace(xbounds[0],xbounds[1],100)

# fit by a powerlaw
coeffs,err_coeffs = powerlaw_fit(m['f'], alpha, err_alpha)
yth = coeffs[1]*x_fit**coeffs[0]

label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'

fig, ax = plt.subplots()
ax.errorbar(m['f'],alpha,yerr = err_alpha,xerr = m['err_f'],fmt = 'o')
ax.plot(x_fit,yth,linewidth = 2,label = label_th)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_from_time_detection_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  

# Save data 
m['alpha'] = alpha
m['err_alpha'] = err_alpha
m['power_law'] = {}
m['power_law']['B'] = coeffs[1]
m['power_law']['err_B'] = err_coeffs[1]
m['power_law']['beta'] = coeffs[0]
m['power_law']['err_beta'] = err_coeffs[0]

file2save = f'{fig_folder}attenuation_data_time_detection_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)