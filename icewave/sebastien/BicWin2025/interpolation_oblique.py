# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 17:44:09 2026

@author: sebas
"""

#%% 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean

import scipy
import cv2 as cv

import h5py
import pickle
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.datafolders as df
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.drone.drone_projection as dp
import icewave.drone.drone_tools as drone_tools
import icewave.tools.weather as weather
import icewave.tools.rw_data as rw
import icewave.drone.attenuation_module as att_mod

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# PARULA COLORMAP 
parula_map = matcmaps.parula()

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



#%% Load data

base = 'U:/Data/'
date = '0211'
drone_ID = 'mesange'
exp_ID = '05-waves_001'
suffixe = f'{date}_{drone_ID}_{exp_ID}'

fig_folder = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load matfile 

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
filelist = glob.glob(f'{path2data}*scaled.mat')
print(filelist)

idx_file = 0
file2load = filelist[idx_file]

# load file 
with h5py.File(file2load, 'r') as fmat:
    S = {}

    print('Top-level keys : ', list(fmat.keys()))

    S = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    S = mat2py.transpose_PIVmat_fields(S)

#%% 

print(S['DRONE']['alpha_0']*180/np.pi)

#%% Supress quadratic noise

Vx = FT.supress_quadratic_noise(np.transpose(S['Vx'],(1,0,2)),S['x'],S['y'])
Vy = FT.supress_quadratic_noise(np.transpose(S['Vy'],(1,0,2)),S['x'],S['y'])
Vx = np.transpose(Vx,(1,0,2))
Vy = np.transpose(Vy,(1,0,2))

#%% Show apparent velocity fields

extents_pix = np.array([S['PIXEL']['x_pix'].min(),S['PIXEL']['x_pix'].max(),
                    S['PIXEL']['y_pix'].min(),S['PIXEL']['y_pix'].max()])

frame = 3600

fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,8))
imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.imshow(Vx[:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_x$')
    else:
        imsh.append(ax.imshow(Vy[:,:,frame].T,cmap = parula_map))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_y$')
   
    ax.set_xlabel(r'$x_{pix}$')     

axs[0].set_ylabel(r'$y_{pix}$')

plt.tight_layout()

#%% Show georectified apparent velocity fields

cmap = cmocean.cm.balance
frame = 3000

X_bounds = np.array([-50,50])
Y_bounds = np.array([-40,40])
Vx_values = np.array([-1,1])
Vy_values = np.array([-1,2])

fig, axs = plt.subplots(ncols = 2,sharey = True,figsize = (12,8))
imsh = []
for i,ax in enumerate(axs):
    if i == 0:
        imsh.append(ax.pcolormesh(S['X'],S['Y'],Vx[:,:,frame],shading = 'gouraud',cmap = cmap,
                                 vmin = Vx_values[0],vmax = Vx_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_X$')
    else:
        imsh.append(ax.pcolormesh(S['X'],S['Y'],Vy[:,:,frame],shading = 'gouraud',cmap = cmap,
                                 vmin = Vy_values[0],vmax = Vy_values[1]))
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh[i],cax = cax)
        cbar.set_label(r'$V_Y$')

    ax.set_aspect(1)

#%% Check values of coefficients

alpha = S['DRONE']['alpha_0']
h = S['DRONE']['h_drone']
f = S['DRONE']['focale']
Y = S['Y']
X = S['X']

z_star = dp.get_zstar(h,alpha,Y)

coeff = (np.cos(alpha) + (np.sin(alpha)**2)*Y/z_star)*f/z_star

fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],coeff,shading = 'gouraud',cmap = parula_map)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
ax.set_aspect(1)

#%% Try inversion of uz

frame = 1
cmap = cmocean.cm.balance
uz = Vy/coeff[:,:,None]

uz_values = np.array([-0.08,0.08])
fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],uz[:,:,frame],shading = 'gouraud',cmap = cmap,
                    vmin = uz_values[0],vmax = uz_values[1])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_z$')
ax.set_aspect(1)

#%% Deduce ux
cmap = cmocean.cm.balance
ux = Vx*z_star[:,:,None]/f - uz*np.sin(alpha)*X[:,:,None]/z_star[:,:,None] 

fig, ax = plt.subplots()
imsh = ax.pcolormesh(S['X'],S['Y'],ux[:,:,frame],shading = 'gouraud',cmap = cmap)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_x$')
ax.set_aspect(1)

#%% Define a regular grid

minx = -100
maxx = 100
miny = -55
maxy = S['Y'].max()

artificial_facq_x = 1/0.9 # artificial spatial frequency in box/meter
grid_x, grid_y = np.meshgrid(
    np.linspace(minx, maxx, int((maxx - minx)*artificial_facq_x)),  
    np.linspace(miny, maxy, int((maxy - miny)*artificial_facq_x))
)


#%% Compare spatial sampling 
diff_y = np.diff(S['Y'][0,:])
diff_x = np.mean(np.diff(S['X'],axis = 1),axis = 1)

fig, ax = plt.subplots()
ax.plot(abs(diff_y),'.',label = '$\Delta_y$')
ax.plot(abs(diff_x),'.', label = '$\Delta_x$')
ax.set_ylim([0,2])
ax.legend()

#%% Interpolate uz over the regular grid

points = np.array([S['X'].ravel(),S['Y'].ravel()]).T
interp_uz = np.zeros((grid_x.shape[0],grid_x.shape[1],uz.shape[2]))
for frame in range(uz.shape[2]):
    print(frame)
    interp_uz[:,:,frame] = dp.interpolate_field(points,uz[:,:,frame],grid_x,grid_y) 

#%% Interpolate ux over the regular grid

points = np.array([S['X'].ravel(),S['Y'].ravel()]).T
interp_ux = np.zeros((grid_x.shape[0],grid_x.shape[1],ux.shape[2]))
for frame in range(ux.shape[2]):
    print(frame)
    interp_ux[:,:,frame] = dp.interpolate_field(points,ux[:,:,frame],grid_x,grid_y)


#%% Show interpolated uz field 

frame = 3000
fig, ax = plt.subplots()
imsh = ax.pcolormesh(grid_x.T,grid_y.T,interp_uz[:,:,frame].T,shading = 'gouraud',cmap = parula_map,
                    vmin = -0.3, vmax = 0.3)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$V_y$')
ax.set_aspect(1)


#%% Save interpolated field 

data = {}
data['DRONE'] = S['DRONE']
data['GPS'] = S['GPS']
data['uz'] = uz
data['interp_uz'] = interp_uz
data['interp_ux'] = interp_ux
data['grid_x'] = grid_x
data['grid_y'] = grid_y
data['SCALE'] = S['SCALE']
data['artificial_facq_x'] = artificial_facq_x
data['t'] = S['t']

file2save = f'{path2data}interpolated_u_quadratic_correction{date}_{drone_ID}_{exp_ID}.h5'
rw.save_dict_to_h5(data, file2save)



#%% Load interpolated field 

path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = f'{path2data}interpolated_u_quadratic_correction{date}_{drone_ID}_{exp_ID}.h5'
S = rw.load_dict_from_h5(file2load)

#%% Set fig_folder 

fig_folder = f'{path2data}Figures_attenuation_uz/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Compute space-time spectrum for uz

N = S['interp_uz'].shape[2]
Efk = FT.space_time_spectrum(S['interp_uz'],1/S['artificial_facq_x'],
                             S['SCALE']['facq_t'],add_pow2 = [0,0,0])

#%% Plot FK spectrum

fig,ax,c,cbar = att_mod.plot_FK_spectrum(Efk)
c.set_clim(vmin = 1e-6,vmax = 1e-3)
cbar.set_label('$|\hat{u}_z| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}uz_spacetime_spectrum_raw_{suffixe}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Detect peaks from FK spectrum, for fixed k, fit over frequencies

gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.01,2.5] # range of wavevector spanned
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

gaussian_width = 3
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
    real_hw = weather.get_water_height(GPS_coords,UTC_t0,disk = 'Backup25',year = '2025')
    main_results['real_hw'] = real_hw

# main_results['Efk'] = Efk
main_results['attenuation_uz'] = m

file2save = f'{path2data}main_results_{suffixe}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(main_results,pf)

print(f'{file2save} file saved !')



# =============================================================================
# %% Compare attenuation results using both detection methods
# =============================================================================

#%% Compare dispersion relations 

m = main_results['attenuation_uz']
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
