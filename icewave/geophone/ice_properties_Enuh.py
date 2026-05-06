# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 10:56:17 2026

@author: sebas

This script enables the quick computation of sea ice mechanical properties from a line of geophones : 
    - Young's modulus E
    - Poisson coefficient nu
    - ice thickness h

First velocities of QS0 modes (longitudinal/acoustic) and SH0 modes (shear) must have been computed 
We also need to extract coordinates (f_mode,k_mode) of flexural waves dispersion relation. 

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pickle 
import scipy

import icewave.geophone.package_geophone as geopack
import icewave.sebastien.set_graphs as set_graphs

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]
#%% Import data 

year = '2025'
date = '0204' #date format, 'mmdd'
acqu_numb = '0005' #acquisition number 

path2data = os.path.join('E:/Rimouski_2025/Data/',date,'Geophones/')
signal_length = 1

fig_folder = path2data + acqu_numb + '/Results/' # folder where figures are saved 

if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

# load c_QS0 and c_SH0
file2load = f'{path2data}Phase_velocity_dictionnary_acqu_{acqu_numb}_sig_length_' + str(signal_length).replace('.','p') + '.pkl'
with open(file2load,'rb') as pfile:
    s = pickle.load(pfile)

# load (f,k) flexural_mode
for direction in ([1,2]):
    key_dir = f'dir{str(direction)}'
    file2load =  f'{path2data}{year}_{date}_acq{acqu_numb}disp_QS_dir{str(direction)}.pkl'
    with open(file2load,'rb') as pfile:
        f_mode,k_mode = pickle.load(pfile)
        s[key_dir]['f'] = f_mode
        s[key_dir]['k'] = k_mode

#%% Combine velocities from different directions 

results = {}
results['C_shear'] = 0.5*(s['dir1']['C_shear'] + s['dir2']['C_shear'])
results['C_longi'] = 0.5*(s['dir1']['C_longi'] + s['dir2']['C_longi'])
results['uC_shear'] = 0.5*(s['dir1']['uC_shear'] + s['dir2']['uC_shear'])
results['uC_longi'] = 0.5*(s['dir1']['uC_longi'] + s['dir2']['uC_longi'])


###################################################
#%% Computes Young modulus & Poisson coefficient 
###################################################

results['rho_ice'] = 917


def compute_nu(C_shear,C_longi):
    nu = 1-2*(C_shear/C_longi)**2
    return nu

def compute_E(C_longi,nu,rho_ice):
    E = rho_ice*C_longi**2*(1-nu**2)
    return E

def compute_unu(C_shear,C_longi,uC_shear,uC_longi):
    unu = 4*C_shear/(C_longi**2) * uC_shear + 4*C_shear**2/(C_longi**3) * uC_longi
    return unu

def compute_uE(C_longi,nu,rho_ice,uC_longi,unu):
    uE = 2*rho_ice * (C_longi*(1-nu**2)*uC_longi + (C_longi**2)*nu*unu)
    return uE

results['nu'] = compute_nu(results['C_shear'],results['C_longi'])
results['E'] = compute_E(results['C_longi'],results['nu'],results['rho_ice'])
results['unu'] = compute_unu(results['C_shear'],results['C_longi'],results['uC_shear'],results['uC_longi'])
results['uE'] = compute_uE(results['C_longi'],results['nu'],results['rho_ice'],
                           results['uC_longi'],results['unu'])

###################################
#%% --- Invert ice thickness -- 
###################################

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 

# direction = 1
# key_dir = f'dir{str(direction)}'
# f_mode = s[key_dir]['f']
# k_mode = s[key_dir]['k']

# interpolate (f,k) points from dispersion relation in both directions
mini = max(min(s['dir1']['f']),min(s['dir2']['f']))
maxi = min(max(s['dir1']['f']),max(s['dir2']['f']))

freq = np.linspace(mini,maxi , 30)

# Interpolate kQS1 at the new frequency points, without extrapolation
interpolator1 = scipy.interpolate.interp1d(s['dir1']['f'], s['dir1']['k'])
kQS1_interp = interpolator1(freq)

# Interpolate kQS2 at the new frequency points, without extrapolation
interpolator2 = scipy.interpolate.interp1d(s['dir2']['f'], s['dir2']['k'])
kQS2_interp = interpolator2(freq)

# Calculate the average values, ignoring NaNs
kQS = np.nanmean([kQS1_interp, kQS2_interp], axis=0)
# kQS = kQS1_interp

# Plot the results
fig, ax = plt.subplots()
ax.plot(s['dir1']['k'], s['dir1']['f'], '.',color = 'tab:blue', label='kQS1')
ax.plot(s['dir2']['k'], s['dir2']['f'],  '.',color = 'tab:orange', label='kQS2')
ax.plot(kQS1_interp, freq,  '--',color = 'tab:blue', label='kQS1 interpolated')
ax.plot(kQS2_interp, freq, '--', color = 'tab:orange',label='kQS2 interpolated')
ax.plot(kQS, freq, 'r^', label='Average kQS')
ax.set_xlabel(r'$k_{QS} \: \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \: \mathrm{(Hz)}$')
ax.legend()

#%% Find ice thickness
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h[i], 
                                                              results['E'], results['nu'],freq,c_w,rho_w)
    error = np.nansum(np.power((kQS-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
print(h_ice)

#%% Save clean graph

k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h_ice, 
                                                          results['E'], results['nu'],freq,c_w,rho_w)

set_graphs.set_matplotlib_param('single')

label_th = r'$h = ' + f'{h_ice:.2f}' '\; \mathrm{m}$'
fig, ax = plt.subplots()
ax.plot(k_mode_synthetic,freq,color = 'tab:orange',label = label_th)
ax.plot(kQS,freq,'o', color = 'tab:blue',mec = 'k')
ax.set_xlabel('$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel('$f \; \mathrm{(Hz)}$')
ax.legend()

# ax.set_xscale('log')
# ax.set_yscale('log')

figname = f'{fig_folder}hice_inversion_{date}_acq_{acqu_numb}'
# plt.savefig(figname + '.pdf', bbox_inches = 'tight')
# plt.savefig(figname + '.png', bbox_inches = 'tight')

results['h_ice'] = h_ice

#%% Save results 

file2save = f'{path2data}{year}_{date}_acq{acqu_numb}_results_inversion.pkl'
with open(file2save,'wb') as pfile:
    pickle.dump(results,pfile)

#%% Check evolution of fit with ice thickness
fig, ax = plt.subplots()
ax.plot(h,l2_norm)
ax.set_xlabel('$h \; \mathrm{(m)}$')
ax.set_ylabel('$||k_{th} - k_{exp}||_2$')



# =============================================================================
#%% Create plots superposing KF spectrum and theory
# =============================================================================

#%% Load FK dictionnary 

file_FK = f'{path2data}{year}_{date}_acq{acqu_numb}_FK_dict.pkl'
with open(file_FK,'rb') as pf:
    FK_dict = pickle.load(pf)


#%% Create subplot 3 panels

direction = 1
c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 

# compute best thickness that matches detected points on flexural wave dispersion relation
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h[i], 
                                                              results['E'], results['nu'],s[key_dir]['f'],c_w,rho_w)
    error = np.nansum(np.power((s[key_dir]['k']-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
print(h_ice)

#%%
set_graphs.set_matplotlib_param('square')
fig, axs = plt.subplots(ncols = 3, figsize = (14,6), constrained_layout = True)

for i in range(len(axs)):
    if i == 0:
        ax = axs[0]
        composante = 'N'
        channel = 1
        key = f'composante_{composante}_channel_{channel}_dir_{direction}'

        ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
                  extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
        ax.set_ylim([0,400])

        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

        # plot theory
        x_th = np.linspace(FK_dict[key]['k'].min(),FK_dict[key]['k'].max(),100)
        y_th = results['C_longi']*x_th/2/np.pi

        ax.plot(x_th,y_th,'--', color = 'tab:red',lw = 3)

    elif i == 1:
        ax = axs[1]
        
        composante = 'E'
        channel = 0
        key = f'composante_{composante}_channel_{channel}_dir_{direction}'

        ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
                  extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
        ax.set_ylim([0,400])

        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

        # plot theory
        x_th = np.linspace(FK_dict[key]['k'].min(),FK_dict[key]['k'].max(),100)
        y_th = results['C_shear']*x_th/2/np.pi

        ax.plot(x_th,y_th,'--', color = 'tab:red',lw = 3)
        
    else:
        ax = axs[i]
        composante = 'Z'
        channel = 2
        key = f'composante_{composante}_channel_{channel}_dir_{direction}'
        
        ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
                  extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
        ax.set_ylim([0,200])

        ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
        ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
        
        # # plot detected points
        # key_dir = f'dir{direction}'
        # ax.plot(s[key_dir]['k'],s[key_dir]['f'],'.', color = 'tab:red')
        
        # plot theory
        freq_th = np.linspace(s[key_dir]['f'].min(),s[key_dir]['f'].max(),100)
        k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h_ice, 
                                                                  results['E'], results['nu'],freq_th,c_w,rho_w)
        ax.plot(k_mode_synthetic,freq_th,'--',color = 'tab:red',lw = 3,label = label_th)

figname = f'E:/PhD_manuscript/ch5/geophone_BicWin/{year}_{date}_acq_{acqu_numb}_dispersion_relation_QS0_SH0_QS'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')



#%% Create plot for flexural wave

composante = 'Z'
channel = '2'
direction = 1
key = f'composante_{composante}_channel_{channel}_dir_{direction}'

c_w = 1450 # sound celerity in water 
rho_w = 1027 # density of water 

fig, ax = plt.subplots()
ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
          extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
ax.set_ylim([0,200])

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

# plot detected points
key_dir = f'dir{direction}'
ax.plot(s[key_dir]['k'],s[key_dir]['f'],'.', color = 'tab:red')

# compute best thickness value that matches current curve
h_precision = 0.01 # precision on ice thickness (in meter)
h = np.arange(0.1,1.0,h_precision) # array of height tested 

# find h_ice that minimzes distance to data points 
l2_norm = np.zeros(len(h))
for i in range(len(h)):
    k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h[i], 
                                                              results['E'], results['nu'],s[key_dir]['f'],c_w,rho_w)
    error = np.nansum(np.power((s[key_dir]['k']-k_mode_synthetic),2))
    l2_norm[i] = np.sqrt(error)
    #plt.plot(f_mode, k_mode_synthetic, color='b', label='Line between points')  
h_ice = h[np.argmin(l2_norm)] # keep thickness of ice that minimizes the error
print(h_ice)

# plot theory
freq_th = np.linspace(s[key_dir]['f'].min(),s[key_dir]['f'].max(),100)
k_mode_synthetic, k_QS0, k_SH0, cphQS = geopack.wavenumbers_stein(results['rho_ice'], h_ice, 
                                                          results['E'], results['nu'],freq_th,c_w,rho_w)
ax.plot(k_mode_synthetic,freq_th,'--',color = 'tab:orange',lw = 3,label = label_th)


#%% Create plot for SH0 

composante = 'E'
channel = 0
direction = 1
key = f'composante_{composante}_channel_{channel}_dir_{direction}'

fig, ax = plt.subplots()
ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
          extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
ax.set_ylim([0,400])

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

# plot theory
x_th = np.linspace(FK_dict[key]['k'].min(),FK_dict[key]['k'].max(),100)
y_th = results['C_shear']*x_th/2/np.pi

ax.plot(x_th,y_th,'--', color = 'tab:red',lw = 3)

#%% Create plot for QS0

composante = 'N'
channel = 1
direction = 1
key = f'composante_{composante}_channel_{channel}_dir_{direction}'

fig, ax = plt.subplots()
ax.imshow(np.transpose(FK_dict[key]['FK']),aspect = 'auto', origin = 'lower',cmap = 'gray_r',
          extent = extents(FK_dict[key]['k']) + extents(FK_dict[key]['f']),vmin = 0,vmax = 1)
ax.set_ylim([0,400])

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

# plot theory
x_th = np.linspace(FK_dict[key]['k'].min(),FK_dict[key]['k'].max(),100)
y_th = results['C_longi']*x_th/2/np.pi

ax.plot(x_th,y_th,'--', color = 'tab:red',lw = 3)







    
