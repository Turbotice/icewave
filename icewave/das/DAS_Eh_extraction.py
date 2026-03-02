# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 10:50:14 2025

@author: sebas

This script is an attempt to extract both Young modulus and ice thickness by combining both active and passive methods 

"""


import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors 
import matplotlib.cm as cm

import h5py 
import glob
import os 
import pickle

from datetime import datetime
import pytz

import scipy

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.das.DAS_package as DS
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw

# PARULA COLORMAP 
parula_map = matcmaps.parula()

# plt.rcParams.update({
#     "text.usetex": True}) # use latex

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

global g
g = 9.81

#%% Function section 

def load_DAS_water_height(file2water):
    DAS_water_height = rw.load_dict_from_h5(file2water)
    print(DAS_water_height.keys())
    # convert UTC string to datetime 
    format_date = '%Y-%m-%dT%H-%M-%S.%f'
    # convert UTC string to datetime 
    DAS_water_height['UTC_datetime'] = []
    for date_txt in DAS_water_height['UTC_t']:
        if date_txt != 'None' :
            datetime_obj = datetime.strptime(date_txt,format_date)
            datetime_obj = datetime_obj.replace(tzinfo = pytz.timezone('UTC'))
        # else :
        #     datetime_obj = None
        DAS_water_height['UTC_datetime'].append(datetime_obj)
        
    DAS_water_height['UTC_t'] = np.array(DAS_water_height['UTC_datetime'])
    del DAS_water_height['UTC_datetime']
    
    return DAS_water_height


def get_water_height(DAS_water_height,UTC_t,selected_x):
    """ Compute water height for a given UTC datetime and position along the fiber
    Input : - DAS_water_height, dictionnary containing recording of water depth evolution during the whole experiment 
            - UTC_t, datetime object, time at which we look for water height
            - selected_x, float, position along optical fiber in meter 
    Output : H, float, water height at the selected time and position 
    """

    closest_idx = np.argmin(abs(UTC_t - DAS_water_height['UTC_t']))
    closest_idx_pos = np.argmin(abs(selected_x - DAS_water_height['s']))
    H = DAS_water_height['water_height'][closest_idx_pos,closest_idx]
    return H

def shallow_hydroelastic(k,D,rho_w,H):
    """ Compute shallow hydroelastic dispersion relation
    Inputs: - k, numpy array, wavevector array
            - D, float  or numpy array, flexural modulus
            - rho_w, float or numpy array, water density
            - H, float or numpy array, water depth 
            
    Outputs : - omega, numpy array, pulsation given by shallow hydroelastic dispersion relation """
    
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega



def wavenumbers_stein(rho_ice, h, E, nu,freq,c_w,rho_w):
    """ This function computes the wave vectors associated to a given array of frequencies
    It takes as arguments : 
        - rho_ice : ice density 
        - h : a given thickness of ice 
        - E : Young modulus of ice
        - nu : Poisson coefficient of ice
        - freq : an array of frequencies, to which will correspond wave vectors 
        - c_w : waves phase velocity
        - rho_w : water density 
        
    The function returns : 
        - k_QS : wave vectors of the flexural mode
        - k_QS0 : wave vecotrs of the acoustic mode
        - k_SH0 : wave vectors of the shear mode
        - cphQS : phase velocity of the flexural mode"""

    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    # k = np.linspace(1e-3,2.0,1000)
    k = np.linspace(1e-2,2,10000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for kf in range(len(freq)):
        omeg = 2*np.pi*freq[kf]
        if omeg == 0:
            flag = 1;
            idx_flag = kf;
        else:
            cph = omeg/k # phase velocity
            # Ludovic version
            # func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            # Sebastien version (Stein 1998)
            func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_ice/D + pow(omeg/cph,4)
            
            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[kf] = (np.where(np.diff(np.signbit(func)))[0].item()) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero] # wave vector associated to flexural mode 
    if flag:
        k_QS[idx_flag] = 0
    
    return k_QS


def Stein(k,args):
    """ args = (rho_ice,h,D,omega,c_w,rho_w) """
    
    rho_ice = args[0]
    h = args[1]
    D = args[2]
    omega = args[3]
    c_w = args[4]
    rho_w = args[5]
    
    cph = omega/k
    
    func = rho_w/D*(g-omega/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omega**2*rho_ice/D + pow(omega/cph,4)
    return func
    

def wavenumbers_stein_seb(rho_ice, h, E, nu,freq,c_w,rho_w):
    
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus
    
    k_array = np.zeros(len(freq))
    for i in range(len(freq)):
        omega = 2*np.pi*freq[i]
        args = [rho_ice,h,D,omega,c_w,rho_w]
        
        x0 = 0.5
        k_QS = scipy.optimize.fsolve(Stein,x0,args = args)
        print(k_QS)
        
        k_array[i] = k_QS.item()

    return k_array 

def residuals(params,model_func,data_x,data_y):
    
    return model_func(params,data_x) - data_y

def shallow_water_detailed(k,E,h,nu,rho_w,H):
    
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus
    omega = np.sqrt((g*k + D/rho_w*k**5)*np.tanh(H*k))
    
    return omega
#%% Load dispersion relation obtained from passive method 
date = '0211'
main_path = 'U:/'
path2data = f'{main_path}Data/{date}/DAS/'

# save results dictionnary
file2load = f'{path2data}avg_disp_curv_from_CWT_{date}_with_swell_correction.h5'
CWT_results = rw.load_dict_from_h5(file2load)

#%% Load dispersion relation obtained from active method 

file2load = os.path.join(path2data, f'new_disp_curv_{date}.h5')
active_results = rw.load_dict_from_h5(file2load)

#%% Load DAS water height 

file2water = glob.glob(f'{path2data}fiber_water_height_GPS_structure_{date}.h5')[0]
print(file2water)    
DAS_water_height = load_DAS_water_height(file2water)

#%% Define fig_folder

fig_folder = f'{path2data}Eh_extraction/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Plot dispersion relation for same x - position 

xvals = np.array(list(active_results['xposition'].values()))
active_xkeys = list(active_results['xposition'].keys())

# select a given position
idx = 0
xpos = xvals[idx]
print(xpos)

cwt_freq = CWT_results[str(xpos)]['freq']
cwt_k = CWT_results[str(xpos)]['mean_k']

active_freq = active_results['freq'][active_xkeys[0]]
active_k = active_results['k_QS'][active_xkeys[0]]

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(cwt_k,cwt_freq,'o', label = 'passive')
ax.plot(active_k,active_freq,'o',label = 'active')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

ax.set_title(f'x = {xpos:.1f} m ')
ax.legend()

pos_txt = f'{xpos:.1f}'
pos_txt = pos_txt.replace('.','p')
figname = f'{fig_folder}disp_relation_superposition_active_passive_{date}_x{pos_txt}'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%%    

nu = 0.3
rho_ice = 917
rho_w = 1000
c_w = 1450 # sound celerity in water 

N = len(DAS_water_height['UTC_t'])
UTC_mid = DAS_water_height['UTC_t'][N//2]

# fit CWT results by hydroelastic dispersion relation 
H = get_water_height(DAS_water_height, UTC_mid, xpos)
print(H)

# compute best parameters using passive sources

model_func = lambda params,k : shallow_water_detailed(k,params[0],params[1],nu,rho_w,H)/2/np.pi

x0 = (3e9,0.5)
bounds = [(5e8,0.1),(1e10,2)]
test = model_func(x0,cwt_k)
result = scipy.optimize.least_squares(residuals, x0, args=(model_func, cwt_k, cwt_freq),bounds = bounds)
best_params = result.x

# show dispersion relation 
theory = {'passive':{},'active':{}}
theory['passive']['k'] = np.linspace(2e-2,2,1000)
theory['passive']['f'] = shallow_water_detailed(theory['passive']['k'],best_params[0],best_params[1],nu,rho_w,H)/2/np.pi


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(cwt_k,cwt_freq,'o', label = 'passive')
ax.plot(active_k,active_freq,'o',label = 'active')
ax.plot(theory['passive']['k'],theory['passive']['f'],'-', color = 'tab:red',label = 'hydroelastic')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

ax.set_title(f'x = {xpos:.1f} m ')
ax.legend()

figname = f'{fig_folder}disp_relation_superposition_active_passive_{date}_x{pos_txt}_with_hydroelastic'
plt.savefig(f'{figname}.png',bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')

#%% compute best parameters using active sources 

model_func_active = lambda params,freq : wavenumbers_stein_seb(rho_ice, params[1], params[0], nu, freq, c_w, rho_w)

x0 = (1e9,0.5)
bounds = [(5e8,0.1),(1e10,2)]
result_active = scipy.optimize.least_squares(residuals, x0, args=(model_func_active, active_freq, active_k),
                                             bounds = bounds)
best_params_active = result_active.x

theory['active']['f'] = np.linspace(1e-2,1e2,1000)
theory['active']['k'] = wavenumbers_stein_seb(rho_ice,best_params_active[1],best_params_active[0],nu,
                                              theory['active']['f'],c_w,rho_w)


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(cwt_k,cwt_freq,'o', label = 'passive')
ax.plot(active_k,active_freq,'o',label = 'active')
ax.plot(theory['active']['k'],theory['active']['f'],'-', color = 'tab:red',label = 'Stein')
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')

ax.set_title(f'x = {xpos:.1f} m ')
ax.legend()

