# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 09:40:03 2026

@author: sebas
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy

import os
import glob
import pickle
import sys
sys.path.append('C:/Users/sebas/git/')

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df
import icewave.sebastien.set_graphs as set_graphs
import icewave.sebastien.theory.module_bilayer_viscous_dimensionless as theory
import icewave.drone.attenuation_module as att_mod

global g
g = 9.81

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Function section 

def get_date_drone_exp(file2load):
    
    chain = file2load.split('\\')
    date = chain[1]
    drone = chain[3]
    exp = chain[5]
    
    return date, drone, exp

def get_subset_dict(data, selection):
    # On ne garde que les expériences qui valident TOUS les critères de sélection
    return {
        key: exp_data 
        for key, exp_data in data.items()
        if all(exp_data.get(param) == value for param, value in selection.items())
    }

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

#%% Define fig_folder

fig_folder = 'F:/PhD_Manuscript/ch4/'
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)


#%% load data

path = 'F:/PhD_Manuscript/ch3/Attenuation/'
file2load = f'{path}attenuation_field_main_data.pkl'

with open(file2load,'rb') as pf:
    data = pickle.load(pf)

#%% Load lab experiments

base =  'F:/Stagiaires/Mariya/'
path2data = f'{base}attenuation_results.pkl'

with open(path2data,'rb') as pf:
    data_lab = pickle.load(pf)
    
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

ax.set_xlim([1e-1,1.2e0])
ax.set_ylim([2e-3,5e-1])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)

figname = f'{fig_folder}attenuation_all_field_observation'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
# %% Plot keeping only trustful points, also plot error
# =============================================================================

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

set_graphs.set_matplotlib_param('single')
fig, axs = plt.subplots(nrows = 1,ncols = 2,layout = 'constrained',figsize = (14,6))

for key,m in data[key_process].items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    mask = m['d'] < 0.14
    
    for i,ax in enumerate(axs):
        if i == 0:
            
            ax.errorbar(x[mask],y[mask],yerr = yerr[mask],xerr = xerr[mask],fmt = '.',label = key)
                
            ax.legend(fontsize = 10)
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            ax.set_xlim([1e-1,1.2e0])
            ax.set_ylim([2e-3,5e-1])
            
            ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
            ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
            ax.grid(True, linestyle='--', alpha=0.3)
        
        else:
            ax.plot(x,m['d'],'.')
            
# =============================================================================
# %% Plot attenuation data with only trustful points
# =============================================================================

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

ms = 10
set_graphs.set_matplotlib_param('powerpoint')
fig, ax = plt.subplots()


for key,m in data[key_process].items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    mask = m['d'] < 0.14
    
    if key == '2024_0226_mesange_10-waves_005':
        mask_freq = x < 0.74
        mask = np.logical_and(mask,mask_freq)
        
    ax.errorbar(x[mask],y[mask],yerr = yerr[mask],xerr = xerr[mask],fmt = '.',label = key,ms = ms)
        
# ax.legend(fontsize = 10)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1.5e0])
ax.set_ylim([2e-3,1e0])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)
        
figname = f'{fig_folder}attenuation_all_field_observation'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
# %% Fit attenuation data with power law 
# =============================================================================

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

ms = 10
xth = np.linspace(1e-2,1e1,200)
set_graphs.set_matplotlib_param('powerpoint')
fig, ax = plt.subplots()

x2fit = []
y2fit = []
for key,m in data[key_process].items():
    mask = m['d'] < 0.14
    
    if key == '2024_0226_mesange_10-waves_005':
        mask_freq = m['f'] < 0.74
        mask = np.logical_and(mask,mask_freq)
        
    x = m['f'][mask]
    xerr = m['err_f'][mask]
    
    y = m['alpha'][mask]
    yerr = m['err_alpha'][mask]
    
    for x_elem in x:
        x2fit.append(x_elem)
    for y_elem in y:
        y2fit.append(y_elem)
        
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',ms = ms)
        
x2fit = np.array(x2fit)
y2fit = np.array(y2fit)

# label_th = r'$\alpha = ' + f'{coeffs[1]:.1f}' + r'f^{' + f'{coeffs[0]:.1f}' + r'}$'
# ax.plot(xth,y_th,'r-.',label = label_th)

# fit by a power law 2 
beta = 2

popt,pcov = scipy.optimize.curve_fit(lambda x,b : affine(x, beta, b),np.log(x2fit*2*np.pi),np.log(y2fit),
                                     bounds = (np.log(1e-5),np.log(1e2)))
print(popt)
coeff = popt[0]
err_coeff = np.sqrt(np.diag(pcov))[0]

yth = np.exp(affine(np.log(xth*2*np.pi),beta,coeff))

B = np.exp(coeff)*(2*np.pi)**beta
label_th = r'$\alpha = ' + f'{B:.2f}' + r'f^2$'
ax.plot(xth,yth,'k-',label = label_th)


# build a cone with which data are contained
p = 3
B = 2.9*1e-3*(0.1)**(-p)
y_high = B*xth**p

p = 1.5
B = 2.3*1e-3*(0.1)**(-p)
y_low = B*xth**p

ax.fill_between(xth,y_high,y_low,color = 'grey',alpha = 0.6)


ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1.5e0])
ax.set_ylim([2e-3,1e0])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)
ax.legend()
    
# figname = f'{fig_folder}attenuation_all_field_observation_with_cone'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
#%% Rescale attenuations data  
# =============================================================================

rho_ice = 917 # in kg/m3
rho_w = 1035 # in kg/m3
H = 0.10 # in meter

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

ms = 10
set_graphs.set_matplotlib_param('powerpoint')
fig, ax = plt.subplots()

x2fit = []
y2fit = []
for key,m in data[key_process].items():
    mask = m['d'] < 0.14
    
    if key == '2024_0226_mesange_10-waves_005':
        mask_freq = m['f'] < 0.74
        mask = np.logical_and(mask,mask_freq)
        
    x = m['f'][mask]*2*np.pi
    xerr = m['err_f'][mask]*2*np.pi
    
    y = m['alpha'][mask] / H**(0.75) / ((rho_w - rho_ice)/rho_ice)**0.5 / g**(-1.5)
    yerr = m['err_alpha'][mask] / H**(0.75) / ((rho_w - rho_ice)/rho_ice)**0.5 / g**(-1.5)
    
    for x_elem in x:
        x2fit.append(x_elem)
    for y_elem in y:
        y2fit.append(y_elem)
        
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',ms = ms)
        
x2fit = np.array(x2fit)
y2fit = np.array(y2fit)

# label_th = r'$\alpha = ' + f'{coeffs[1]:.1f}' + r'f^{' + f'{coeffs[0]:.1f}' + r'}$'
# ax.plot(xth,y_th,'r-.',label = label_th)

# fit by a power law 
beta = 3

popt,pcov = scipy.optimize.curve_fit(lambda x,b : affine(x, beta, b),np.log(x2fit),np.log(y2fit),
                                     bounds = (np.log(1e-5),np.log(1e2)))
print(popt)
coeff = np.exp(popt[0])
err_coeff = np.exp(np.sqrt(np.diag(pcov))[0])

xth = np.linspace(1e-2,1e1,200)*2*np.pi
yth = coeff*xth**beta

label_th = r'$\alpha = ' + f'{coeff:.2f}' + r'\omega^{' + f'{beta:.1f}' + r'}$'
ax.plot(xth,yth,'k-',label = label_th)


ax.set_xscale('log')
ax.set_yscale('log')
xlim = np.array([6e-1,1e1])
ax.set_xlim(xlim)
ax.set_ylim([1e-1,1e3])

ax.set_xlabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')
ax.set_ylabel(r'$\alpha/ \left[ \left(\frac{\Delta \rho}{\rho}\right)^{1/2} h^{3/4} g^{-3/2} \right] \; \mathrm{(m^{-1/4}.s^{-3})}$')
ax.grid(True, linestyle='--', alpha=0.3)
ax.legend()

# compute wave amplitude 
A = 0.25*(coeff)**(-4)
print(A)
    
# figname = f'{fig_folder}attenuation_all_field_observation_rescaling'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

# =============================================================================
# %% Compare attenuation with models 
# =============================================================================

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

xth = np.linspace(1e-2,1e1,200)

set_graphs.set_matplotlib_param('double')
fig, ax = plt.subplots()

x2fit = []
y2fit = []

for key,m in data[key_process].items():
    
    mask = m['d'] < 0.14
    x = m['f'][mask]
    xerr = m['err_f'][mask]
    
    y = m['alpha'][mask]
    yerr = m['err_alpha'][mask]
    
    for x_elem in x:
        x2fit.append(x_elem)
    for y_elem in y:
        y2fit.append(y_elem)
    
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.')
    
x2fit = np.array(x2fit)
y2fit = np.array(y2fit)

# fit by a power law
coeffs,err_coeffs = powerlaw_fit(x2fit, y2fit)
y_th = coeffs[1]*xth**coeffs[0]

# label_th = r'$\alpha = ' + f'{coeffs[1]:.1f}' + r'f^{' + f'{coeffs[0]:.1f}' + r'}$'
# ax.plot(xth,y_th,'r-.',label = label_th)

# fit by a power law 2 
beta = 2

popt,pcov = scipy.optimize.curve_fit(lambda x,b : affine(x, beta, b),np.log(x2fit*2*np.pi),np.log(y2fit),
                                     bounds = (np.log(1e-5),np.log(1e2)))
print(popt)
coeff = popt[0]
err_coeff = np.sqrt(np.diag(pcov))[0]
yth = np.exp(affine(np.log(xth*2*np.pi),beta,coeff))

B = np.exp(coeff)
err_B = np.sqrt((B*err_coeff)**2)
label_th = r'$\alpha = ' + f'{B:.2f}' + r'f^2$'
# ax.plot(xth,yth,'k-',label = label_th)

# ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1.5e0])
ax.set_ylim([2e-3,1e0])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)

# compute Gpp/h
rho_w = 1e3

Gpp_h = 0.5*B*rho_w*g**2
print(Gpp_h)

Gpp = 30
err_G = 0
h = Gpp/Gpp_h

err_h = h*np.sqrt((err_G/Gpp)**2 + (err_B/B)**2)
print(f'h ={h:.3f} ± {err_h:.3f}')

h_min = 0.05
h_max = 0.20

B_min = 2*Gpp/h_min/rho_w/(g**2)
B_max = 2*Gpp/h_max/rho_w/(g**2)
y_min = B_min*(xth*2*np.pi)**2
y_max = B_max*(xth*2*np.pi)**2

# ax.fill_between(xth,y_min,y_max,color = 'k',alpha = 0.2)

# figname = f'{fig_folder}attenuation_all_field_observation'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Plot with comparison with bilayer model

# physical properties
f_th = np.linspace(1e-1,1e1,500)
rho_1 = 1e3 # water density
r = 1
rho_2 = rho_1/r
nu_2 = 1e-6 # water viscosity
h_meter = 0.3 # in meter 

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

set_graphs.set_matplotlib_param('double')
fig, ax = plt.subplots()

x2fit = []
y2fit = []

for key,m in data[key_process].items():
    x = m['f']
    xerr = m['err_f']
    
    y = m['alpha']
    yerr = m['err_alpha']
    
    for x_elem in x:
        x2fit.append(x_elem)
    for y_elem in y:
        y2fit.append(y_elem)
    
    ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)
    
x2fit = np.array(x2fit)
y2fit = np.array(y2fit)

# fit h,nu_1,nu_2
# popt,pcov = scipy.optimize.curve_fit(lambda freq,h_meter,nu_1,nu_2 : 
#                         np.imag(theory.get_wavevector_theory(freq, h_meter, rho_1, 
#                                                                 rho_2, nu_1, nu_2,display = False)),
#                         x2fit,y2fit,p0 = (1e-1,1e-3,1e-4), bounds = ([1e-2,1e-6,1e-6],[1,1e0,1e-2]))
    
# print(popt)
# print(np.sqrt(np.diag(pcov)))
# nu_2 = popt[2]
# nu_1 = popt[1]
# h_meter = popt[0]
    
# fit nu_1,nu_2
popt,pcov = scipy.optimize.curve_fit(lambda freq,nu_1,nu_2 : 
                        np.imag(theory.get_wavevector_theory(freq, h_meter, rho_1, 
                                                                rho_2, nu_1, nu_2,display = False)),
                        x2fit,y2fit,p0 = (1e-3,1e-4), bounds = ([1e-6,1e-6],[1e0,1e-2]))
    
print(popt)
print(np.sqrt(np.diag(pcov)))
nu_2 = popt[1]
nu_1 = popt[0]
    
label_th = r'Bilayer'

# plot theory 
k_th = theory.get_wavevector_theory(f_th, h_meter, rho_1, rho_2, nu_1, nu_2,display = False)
ax.plot(f_th,np.imag(k_th),'k-.',label = label_th)

# fit by a power law
coeffs,err_coeffs = powerlaw_fit(x2fit, y2fit)
y_th = coeffs[1]*f_th**coeffs[0]

label_th = r'$\alpha = ' + f'{coeffs[1]:.1f}' + r'f^{' + f'{coeffs[0]:.1f}' + r'}$'
ax.plot(f_th,y_th,'r-.',label = label_th)

# fit by a second power law, only in interval f = [0.47,0.75]
# fmin = 0.47
# fmax = 0.75
# mask = np.logical_and(x2fit > fmin, x2fit < fmax)

# # ax.plot(x2fit[mask],y2fit[mask],'.')
# # fit by a power law
# coeffs,err_coeffs = powerlaw_fit(x2fit[mask], y2fit[mask])
B = 1.2
beta = 5
y_th = B*f_th**beta

label_th = r'$\alpha = ' + f'{B:.1f}' + r'f^{' + f'{beta:.1f}' + r'}$'
ax.plot(f_th,y_th,'-.',color = 'darkorchid',label = label_th)


ax.legend(fontsize = 10)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-1,1.5e0])
ax.set_ylim([2e-3,1e0])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)

figname = f'{fig_folder}attenuation_all_field_observation_with_theory'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Plot lab experiments + field observations

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

ms = 4

set_graphs.set_matplotlib_param('double')
fig, ax = plt.subplots()

x2fit = []
y2fit = []

for key,m in data[key_process].items():
    mask = m['d'] < 0.14
    x = m['f'][mask]
    xerr = m['err_f'][mask]
    
    y = m['alpha'][mask]
    yerr = m['err_alpha'][mask]
    
    for x_elem in x:
        x2fit.append(x_elem)
    for y_elem in y:
        y2fit.append(y_elem)
    
    # ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)
    
x2fit = np.array(x2fit)
y2fit = np.array(y2fit)

ax.plot(x2fit,y2fit,'.',color = 'tab:blue',ms = ms,alpha = 1)

# plot lab experiments 

# current_dict = data_lab

# alpha = np.array([current_dict[key]['alpha'] for key in current_dict.keys()])
# err_alpha = np.array([current_dict[key]['err_alpha'] for key in current_dict.keys()])
# f = np.array([current_dict[key]['f_demod'] for key in current_dict.keys()])

rho = 1.03
selection = {'rho':rho}
current_dict = get_subset_dict(data_lab, selection)

alpha = []
err_alpha = []
f = []
for key in current_dict.keys():
    if current_dict[key]['h'] > 5.0:

        alpha.append(current_dict[key]['alpha'])
        err_alpha.append(current_dict[key]['err_alpha'])
        f.append(current_dict[key]['f_demod'])

f = np.array(f)
alpha = np.array(alpha)

ax.plot(f,alpha,'.',color = 'tab:orange',ms = ms,alpha = 1)

# perform fit 
x_general = np.concatenate((x2fit,f))
y_general = np.concatenate((y2fit,alpha))

# free power law
# coeffs,err_coeffs = powerlaw_fit(x_general,y_general)
# xth = np.linspace(1e-1,1e1,200)
# yth = coeffs[1]*xth**coeffs[0]

# imposed power law
beta = 3.0

popt,pcov = scipy.optimize.curve_fit(lambda x,b : affine(x, beta, b),np.log(x_general*2*np.pi),
                                     np.log(y_general),
                                     bounds = (np.log(1e-5),np.log(1e-1)))
print(popt)
coeff = popt[0]
err_coeff = np.sqrt(np.diag(pcov))[0]
xth = np.linspace(1e-1,1e1,200)
yth = np.exp(affine(np.log(xth*2*np.pi),beta,coeff))

label = r'$\alpha = ' + f'{np.exp(coeff)*(2*np.pi)**3:.2f}' + 'f^3$'
ax.plot(xth,yth,'k-',label = label)
ax.legend()

# fill between
# y_low = np.exp(affine(np.log(xth*2*np.pi),beta,coeff - err_coeff))
# y_high = np.exp(affine(np.log(xth*2*np.pi),beta,coeff + err_coeff))
# ax.fill_between(xth,y_high,y_low,color = 'k',alpha = 0.6)


ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True, linestyle='--', alpha=0.3)

ax.set_xlim([1e-1,1e1])
ax.set_ylim([1e-3,3e2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
# ax.legend()

figname = f'{fig_folder}comparison_attenuation_lab_field_power_law_3'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# compute Gpp
Gpp = 2*pow(np.exp(coeff),2)*pow(g,4)*rho*1e3
err_Gpp = 4*pow(g,4)*np.exp(coeff)*np.exp(coeff)*err_coeff*rho*1e3

print(f'Gpp = {Gpp:.1f} ± {err_Gpp:.1f} Pa')

#%% Plot dispersion relation for all field observations 

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for key,m in data[key_process].items():
    x = m['k']
    y = m['f']
    
    ax.plot(x,y,'.',label = key)

ax.legend(fontsize = 10)
ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xlim([1e-1,1.2e0])
# ax.set_ylim([2e-3,5e-1])

ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.grid(True, linestyle='--', alpha=0.3)

figname = f'{fig_folder}disp_relation_all_field_observation'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')



#%% Plot hw VS real_hw

component = 'ux'
dim = 'time'
key_process = f'{component}_{dim}'

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

for key,m in data[key_process].items():
    print(m.keys())
    # x = m['real_hw']
    # xerr = m['err_f']
    
    # y = m['alpha']
    # yerr = m['err_alpha']
    
    # ax.errorbar(x,y,yerr = yerr,xerr = xerr,fmt = '.',label = key)

