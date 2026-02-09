#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:47:56 2024

@author: moreaul
"""
#%% import libraries
import numpy as np
from random import random
import matplotlib.pyplot as plt
import pickle
import os , sys
import time

#%% code

year = '2024'
date = '0210' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 
direction = 1 # 1 ou 2 
equation = 'squire'

#path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/' +year+'_BICWIN/'
path2data = f'D:/Startup_kit_Stage_MSIM/data/{year}_BICWIN/' # arborescence sur ssd vasco

filename = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
filename_synthetics = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'synth.pkl'
file2load1 = path2data  + date + '/Geophones/' + filename
file2load2 = path2data  + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ '_cQS0_dir' +str(direction) +'.pkl'
file2load3 = path2data  + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_dir' +str(direction) +'.pkl'
file2save = path2data  + date + '/Geophones/' + year + '_' + date + '_acq'+acqu_numb+ '_SA_inversion' '.pkl'
fig_save_path = path2data+'/'+date+'/Geophones/'+acqu_numb



with open(file2load1, "rb") as f:
    data = pickle.load(f)
freq = data[0] 
kQS = data[1] 


with open(file2load2, "rb") as f:
    data = pickle.load(f)
cQS0 = data 
with open(file2load3, "rb") as f:
    data = pickle.load(f)
cSH0 = data 


data = {
    'freq':freq,
    'kQS': kQS,
    'cQS0': cQS0,
    'cSH0': cSH0
}



def tic():
    global start_time
    start_time = time.perf_counter()

def toc():
    elapsed_time = time.perf_counter() - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")





tic()


def candidate(X,delta_X,min_X,max_X): 
    cand = np.zeros(X.shape[0])   
    
    for k in range(X.shape[0]):
        # force jump to remain within domain bounds
        if X[k] + delta_X[k] >= max_X[k]:
            bound_max = max_X[k]            
        else:
            bound_max = X[k] + delta_X[k]
                           
        if X[(k)] - delta_X[(k)] <= min_X[k]:
            bound_min = min_X[k]
        else:
            bound_min = X[k] - delta_X[k]
        # random candidate generation in multispace delat_X box around position  
        cand[k] = bound_min + (bound_max-bound_min)*random()

    return cand    

def wavenumbers_stein_squire( rho_ice, h, H, E, nu,freq,c_w,rho_w,equation):
    
    g = 9.81
    # G = E/(2*(1+nu))
    # cS0 = np.sqrt(E/(rho_ice*(1-nu**2)))
    # cSH0 = np.sqrt(G/rho_ice)
    D = E*pow(h,3)/(12*(1-nu**2))

    k = np.linspace(1e-12,5,100000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for kf in range(len(freq)):
        omeg = 2*np.pi*freq[kf]
        if omeg == 0:
            flag = 1;
            idx_flag = kf;
        else:
            cph = omeg/k
            if equation == 'stein':
                func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4) 
            elif equation == 'squire':
                coth = 1 / np.tanh(k* H)
                func = pow(omeg, 2)*(k*h*rho_ice/rho_w + coth) - D*pow(k, 5)/rho_w - g * k 
            else:
                print('unapropriate equation name: chose between stein or squire ')
            
            func[func.imag != 0] = -1
            func = func.real
            idx_zero[kf] = (np.where(np.diff(np.signbit(func)))[0])
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero]       
    if flag:
        k_QS[idx_flag] = 0
        
    # k_QS0 = freq/cS0*2*np.pi
    # k_SH0 = freq/cSH0*2*np.pi  
    # cphQS = freq/k_QS*2*np.pi

    
    return k_QS





def simulated_annealing(delta_param0, T0, T0param, Tmin, Tminparam, X, MIN_param0, MAX_param0, data, isave):
    H = 10
    c_w = 1450
    rho_w = 1027
    freq = data['freq']
    kQSdata = data['kQS']
    cQS0data = data['cQS0']
    cSH0data = data['cSH0']
    N_SA = X.shape[1]
    T = np.exp(np.linspace(np.log(T0), np.log(Tmin), N_SA))
    variance = T**2
    Tparam = np.exp(np.linspace(np.log(T0param), np.log(Tminparam), N_SA))
    T_critic = Tmin
    it_critic = N_SA
    misfit_accepted = np.zeros(N_SA)
    likelihood = np.zeros(N_SA)
    G = X[1, 0] / (2 * (1 + X[2, 0]))
    cQS0synthetics = np.sqrt(X[1, 0] / (X[3, 0] * (1 - X[2, 0]**2)))
    cSH0synthetics = np.sqrt(G / X[3, 0])
    kQSsynthetics = wavenumbers_stein_squire(X[3, 0], X[0, 0], H, X[1, 0], X[2, 0], freq, c_w, rho_w, equation)
    misfit_accepted[0] = 1/3 * (np.mean(np.abs(kQSdata - kQSsynthetics) / kQSdata) + np.mean(np.abs(cQS0data - cQS0synthetics) / cQS0data) + np.mean(np.abs(cSH0data - cSH0synthetics) / cSH0data))
    likelihood[0] = np.exp(-misfit_accepted[0]**2 / (2 * variance[0]))
    it = 1
    it_tot = 0
    cpt_out = 0
    while it <= N_SA - 1: 
        it_tot += 1
        delta_param = delta_param0 * Tparam[it]
        success = False
        while not success:
            try:
                Xtest = candidate(X[:, it - 1], delta_param, MIN_param0, MAX_param0)
                G = Xtest[1] / (2 * (1 + Xtest[2]))
                cQS0synthetics = np.sqrt(Xtest[1] / (Xtest[3] * (1 - Xtest[2]**2)))
                cSH0synthetics = np.sqrt(G / Xtest[3])
                kQSsynthetics = wavenumbers_stein_squire(Xtest[3], Xtest[0], H, Xtest[1], Xtest[2], freq, c_w, rho_w, equation)
                success = True
            except Exception as e:
                print(f"Error: {e}. Retrying...")
        misfit_cand = 1/3 * (np.mean(np.abs(kQSdata - kQSsynthetics) / kQSdata) + np.mean(np.abs(cQS0data - cQS0synthetics) / cQS0data) + np.mean(np.abs(cSH0data - cSH0synthetics) / cSH0data))
        likelihood_cand = np.exp(-misfit_cand**2 / (2 * variance[it]))
        test_likelihood = min(1, likelihood_cand / likelihood[it - 1])
        if random() < test_likelihood:
            X[:, it] = Xtest
            misfit_accepted[it] = misfit_cand
            likelihood[it] = likelihood_cand
            cpt_out = 0
            it += 1
            print(it)
        else:
            cpt_out += 1
            if cpt_out > 200:
                I = np.argmin(misfit_accepted[0:it])
                T_critic = T[it]
                it_critic = it
                with open(file2save, 'wb') as f:
                    pickle.dump([T_critic, X, misfit_accepted], f)    
                break
        if it % isave == 0 and cpt_out < 1:
            with open(file2save, 'wb') as f:
                pickle.dump([T_critic, X, misfit_accepted], f)
            fig0, ax0 = plt.subplots()
            ax0.plot(misfit_accepted[:it], linestyle='-', color='k')
            ax0.set_ylabel('misfit vs iteration', fontsize=10)
            fig1, ax1 = plt.subplots(2, 2)
            ax1[0, 0].plot(X[0, :it], linestyle='-', color='k')
            ax1[0, 0].set_ylabel('thickness (m)', fontsize=10)
            ax1[0, 1].plot(X[1, :it], linestyle='-', color='k')
            ax1[0, 1].set_ylabel('E (GPa)', fontsize=10)
            ax1[1, 0].plot(X[2, :it], linestyle='-', color='k')
            ax1[1, 0].set_xlabel('iteration', fontsize=10)
            ax1[1, 0].set_ylabel('nu (m)', fontsize=10)
            ax1[1, 1].plot(X[3, :it], linestyle='-', color='k')
            ax1[1, 1].set_xlabel('iteration', fontsize=10)
            ax1[1, 1].set_ylabel('density (kg/m3)', fontsize=10)
            plt.tight_layout()
            fig1.savefig(f"{fig_save_path}/SA_inversion_parameters.png", dpi=300, bbox_inches='tight')
            fig0.savefig(f"{fig_save_path}/SA_inversion_misfit.png", dpi=300, bbox_inches='tight')
            plt.close(fig0)
            plt.close(fig1)
            I = np.argmin(misfit_accepted[0:it - 1])    
            print(['best fit for iter = ' + str(I) + ' -> ' + str(X[:, I])])
    I = np.argmin(misfit_accepted[0:it_critic - 1])         
    return I, T_critic, it_critic, X, likelihood, misfit_accepted






# simulated annealing schedule
isave = 200
N_SA = 30000
T0 = 0.5; Tmin = 0.001
T0param = 0.8; Tminparam = 0.05



min_thickness = 20e-2
max_thickness = 40e-2
min_E = 4.5e9
max_E = 6e9
min_nu = 0.28
max_nu = 0.38
min_rho = 600
max_rho = 950



MAX_param0 = np.zeros(4)
MIN_param0 = np.zeros(4)    
MAX_param0[0] = max_thickness
MIN_param0[0] = min_thickness   
MAX_param0[1] = max_E   
MIN_param0[1] = min_E   
MAX_param0[2] = max_nu   
MIN_param0[2] = min_nu  
MAX_param0[3] = max_rho       
MIN_param0[3] = min_rho    
delta_param0 = (MAX_param0 - MIN_param0)
 
X = np.zeros((4,N_SA))
X[0,0] = 0.15
X[1,0] = 2e9
X[2,0] = 0.35
X[3,0] = 917













# #-----------------------------Plotting--------------------------------------------------
# xmin = 0  # Minimum frequency
# xmax = 20  # Maximum frequency
# ymin = 0  # Minimum wavenumber
# ymax = 4  # Maximum wavenumber
# fig, ax = plt.subplots()
# ax.plot(f_mode, k_mode, linestyle='--', color='r', label='Line between points') 
# ax.plot(f_mode, k_QS, linestyle='--', color='b', label='Line between points') 
# ax.set_title('Flexure wave dispersion curve')
# ax.set_xlabel('Frequency (Hz)', fontsize = 20)
# ax.set_ylabel('Wavenumber (rad/m)', fontsize = 20)
# ax.set_xlim(xmin, xmax)
# ax.set_ylim(ymin, ymax)


I,T_critic, it_critic, X, likelihood  ,misfit_accepted   = simulated_annealing(delta_param0,T0,T0param,Tmin,Tminparam,X,MIN_param0,MAX_param0,data,isave)

toc()
print(['best fit for iter = ' + str(I) + ' -> ' + str(X[:, I])])

# %%
