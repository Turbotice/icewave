#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:06:08 2024

@author: moreaul
"""


import numpy as np
from random import random
import matplotlib.pyplot as plt 
import pickle
import os , sys
import time
from scipy.interpolate import interp1d
import seaborn as sns


#%%-------------------------------LOADING DATA--------------------------------------------------------------------------
plt.close('all')

year = '2025'
date = '0206' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 
equation = 'stein'
# path to dispersion relation data 
path2data = os.path.join('E:/Data/',date,'Geophones/')
# path2data = 'C:/Users/sebas/icewave/icewave/sebastien/geophones/updatescriptspython/0211/Geophones/'
# path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/' +year+'_BICWIN/'

# Load (f,k) points of QS dispersion relation 
direction = 1 # 1 ou 2 
filename = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
# file2load1 = path2data  + date + '/Geophones/' + filename
file2load1  = os.path.join(path2data,filename)
with open(file2load1, "rb") as f:
    data = pickle.load(f)
freq1 = data[0] 
kQS1 = data[1] 

# Load (f,k) points of QS dispersion relation
direction = 2 # 1 ou 2 
filename = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
# file2load1 = path2data  + date + '/Geophones/' + filename
file2load1  = os.path.join(path2data,filename)
with open(file2load1, "rb") as f:
    data = pickle.load(f)
freq2 = data[0] 
kQS2 = data[1] 


mini = max(min(freq1),min(freq2))
maxi = min(max(freq1),max(freq2))
freq = np.linspace(mini,maxi , 15)

# Interpolate kQS1 at the new frequency points, without extrapolation
interpolator1 = interp1d(freq1, kQS1)
kQS1_interp = interpolator1(freq)

# Interpolate kQS2 at the new frequency points, without extrapolation
interpolator2 = interp1d(freq2, kQS2)
kQS2_interp = interpolator2(freq)

# Calculate the average values, ignoring NaNs
kQS = np.nanmean([kQS1_interp, kQS2_interp], axis=0)

# Plot the results
plt.figure(figsize=(12, 6))
plt.plot(kQS1, freq1,  'bo', label='kQS1')
plt.plot( kQS2,freq2, 'go', label='kQS2')
plt.plot(kQS1_interp, freq,  'b--', label='kQS1 interpolated')
plt.plot(kQS2_interp, freq, 'g--', label='kQS2 interpolated')
plt.plot(kQS, freq, 'r--', label='Average kQS')
plt.xlabel(r'$k_{QS} \: \mathrm{(rad.m^{-1})}$')
plt.ylabel(r'$f \: \mathrm{rad.m^{(-1})}$')
plt.legend()
plt.title('Comparison and Average of kQS1 and kQS2')
plt.show()

#%% Load data for inversion 
# file that will be saved 
file2save  = path2data + year + '_' + date + '_acq'+acqu_numb+'_bidir_inversion.pkl'
fig_save_path = path2data + '/Figures_inversion_MCMC_T_0p25/' 
if not os.path.isdir(fig_save_path):
    os.mkdir(fig_save_path)

file2load2 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir.pkl'
file2load3 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir.pkl'


# Load phase velocity of in-plane waves 
with open(file2load2, "rb") as f:
    data = pickle.load(f)
cQS0 = data 
with open(file2load3, "rb") as f:
    data = pickle.load(f)
cSH0 = data 

# load dictionnary of phase velocity 
# file2load2 = path2data +'Phase_velocity_dictionnary_acqu_0001_sig_length_0p3.pkl'
# with open(file2load2,'rb') as pfile:
#     dic_phase_velocity = pickle.load(pfile)
# print('Phase velocities loaded')

data = {
    'freq':freq,
    'kQS': kQS,
    'cQS0': cQS0,
    'cSH0': cSH0
}

#%% Find a first guess of E and nu, assuming a given density for ice 

rho_ice = 917
nu = 1-2*(data['cSH0']/data['cQS0'])**2
E = rho_ice*data['cQS0']**2*(1-nu**2)

print(f'Young modulus, E = {E*1e-9} and Poisson coefficient, nu = {nu}')

#########################################################################
#%%----------------------- PROCEED INVERSION ---------------------------
#########################################################################

#--------------------------------------------------------------------------------------------------------------------------

def tic():
    global start_time
    start_time = time.perf_counter()

def toc():
    elapsed_time = time.perf_counter() - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")



#----------------------------------FUNCTIONS DEFINITION---------------------------------------------------------------------

tic()


def candidate(X,delta_X,min_X,max_X): 
    """ Generate an array of random values between a given range of values. The accessible range is defined by 
    a minimal and a maximal value, but can also be restrained to a smaller range
    Inputs :
        - X : array, generation of candidate will be done along axis = 0
        - delta_X : size of step within which we are generating a new candidate. 
        Corresponds to step in Metropolis algorithm 
        - min_X : array, minimal values of each candidate
        - max_X : array, maximal values of each candidate
    Output :
        - cand : array of random values, where min_X[k] < cand[k] < max_X[k] """
        
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
        # random candidate generation in multispace delta_X box around position  
        cand[k] = bound_min + (bound_max-bound_min)*random()

    return cand    

def wavenumbers_stein_squire( rho_ice, h, H, E, nu,freq,c_w,rho_w,equation):
    """ This function computes the wave vectors associated to a given array of frequencies
    Inputs : 
        - rho_ice : ice density 
        - h : thickness of ice 
        - H : water depth
        - E : Young modulus of ice
        - nu : Poisson coefficient of ice
        - freq : an array of frequencies, to which will correspond wave vectors 
        - c_w : sound waves phase velocity in water
        - rho_w : water density 
        - equation: string, 'stein' : compute wavevectors using stein dispersion relation
                            'squire' : compute wavevectors using squire dispersion relation 
        
    The function returns : 
        - k_QS : wave vectors of the flexural mode
        """
    g = 9.81
    # G = E/(2*(1+nu))
    # cS0 = np.sqrt(E/(rho_ice*(1-nu**2)))
    # cSH0 = np.sqrt(G/rho_ice)
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus 

    k = np.linspace(1e-12,8,200000)
    
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
            a = np.where(np.diff(np.signbit(func)))[0] # indices of the array k at which func(k) = 0
            idx_zero[kf] = a
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero]       
    if flag:
        k_QS[idx_flag] = 0
        
    # k_QS0 = freq/cS0*2*np.pi
    # k_SH0 = freq/cSH0*2*np.pi  
    # cphQS = freq/k_QS*2*np.pi

    
    return k_QS



def simulated_annealing(delta_param0, T0, T0param, Tmin, Tminparam, X, MIN_param0, MAX_param0, data, isave):
    """ Find best parameters using Monte-Carlo-Markof-Chains inversion (MCMC inversion), Metropolis-Hasting algorithm 
    Inputs :
        - delta_param0 : array, initial uncertainty over parameters (h,E,nu,rho)
        - T0 :float, initial temperature at which we perform Metropolis algorithm
        - T0param : float, initial temperature used to determine the step used for Metropolis algorithm
        - Tmin : float, minimal temperature at which we perform Metropolis algorithm
        - Tminparam : float, minimal temperature used to determine the step used for Metropolis algorithm
        - X : array of parameters to fit 
        - MIN_param0 : initial minimal value of each parameter
        - MAX_param0 : initial maximal value of each parameter
        - data : dictionary of data extracted from geophones data, with relevant keys : 'freq','kQS','cQS0','cSH0'
        - isave : int, save generated values every isave iterations
    Outputs : 
        - I : index of the best iteration
        - T_critic : critical temperature 
        - it_critic : maximal iteration reached 
        - X : array of parameters generated at each iteration 
        - likelihood : likelihood of each set of parameters generated
        - misfit_accepted : error to experimental values of kQS, cQS0 and cQSH0 for each generated set of parameters
        """
        
    H = 1 # water depth
    c_w = 1450 # sound waves velocity in water
    rho_w = 1027 # water density 
    freq = data['freq']
    kQSdata = data['kQS']
    cQS0data = data['cQS0']
    cSH0data = data['cSH0']
    N_SA = X.shape[1]
    T = np.exp(np.linspace(np.log(T0), np.log(Tmin), N_SA)) # temperatures used for 'recuit simulé'
    variance = T**2
    Tparam = np.exp(np.linspace(np.log(T0param), np.log(Tminparam), N_SA))
    T_critic = Tmin
    it_critic = N_SA # maximum number of iteration for 'recuit simulé'
    misfit_accepted = np.zeros(N_SA)
    likelihood = np.zeros(N_SA)
    G = X[1, 0] / (2 * (1 + X[2, 0])) # shear modulus
    cQS0synthetics = np.sqrt(X[1, 0] / (X[3, 0] * (1 - X[2, 0]**2))) # phase velocity of acoustic mode
    cSH0synthetics = np.sqrt(G / X[3, 0]) # phase velocity of shear mode 
    # compute wavenumber of flexural mode
    kQSsynthetics = wavenumbers_stein_squire(X[3, 0], X[0, 0], H, X[1, 0], X[2, 0], freq, c_w, rho_w, equation)
    misfit_accepted[0] = 1/3 * (np.mean(np.abs(kQSdata - kQSsynthetics) / kQSdata) + np.mean(np.abs(cQS0data - cQS0synthetics) / cQS0data) 
                                + np.mean(np.abs(cSH0data - cSH0synthetics) / cSH0data))
    likelihood[0] = np.exp(-misfit_accepted[0]**2 / (2 * variance[0])) # gaussian likelihood 
    it = 1 
    it_tot = 0 # iteration counter 
    cpt_out = 0
    while it <= N_SA - 1: 
        it_tot += 1
        delta_param = delta_param0 * Tparam[it]
        success = False
        while not success:
            try:
                Xtest = candidate(X[:, it - 1], delta_param, MIN_param0, MAX_param0) # generate a set of parameters
                G = Xtest[1] / (2 * (1 + Xtest[2]))
                cQS0synthetics = np.sqrt(Xtest[1] / (Xtest[3] * (1 - Xtest[2]**2)))
                cSH0synthetics = np.sqrt(G / Xtest[3])
                kQSsynthetics = wavenumbers_stein_squire(Xtest[3], Xtest[0], H, Xtest[1], Xtest[2], freq, c_w, rho_w, equation)
                success = True
            except Exception as e:
                print(f"Error: {e}. Retrying...")
        misfit_cand = 1/3 * (np.mean(np.abs(kQSdata - kQSsynthetics) / kQSdata) + np.mean(np.abs(cQS0data - cQS0synthetics) / cQS0data) 
                             + np.mean(np.abs(cSH0data - cSH0synthetics) / cSH0data))
        likelihood_cand = np.exp(-misfit_cand**2 / (2 * variance[it])) # likelihood of the new candidate 
        # compute acceptance factor 
        test_likelihood = min(1, likelihood_cand / likelihood[it - 1])
        
        if random() < test_likelihood: # keep candidate if acceptance factor is higher than a random number (Metropolis-Hasting)
            X[:, it] = Xtest # store new candidate
            misfit_accepted[it] = misfit_cand 
            likelihood[it] = likelihood_cand
            cpt_out = 0
            it += 1
            print(it)
        else: # acceptance factor too small -> new candidate has a smaller likelihood
            cpt_out += 1
            # If we have more than 500 consecutive acceptance test, we consider we found critical temperature
            if cpt_out > 500: # 500 failed acceptance test
                I = np.argmin(misfit_accepted[0:it])
                T_critic = T[it]
                it_critic = it # critical iteration
                with open(file2save, 'wb') as pfile:
                    pickle.dump([T_critic, X, misfit_accepted], pfile)    
                break
            
        if it % isave == 0 and cpt_out < 1:
            
            
            with open(file2save, 'wb') as pfile:
                pickle.dump([T_critic, X, misfit_accepted], pfile)
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
            fig1.savefig(f"{fig_save_path}_{acqu_numb}_dir2_SA_inversion_parameters.png", dpi=300, bbox_inches='tight')
            fig0.savefig(f"{fig_save_path}_{acqu_numb}_dir2_SA_inversion_misfit.png", dpi=300, bbox_inches='tight')
            plt.close(fig0)
            plt.close(fig1)
            I = np.argmin(misfit_accepted[0:it - 1])    
            print(['best fit for iter = ' + str(I) + ' -> ' + str(X[:, I])])
            
            
            nonzeroidx = np.nonzero(X[0, :])[0][-1]
            data_to_plot = [X[i, :nonzeroidx] for i in range(4)]
            fig2, axes2 = plt.subplots(2, 2, figsize=(12, 10))
            titles = ['Thickness (m)', 'E (GPa)',
                      'nu ', 'rho (kg/m3)']
            for i, ax in enumerate(axes2.flat):
                ax.hist(data_to_plot[i], bins=30, edgecolor='black')
                ax.set_title(titles[i])
                ax.set_xlabel('Value')
                ax.set_ylabel('Frequency')
            plt.tight_layout()
            figname = fig_save_path + year + '_' + date + '_acq'+acqu_numb+ '_MCMCinversion_dir2.png'
            plt.savefig(figname)
            plt.close(fig2)
            
            
    I = np.argmin(misfit_accepted[0:it_critic - 1])   # best iteration with minimal error to experimental values     
    return I, T_critic, it_critic, X, likelihood, misfit_accepted


#-------------------------INVERSION PARAMETERS-------------------------------------------------------------------------



# simulated annealing schedule
isave = 200 # save data every isave iterations
N_SA = 5000 # maximum number of iterations # initial value = 30 000
T0 = 0.025; Tmin = 0.025 # range of temperature we are using 
#T0param = 0.8; Tminparam = 0.05
T0param = 0.25; Tminparam = 0.25 # minimal and maximal range of the step for Metropolis algorithm 
# If T0param != Tminparam : we are performing a Metropolis algorithm with a progressive step, we can make smaller step 
# after each iteration 


# extremal values of the fitted parameters
min_thickness = 0.2
max_thickness = 1.5
min_E = 1.5e9
max_E = 8e9
min_nu = 0.15
max_nu = 0.45
min_rho = 600
max_rho = 1000


# Initial values of each parameter 
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
 
# Initial values of the different parameters 
X = np.zeros((4,N_SA))
X[0,0] = 0.4 # ice thickness
X[1,0] = 2.5e9 # Young modulus
X[2,0] = 0.3 # Poisson coefficient 
X[3,0] = 900 # ice density 




# #-----------------------------Plotting--------------------------------------------------


print('Start MCMC algorithm...')
I,T_critic, it_critic, X, likelihood  ,misfit_accepted   = simulated_annealing(delta_param0,T0,T0param,Tmin,Tminparam,X,MIN_param0,
                                                                               MAX_param0,data,isave)

toc()
print('MCMC algorithm succeeded')
print(['best fit for iter = ' + str(I) + ' -> ' + str(X[:, I])])



#%%
with open(file2save, 'rb') as file:
    data = pickle.load(file)

T_critic = data[0]
X = data[1]
misfit_accepted = data[2]

nonzeroidx = np.nonzero(X[0, :])[0][-1]

data_to_plot = [X[i, :nonzeroidx] for i in range(4)]

# Step 3: Create a 2x2 subplot and plot histograms with KDE
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

titles = ['Thickness (m)', 'E (Pa)',
          'nu ', 'rho (kg/m3)']

for i, ax in enumerate(axes.flat):
    data = data_to_plot[i]
    
    # Plot the histogram
    sns.histplot(data, bins=30, kde=True, ax=ax, edgecolor='black', color='blue', alpha=0.6)

    # Set titles and labels
    ax.set_title(titles[i])
    ax.set_xlabel('Value')
    ax.set_ylabel('PDF')

plt.tight_layout()
plt.show()

figname = fig_save_path + year + '_' + date + '_acq'+ acqu_numb + '_final_MCMCinversion_dir2'
plt.savefig(figname + '.pdf', dpi = 800, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = 800, bbox_inches = 'tight')

#%% Try to find maximum of gaussian kde
from scipy import stats


kde = stats.gaussian_kde(data_to_plot[0])
Y = kde.evaluate(data_to_plot[0])

# compute argmax
test = np.argmax(Y)



fig, ax = plt.subplots(figsize = (12,10))
ax.plot(data_to_plot[0],Y)