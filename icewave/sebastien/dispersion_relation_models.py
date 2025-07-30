# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:27:35 2025

@author: sebas
"""
import numpy as np
import cmath
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

#%% 


def Lamb_viscous(freq,nu) : 
    
    k = np.linspace(1e-12,5,10000)
    idx_zero  = np.zeros(len(freq))
    for i,f in enumerate(freq):
        omeg = 2*np.pi*f
        
        func = omeg**2 + np.sqrt(nu/2/omeg)*k**2 - 9.81*k
        
        func[func.imag != 0] = -1
        func = func.real
        zero_crossing = np.where(np.diff(np.signbit(func)))[0]
        if zero_crossing.size > 0:
            idx_zero[i] = zero_crossing[0]
        else:
            idx_zero[i] = 0
            
    idx_zero = idx_zero.astype(int)
    k_viscous = k[idx_zero]
    
    return k_viscous


def shallow_water(k,H,g = 9.81):
    """ Shallow water dispersion relation """

    omega = np.sqrt(g*k*np.tanh(k*H))
    return omega

def heavy_water(k,h_ice,H,ratio,g = 9.81):
    """ Shallow water dispersion relation, taking into account sea ice weight"""
    
    omega = np.sqrt(g*k*np.tanh(k*H)/(1 + h_ice*k*ratio))
    return omega

def hydroelastic_squire(k,h_ice,H,E,nu,g = 9.81,rho_w = 1000):
    """ Hydroelastic wave disperion relation, Squire version """
    
    D = E*h_ice**3/12/(1-nu**2)
    omega = np.sqrt((g*k + D*k**5/rho_w)*np.tanh(k*H))
    
    return omega

def heavy_hydroelastic(k,h_ice,H,E,nu,rho_ice,g = 9.81,rho_w = 1000):
    """ Hydroelastic wave dispersion relation, taking into account ice weight"""
    
    D = E*h_ice**3/(12*(1-nu**2))
    ratio = rho_ice/rho_w
    omega = np.sqrt((g*k + D*k**5/rho_w)*np.tanh(k*H)/(1 + h_ice*k*ratio*np.tanh(k*H)))
    
    return omega
    
    
#%% Plot all disperion relations for a given set of parameters

E = 3.5e9 # ice Young modulus
nu = 0.3 # ice Poisson coefficient
rho_ice = 917 # ice density
ratio = rho_ice/1e3
h_ice = 0.1 # ice thickness

H = 100 # water depth 

k = np.linspace(0.1,3,100)
omega_water = shallow_water(k, H)/2/np.pi
omega_heavy_water = heavy_water(k, h_ice, H, ratio)/2/np.pi
omega_squire = hydroelastic_squire(k,h_ice,H,E,nu)/2/np.pi
omega_heavy_squire = heavy_hydroelastic(k, h_ice, H, E, nu, rho_ice)/2/np.pi

fig, ax = plt.subplots()
ax.plot(k,omega_water,label = 'shallow water')
ax.plot(k,omega_squire,label = 'squire')
ax.plot(k,omega_heavy_water,label = 'heavy water')
ax.plot(k,omega_heavy_squire,label = 'heavy squire')


ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$')

ax.set_ylim([0,2]) # set range of frequencies 

ax.legend()




#%%

freq = np.linspace(0.2,1,100)
nu = 1e-2 # in m^2/s

           
k_viscous = Lamb_viscous(freq, nu)
k_deep_water = (2*np.pi*freq)**2/9.81


fig,ax = plt.subplots()
ax.loglog(k_deep_water,freq)
ax.loglog(k_viscous,freq,'--')



#%% Plot tanh(x)




x = np.linspace(0,10)
y = np.tanh(x)

fig, ax = plt.subplots()
ax.plot(x,y)




