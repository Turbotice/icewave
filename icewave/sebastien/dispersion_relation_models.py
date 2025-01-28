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

#%%

freq = np.linspace(0.2,1,100)
nu = 1e-2 # in m^2/s

           
k_viscous = Lamb_viscous(freq, nu)
k_deep_water = (2*np.pi*freq)**2/9.81


fig,ax = plt.subplots()
ax.loglog(k_deep_water,freq)
ax.loglog(k_viscous,freq,'--')
