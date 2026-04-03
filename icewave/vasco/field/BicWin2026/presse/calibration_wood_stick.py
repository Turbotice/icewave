#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os, pickle


from load_press_data import load_single_test_curve


#%% Load data for one acquisition
disk = 'H:'
date = '0223'
outputdir = f'{disk}/data/{date}/Tests_Bresiliens/results'

acqnum = 'Wood_Component Compression2 49 04-2-26 14 39 43 PM'
%matplotlib inline
dict_single_test = load_single_test_curve(acqnum=acqnum,plot=True,disk=disk)


# %%
#estimation à la louche
h_mm = 25 # mm
L_mm = 200 
w_mm = 25

h, L, w = 1e-3*h_mm, 1e-3*L_mm, 1e-3*w_mm

Actuator_mm = dict_single_test['Actuator'].values
Load_kN = dict_single_test['Load'].values
Time_sec = dict_single_test['Time'].values

Actuator_m = 1e-3 * Actuator_mm
Load_N = 1e3 * Load_kN

#%% 
def linear_model(x, a, b):
    return a*x+b

Disp_mm_inf_fit = 0.65
Disp_mm_sup_fit = 0.9

indices_fit = np.where((Actuator_mm>=Disp_mm_inf_fit)&(Actuator_mm<=Disp_mm_sup_fit))[0]

sigma_yy = Load_N /(w*L)
epsilon_xx = Actuator_m /h

xdata2fit = epsilon_xx[indices_fit]
ydata2fit = sigma_yy[indices_fit]
popt, pcov = curve_fit(linear_model, xdata2fit, ydata2fit)

plt.figure()
plt.plot(epsilon_xx, sigma_yy)
plt.plot(xdata2fit, linear_model(xdata2fit,popt[0],popt[1]),label='slope='+str(np.round(popt[0]*1e-9,2))+' GPa')
plt.legend()
# %%
