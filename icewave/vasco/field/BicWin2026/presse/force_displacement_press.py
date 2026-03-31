#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os, pickle


from load_press_data import load_single_test_curve
sys.path.append('../../../tools')
from clickonfigures import get_n_points



#%% Load data for one acquisition
disk = 'H:'
date = '0223'
outputdir = f'{disk}/data/{date}/Tests_Bresiliens/results'


acqnum = '1-B'
%matplotlib inline
dict_single_test = load_single_test_curve(acqnum=acqnum,plot=True,disk=disk)

Actuator_mm = dict_single_test['Actuator'].values
Load_kN = dict_single_test['Load'].values
Time_sec = dict_single_test['Time'].values

Actuator_m = 1e-3 * Actuator_mm
Load_N = 1e3 * Load_kN

#%% click to get the limits between which to linfit

n_points = int(input('numbers of clicks ?'))
%matplotlib qt
coords = get_n_points(Actuator_mm, Load_kN, n_points=n_points)

#%% store coords limits in dict for linear fits

interp_LoadkN_vs_Actuatormm = interp1d(Actuator_mm, Load_kN)

dict_coords = {}
dict_coords['x'] = {}
dict_coords['y_interp'] = {}

for i in range(int(n_points)):
    if i%2==0:
        x = coords[i][0]
        dict_coords['x'][f'x{(i+1+1)//2}inf'] = x
        dict_coords['y_interp'][f'y{(i+1+1)//2}inf'] = float(interp_LoadkN_vs_Actuatormm(x))

    elif i%2==1:
        x = coords[i][0]
        dict_coords['x'][f'x{(i+1+1)//2}sup'] = x
        dict_coords['y_interp'][f'y{(i+1+1)//2}sup'] = float(interp_LoadkN_vs_Actuatormm(x))


# %% linear fits between the limits and plots

def linear_model(x, a, b):
    return a*x+b

plt.figure()
# 🔹 Plot raw data once
plt.scatter(Actuator_mm, Load_kN, s=10, label='Data')

dict_fits = {}

dict_fits['Actuator_mm'] = Actuator_mm
dict_fits['Load_kN'] = Load_kN

for key in dict_coords['x']:
    if 'inf' in key:
        idx = key.replace('x', '').replace('inf', '')
        
        x_inf = dict_coords['x'][f'x{idx}inf']
        x_sup = dict_coords['x'][f'x{idx}sup']
        
        x_min, x_max = sorted([x_inf, x_sup])
        
        # Select interval data
        mask = (Actuator_mm >= x_min) & (Actuator_mm <= x_max)
        x_data = Actuator_mm[mask]
        y_data = Load_kN[mask]
        
        if len(x_data) > 1:
            try:
                # 🔹 Fit
                popt, pcov = curve_fit(linear_model, x_data, y_data)
                a, b = popt
                
                dict_fits[f'fit_{idx}'] = {'a': a, 'b': b, 'cov': pcov}
                
                # 🔹 Plot fit immediately
                x_fit = np.linspace(x_min, x_max, 50)
                y_fit = linear_model(x_fit, a, b)
                
                dict_fits[f'fit_{idx}']['x_data'] = x_data
                dict_fits[f'fit_{idx}']['y_data'] = y_data
                
                plt.plot(x_fit, y_fit,color='tab:orange')
                
            except RuntimeError:
                continue

# 🔹 Cosmetics
plt.xlabel('Actuator (mm)')
plt.ylabel('Load (kN)')
plt.grid()

plt.show()

#%% save in dictionary of all acqs for this acqnum value


os.makedirs("results", exist_ok=True)
fp = f"{outputdir}/fits_force_displacement.pkl"

all_results = pickle.load(open(fp, "rb")) if os.path.exists(fp) else {}
all_results[acqnum] = dict_fits
ecr = input('Ecraser les données précédentes pour cette acquisition? (y/n)')
if ecr=='y':
    pickle.dump(all_results, open(fp, "wb"))
else:
    pass
# %%
