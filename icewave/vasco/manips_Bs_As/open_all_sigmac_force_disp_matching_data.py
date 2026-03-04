#%%
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import pandas as pd
import csv
#%%
disk = 'D:'
dict_dir = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/"
# load all matching data
matched_data_path = f"{dict_dir}all_sigmac_force_disp_matching_data.pkl"

with open(matched_data_path, 'rb') as f:
    all_sigmac_force_disp_matching_data = pickle.load(f)

#%%
# plotting sigma_c vs Young's modulus for all matched data


dict_all = {
    'date':[],
    'acq':[],
    'serie':[],
    'idx_frac_common':[],
    'prefix_filename_thickness':[],
    'force_newtons_comidx':[],
    'displacement_meters_comidx':[],
    'maxforcefit':[],
    'qualite_fit':[],
    'xdata2fit':[],
    'ydata2fit':[],
    'w':[],
    'L':[],
    'E':[],
    'E_err':[],
    'formula E':[],
    'slope':[],
    'slope_err':[],
    'slope info':[],
    'intercept':[],
    'intercept_err':[],
    'h_avg':[],
    'h_std':[],
    'time2break_sec':[],
    'thicknesses_mm_avg':[],
    'thicknesses_mm_std':[],
    'mc_avg':[],
    'mc_err':[],
    'dict_epaisseurs':[]
    }

for key in all_sigmac_force_disp_matching_data:
    single_sample_dict = all_sigmac_force_disp_matching_data[key]
    for KEY in dict_all:
        if (KEY == 'thicknesses_mm_avg')|(KEY == 'thicknesses_mm_std'):
            NEWKEY = KEY.replace('_mm','')
            dict_all[KEY].append(all_sigmac_force_disp_matching_data[key][NEWKEY])    
        else:
            dict_all[KEY].append(all_sigmac_force_disp_matching_data[key][KEY])

dict_all['mc_avg'] = np.array(dict_all['mc_avg'])
dict_all['mc_err'] = np.array(dict_all['mc_err'])
dict_all['E'] = np.array(dict_all['E'])
dict_all['E_err'] = np.array(dict_all['E_err'])
dict_all['slope'] = np.array(dict_all['slope'])
dict_all['slope_err'] = np.array(dict_all['slope_err'])
dict_all['intercept'] = np.array(dict_all['intercept'])
dict_all['intercept_err'] = np.array(dict_all['intercept_err'])
dict_all['h_avg'] = np.array(dict_all['h_avg'])
dict_all['h_std'] = np.array(dict_all['h_std'])
dict_all['thicknesses_mm_avg'] = np.array(dict_all['thicknesses_mm_avg'])
dict_all['thicknesses_mm_std'] = np.array(dict_all['thicknesses_mm_std'])
dict_all['time2break_sec'] = np.array(dict_all['time2break_sec'])
dict_all['L'] = np.array(dict_all['L'])
dict_all['w'] = np.array(dict_all['w'])

dict_all['sigma_c_avg'] = (3/2) * (1e-3* dict_all['mc_avg']*9.81 * dict_all['L'])/(dict_all['w'] * (dict_all['h_avg']**2))
dict_all['sigma_c_std'] = 3 * (dict_all['mc_avg']*1e-3*9.81 * dict_all['L'] * dict_all['h_std'])/(dict_all['w'] * (dict_all['h_avg']**3))


#dict_all['epsilon_c_avg'] = 

# %%

# plot sigmac vs E
plt.figure()
plt.errorbar(dict_all['E'], dict_all['sigma_c_avg'], yerr=dict_all['sigma_c_std'], linestyle='', marker='o', ecolor='k', color='b')
plt.show()

epsilon_c_avg = dict_all['sigma_c_avg']/dict_all['E']

# plot epsilonc vs E
plt.figure()
plt.plot(dict_all['h_avg'], epsilon_c_avg, 'o')
plt.xlim(0, 8e-3)
plt.ylim(0,1e-3)
plt.show()

# %%
plt.figure()
plt.plot(epsilon_c_avg/dict_all['time2break_sec'], 1e-9*dict_all['E'],'o')
plt.title('Youngs modulus vs strain rate')
plt.xlabel('$\dot{\epsilon}$ [$s^{-1}$]',fontsize=15)
plt.ylabel('Fitted Youngs modulus [GPa]',fontsize=15)
plt.xlim(0,0.0007)
plt.ylim(0,15)
plt.show()

plt.figure()
plt.title('Critical stress at rupture vs strain rate')
plt.plot(epsilon_c_avg/dict_all['time2break_sec'], 1e-6*dict_all['sigma_c_avg'],'o')
plt.xlabel('$\dot{\epsilon}$ [$s^{-1}$]',fontsize=15)
plt.ylabel('Critical stress [MPa]',fontsize=15)
plt.show()
