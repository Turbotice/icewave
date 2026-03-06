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
    'dict_epaisseurs':[],
    'temperatures_samples':[]
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
dict_all['temperatures_samples'] = np.array(dict_all['temperatures_samples'])

dict_all['sigma_c_avg'] = (3/2) * (1e-3* dict_all['mc_avg']*9.81 * dict_all['L'])/(dict_all['w'] * (dict_all['h_avg']**2))
dict_all['sigma_c_std'] = 3 * (dict_all['mc_avg']*1e-3*9.81 * dict_all['L'] * dict_all['h_std'])/(dict_all['w'] * (dict_all['h_avg']**3))


#dict_all['epsilon_c_avg'] = 

#%%
Summary_Gre25_samples_path = f"{disk}/Grenoble/Gre25/Summary/flexion_3pts/Summary_flexion3pts_2025.csv"

dict_all['Gre25_samples'] = {}
df = pd.read_csv(Summary_Gre25_samples_path,sep=';',encoding='latin-1')


dict_all['Gre25_samples']['Date'] = df['Date'].values
dict_all['Gre25_samples']['hour'] = df['hour'].values
dict_all['Gre25_samples']['Sample'] = df['Sample'].values
dict_all['Gre25_samples']['Protocole'] = df['Protocole'].values
dict_all['Gre25_samples']['Thickness'] = df['Thickness'].values
dict_all['Gre25_samples']['Thickness_err'] = df['Thickness_err'].values
dict_all['Gre25_samples']['Force_kg'] = df['Force_kg'].values
dict_all['Gre25_samples']['Total_force_N'] = df['Total_force_N'].values
dict_all['Gre25_samples']['Critical_stress_MPa'] = df['Critical_stress_MPa'].values
dict_all['Gre25_samples']['plot'] = df['plot'].values
dict_all['Gre25_samples']['time2break_sec_approx'] = df['time2break_sec_approx'].values

def str_to_float(value):
    if type(value) == float:
        return value  # Return the original value if it's already a float
    elif value=='':
        return np.nan  # Return NaN for empty strings
    else:
        # Replace comma with dot (European decimal format)
        print(value)
        value = value.replace(',', '.')
        # Remove spaces just in case
        value = value.strip()
        # Convert to float
        print(value)
        return float(value)
def convert_array2floatarray(arr):
    for i in range(len(arr)):
        arr[i] = str_to_float(arr[i])
    return arr

dict_all['Gre25_samples']['Thickness'] = convert_array2floatarray(dict_all['Gre25_samples']['Thickness'])
dict_all['Gre25_samples']['Thickness_err'] = convert_array2floatarray(dict_all['Gre25_samples']['Thickness_err'])
dict_all['Gre25_samples']['Force_kg'] = convert_array2floatarray(dict_all['Gre25_samples']['Force_kg'])
dict_all['Gre25_samples']['Total_force_N'] = convert_array2floatarray(dict_all['Gre25_samples']['Total_force_N'])
dict_all['Gre25_samples']['Critical_stress_MPa'] = convert_array2floatarray(dict_all['Gre25_samples']['Critical_stress_MPa'])
dict_all['Gre25_samples']['time2break_sec_approx'] = convert_array2floatarray(dict_all['Gre25_samples']['time2break_sec_approx'])

mask = dict_all['Gre25_samples']['plot']==1

plt.figure()
plt.plot(dict_all['Gre25_samples']['Thickness'][mask], dict_all['Gre25_samples']['Critical_stress_MPa'][mask], 'o')
plt.xlabel('Thickness (mm)')
plt.ylabel('Critical stress (MPa)')
plt.ylim(0, 5)
plt.xlim(0, 6)
plt.show()

#%%
plt.figure()
plt.plot(dict_all['Gre25_samples']['Thickness'][mask] * 1e-3, dict_all['Gre25_samples']['Critical_stress_MPa'][mask]*1e6, 'o',label='points Gre25')
plt.plot(dict_all['h_avg'], dict_all['sigma_c_avg'], 'o',label='points BsAs')
plt.xlabel('Thickness (m)')
plt.ylabel('Critical stress (Pa)')
plt.ylim(0, 5e6)
plt.xlim(0, 6e-3)
plt.show()
# %%

# plot sigmac vs E
plt.figure()
plt.errorbar(dict_all['E'], dict_all['sigma_c_avg']/dict_all['E'], xerr=dict_all['E_err'], yerr=dict_all['sigma_c_std']/dict_all['E'], linestyle='', marker='o', ecolor='k', color='b')
plt.xlabel('E')
plt.ylabel('$\epsilon_c$')
plt.show()

# plot sigmac vs E
plt.figure()
plt.errorbar(dict_all['E'], dict_all['sigma_c_avg'], yerr=dict_all['sigma_c_std']/dict_all['E'], linestyle='', marker='o', ecolor='k', color='b')
plt.xlabel('E')
plt.ylabel('$\sigma_c$')
plt.show()


epsilon_c_avg = dict_all['sigma_c_avg']/dict_all['E']

# plot epsilonc vs E
plt.figure()
plt.plot(dict_all['h_avg'], epsilon_c_avg, 'o')
plt.xlim(0, 8e-3)
plt.ylim(0,1e-3)
plt.show()

# 
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
plt.plot(epsilon_c_avg/dict_all['time2break_sec'], epsilon_c_avg,'o')
plt.xlabel('$\dot{\epsilon}$ [$s^{-1}$]',fontsize=15)
plt.ylabel('Critical strain $\epsilon_c$',fontsize=15)
plt.show()

plt.figure()
plt.title('Critical stress at rupture vs strain rate')
plt.plot(epsilon_c_avg/dict_all['time2break_sec'], 1e-6*dict_all['sigma_c_avg'],'o')
plt.xlabel('$\dot{\epsilon}$ [$s^{-1}$]',fontsize=15)
plt.ylabel('Critical stress [MPa]',fontsize=15)
plt.show()

plt.figure()
plt.plot((1e3*epsilon_c_avg*dict_all['L']**2)/(6*dict_all['h_avg'])/(dict_all['time2break_sec']) , 1e-9*dict_all['E'],'o')
plt.xlabel('vitesse sollicitation (mm/s)')
plt.ylabel('Youngs modulus (GPa)')

# donnée indicatrice : 
# vitesse d'un échantillon à 0deg posé sur un plot à 18deg : 0.24mm/min (mesuré le 30/09)

plt.figure()
plt.plot((60*1e3*epsilon_c_avg*dict_all['L']**2)/(6*dict_all['h_avg'])/(dict_all['time2break_sec']) , 1e-9*dict_all['E'],'o')
plt.vlines(0.24, 0, 15)
#plt.xlim(0, 13)
#plt.ylim(0,15)
plt.loglog()
plt.xlabel('vitesse sollicitation (mm/min)')
plt.ylabel('Youngs modulus (GPa)')

# plot the equivalent defect size

# fracture toughness ice : between 50 and 150 kPa * m^1/2
K_Ic = 1.5 * 1e5

plt.figure()
plt.plot(dict_all['E'], (1/(2*np.pi))*(K_Ic**2/dict_all['sigma_c_avg']**2),'o')
plt.show()


