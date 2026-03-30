#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
disk = "H"
date = '0223'
# Load the data
dir_path = f'{disk}:/data/{date}/Tests_Bresiliens/'
data_path = f'{dir_path}Summary_samples_tests.csv'

df = pd.read_csv(data_path, sep=';', encoding='latin-1')
#%%
def convert_to_float(value):
    if type(value) != str:
        return value
    if type(value) == str and value!='A NOTER':
        return float(value.replace(',', '.'))
    else:
        return np.nan
def convert_to_arrayoffloats(values):
    return np.array([convert_to_float(value) for value in values])


# Display the first few rows of the DataFrame
print(df.head())

arr_L_mm = df['L_mm'].values
arr_L_err_mm = df['L_err_mm'].values
arr_qualite_test = df['qualite_test'].values
arr_Fc_kN = convert_to_arrayoffloats(df['Fc_kN'].values)
arr_D_mm = df['D_mm'].values
arr_T_core_celsius = convert_to_arrayoffloats(df['T_core_celsius'].values)
arr_mass_kg = 1e-3 * df['mass_g'].values
arr_loc_core = df['loc_core'].values

#indices2plot = np.where((arr_qualite_test == 'ok')|(arr_qualite_test=='perfect'))[0]
mask_quality_perfect = (arr_qualite_test=='perfect')
mask_quality_ok = np.array(['ok' in x for x in arr_qualite_test])
mask_quality = mask_quality_ok|mask_quality_perfect

arr_sigma_c = 1e-6 * (2*arr_Fc_kN*1e3)/(np.pi*arr_D_mm*1e-3*arr_L_mm*1e-3)


print(arr_Fc_kN)
plt.figure(figsize=(10, 6))
plt.scatter(arr_L_mm[mask_quality], arr_sigma_c[mask_quality], color='blue', label='Qualité OK')
plt.xlabel('L (mm)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(0,np.max(arr_L_mm[mask_quality])*1.1)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs L')
plt.savefig(f'{dir_path}sigmac_vs_L.pdf', dpi=300)
plt.show()


plt.figure(figsize=(10, 6))
plt.scatter(arr_T_core_celsius[mask_quality], arr_sigma_c[mask_quality], color='orange', label='Qualité OK')
plt.xlabel('T (°C)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(-5, 0)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs T')
plt.savefig(f'{dir_path}sigmac_vs_T.pdf', dpi=300)
plt.show()

# %%


mask_HaHa = np.array(['HaHa' in x for x in arr_loc_core])
mask_Hatee = np.array(['Hatee' in x for x in arr_loc_core])
mask_Capelans_beach = np.array(['Capelans beach' in x for x in arr_loc_core])
mask_Capelans_floating_ice = np.array(['Capelans floating' in x for x in arr_loc_core])


plt.figure(figsize=(10, 6))
plt.scatter(arr_T_core_celsius[mask_quality&mask_HaHa], arr_sigma_c[mask_quality&mask_HaHa], color='tab:blue', label='Ha!Ha!')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Hatee], arr_sigma_c[mask_quality&mask_Hatee], color='tab:orange', label='Hatée')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Capelans_beach], arr_sigma_c[mask_quality&mask_Capelans_beach], color='tab:green', label='Capelans beach')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Capelans_floating_ice], arr_sigma_c[mask_quality&mask_Capelans_floating_ice], color='tab:red', label='Capelans floating ice')

plt.xlabel('T (°C)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(-5, 0)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs T')
plt.savefig(f'{dir_path}sigmac_vs_T.pdf', dpi=300)
plt.legend()
plt.show()
