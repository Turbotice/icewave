#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
disk = "H"
date = '0223'
# Load the data
data_path = f'{disk}:/data/{date}/Tests_Bresiliens/Summary_samples_tests.csv'

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

indices2plot = np.where(arr_qualite_test == 'ok')[0]
arr_sigma_c = 1e-6 * (2*arr_Fc_kN*1e3)/(np.pi*arr_D_mm*1e-3*arr_L_mm*1e-3)


print(arr_Fc_kN)
plt.figure(figsize=(10, 6))
plt.scatter(arr_L_mm[indices2plot], arr_sigma_c[indices2plot], color='blue', label='Qualité OK')
plt.xlabel('L (mm)')
plt.ylabel('Sigma_c (MPa)')
plt.title('Sigma_c vs L')
# %%
