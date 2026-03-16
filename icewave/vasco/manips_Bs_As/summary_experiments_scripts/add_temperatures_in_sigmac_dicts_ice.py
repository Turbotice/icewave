#%%
import numpy as np
import pickle
import pandas as pd
import os

#%%
disk = 'D:'
# first step : read the csv file containing temperature info 
csv_path = f'{disk}/manips_BsAs/Summary/temperatures_samples_ice.csv'
df = pd.read_csv(csv_path, sep=';')

def convert_arrint2arrstr(arr):
    """
    converts 1d array of int to an array of str
    """
    arr_str = np.empty(len(arr),dtype=object)
    for i in range(len(arr)):
        arr_str[i] = str(arr[i])
    return arr_str

date_arr = convert_arrint2arrstr(df['date'].values)
T_celsius_arr = df['T_celsius'].values
print(date_arr)
print(T_celsius_arr)

# second step : list all the pkl files in sigmac dir, open them,
#  and add them the corresponding temperature value

sigmac_dir_path = f'{disk}/manips_BsAs/Summary/sigma_c/'
list_pkl_files = [f for f in os.listdir(sigmac_dir_path) if f.endswith('.pkl')]

print(list_pkl_files)


# %%
# for all pkl files, if 'temperatures_samples' exists, don't do anything
# if it doesn't exist, create an array of the same length of the array 
# from the key 'thickness_avg' at the key 'temperatures_samples' and put
# on the same value fir all elements, the value being the temperature value
#  corresponding to the date in the filename
for pkl_file in list_pkl_files:
    date_in_filename = pkl_file.split('mc_')[1].split('.pkl')[0]
    if date_in_filename in date_arr:
        idx_date = np.where(date_arr == date_in_filename)[0][0]
        T_celsius = T_celsius_arr[idx_date]
        # open the pkl file, check if 'temperatures_samples' exists, if not create it
        with open(sigmac_dir_path+pkl_file, 'rb') as f:
            data = pickle.load(f)
        if 'temperatures_samples' not in data:
            length_thickness_avg = len(data['thicknesses_avg'])
            print(pkl_file)
            data['temperatures_samples'] = T_celsius * np.ones(length_thickness_avg)
            with open(sigmac_dir_path+pkl_file, 'wb') as f:
                pickle.dump(data, f)

#%%
# check if it worked for one file
with open(sigmac_dir_path+list_pkl_files[16], 'rb') as f:
    data = pickle.load(f)
print(data['temperatures_samples'])
# %%
