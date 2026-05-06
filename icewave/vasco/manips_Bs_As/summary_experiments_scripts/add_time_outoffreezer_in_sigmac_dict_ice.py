#%%
import numpy as np
import pickle
import pandas as pd
import os

#%%
disk = 'D:'
# first step : read the csv file containing temperature info 
csv_path = f'{disk}/manips_BsAs/Summary/time_out_of_freezer.csv'
df = pd.read_csv(csv_path, sep=';')

def convert_arrint2arrstr(arr):
    """
    converts 1d array of int to an array of str
    """
    arr_str = np.empty(len(arr),dtype=object)
    for i in range(len(arr)):
        arr_str[i] = str(arr[i])
    return arr_str

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
    arr_new = np.zeros(len(arr))
    for i in range(len(arr)):
        arr_new[i] = str_to_float(arr[i])
    return arr_new
# %%
date_arr = convert_arrint2arrstr(df['date'].values)
serie_arr = df['serie'].to_numpy()
acq_arr = df['acq'].to_numpy()

import sys
sys.path.append('../tracking_force_displacement/functions')
from utils import correspond_samplenum_acqnum

samplenum_arr = np.empty(len(acq_arr),dtype=object)
for i in range(len(acq_arr)):
    fd, sd = correspond_samplenum_acqnum(date=date_arr[i], acq=acq_arr[i], serie=serie_arr[i])
    samplenum_arr[i] = f'M{fd}{sd}'

time_at_warm_sec_arr = convert_array2floatarray(df['time_at_warm_sec'].to_numpy())



sigmac_dir_path = f'{disk}/manips_BsAs/Summary/sigma_c/'
list_pkl_files = [f for f in os.listdir(sigmac_dir_path) if f.endswith('.pkl')]

print(list_pkl_files)


for pkl_file in list_pkl_files:
    date_in_filename = pkl_file.split('mc_')[1].split('.pkl')[0]
    if date_in_filename in date_arr:
        idx_date = np.where(date_arr == date_in_filename)[0][0]
        

        # open the pkl file, check if 'temperatures_samples' exists, if not create it
        with open(sigmac_dir_path+pkl_file, 'rb') as f:
            data = pickle.load(f)
        
        time_at_warm_sec_tosave = np.zeros(len(data['dict_epaisseurs'].keys()))

        for i in range(len(data['dict_epaisseurs'].keys())):
            samplenum = list(data['dict_epaisseurs'].keys())[i]
            samplenum_arr_date = samplenum_arr[date_arr==date_in_filename]
            if samplenum in samplenum_arr_date:
                idx = np.where((samplenum_arr==samplenum)&(date_arr==date_in_filename))[0][0]
                time_at_warm_sec_tosave[i] = time_at_warm_sec_arr[idx]
            else:
                time_at_warm_sec_tosave[i] = np.nan

        if 'time_at_warm_sec' not in data: # attention la suite ne se fait que à cette condition
            
            data['time_at_warm_sec'] = time_at_warm_sec_tosave

            with open(sigmac_dir_path+pkl_file, 'wb') as f:
                pickle.dump(data, f)
    elif date_in_filename not in date_arr:
        with open(sigmac_dir_path+pkl_file, 'rb') as f:
            data = pickle.load(f)
        time_at_warm_sec_tosave = np.ones(len(data['dict_epaisseurs'].keys())) * np.nan

        if 'time_at_warm_sec' not in data: # attention la suite ne se fait que à cette condition
            
            data['time_at_warm_sec'] = time_at_warm_sec_tosave
            data['time out of freezer not measured'] = True
            with open(sigmac_dir_path+pkl_file, 'wb') as f:
                pickle.dump(data, f)        

#%% check
for i in range(50):
    try:
        # check if it worked for one file
        with open(sigmac_dir_path+list_pkl_files[i], 'rb') as f:
            data = pickle.load(f)
        print(data['time_at_warm_sec'])
    except:
        pass