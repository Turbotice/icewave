#%%
import numpy as np
import os
import pickle
import csv

#%%
disk = 'D:'
sigmac_dict_dir = f"{disk}/manips_BsAs/Summary/sigma_c/"

def open_single_day_result_dict(date, sigmac_dict_dir=sigmac_dict_dir):

    dict_path = sigmac_dict_dir + f'data_h_vs_mc_{date}.pkl'
    
    with open(dict_path, 'rb') as file:
        results_dict = pickle.load(file)

    return results_dict

def list_experiments_to_load(sigmac_dict_dir=sigmac_dict_dir):
    files = os.listdir(sigmac_dict_dir)
    dates = []
    acqs = []
    series = []
    for file in files:
        if file.startswith('data_h_vs_mc_') and file.endswith('.pkl'):
            print(file)
            dates.append(file.split('_')[4].split('.')[0].replace('.pkl',''))
    return dates, acqs, series

def open_all_sigmac_dicts(sigmac_dict_dir=sigmac_dict_dir):
    
    dates, acqs, series = list_experiments_to_load(sigmac_dict_dir)
    print(dates)
    print(acqs)
    print(series)
    
    all_results_dicts = {}
    for date in dates:
        results_dict = open_single_day_result_dict(date, sigmac_dict_dir)
        all_results_dicts[f'date_{date}'] = results_dict    
    return all_results_dicts

#%%

all_sigmac_dicts = open_all_sigmac_dicts(sigmac_dict_dir=sigmac_dict_dir)



# %%

# save all_sigmac_dicts in a pkl file
output_path = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/all_sigmac_dicts.pkl"
with open(output_path, 'wb') as file:
    pickle.dump(all_sigmac_dicts, file)
