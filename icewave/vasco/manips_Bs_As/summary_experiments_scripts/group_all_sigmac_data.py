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


def load_grain_sizes_byhand(date,sample_name,disk=disk):
    path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{date}/{sample_name}_grains_sizes_mm.txt"
    datagrainsizes = np.loadtxt(path2grainsizes) 

    return datagrainsizes, np.mean(datagrainsizes), np.std(datagrainsizes)

def load_all_available_grainsizes(all_sigmac_dicts = dict, sigmac_dict_dir=sigmac_dict_dir):
    dates,_,_ = list_experiments_to_load(sigmac_dict_dir=sigmac_dict_dir)
    for i in range(len(dates)):
        print(i)
        daily_grain_sizes_avg = np.zeros(len(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys()))
        daily_grain_sizes_std = np.zeros(len(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys()))
        print('length list_samples',len(all_sigmac_dicts['date_'+dates[i]]['list_samples']))
        print('length list(keys)', len(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys()))
        for j in range(len(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys())): # boucle for sur toutes les acquisitions
            #print(all_sigmac_dicts['date_'+dates[i]]['list_samples'])
            #print(all_sigmac_dicts['date_'+dates[i]]['list_samples'][j])
            #sample_name = all_sigmac_dicts['date_'+dates[i]]['list_samples'][j][:3]
            print(list(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys())[j])
            sample_name = list(all_sigmac_dicts['date_'+dates[i]]['dict_epaisseurs'].keys())[j]
            path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{dates[i]}/{sample_name}_grains_sizes_mm.txt"
            print(path2grainsizes)
            if os.path.exists(path2grainsizes):
                print(dates[i], sample_name)
                print('exists')
                _, grainsizeavg, grainsizestd = load_grain_sizes_byhand(dates[i], sample_name)
                daily_grain_sizes_avg[j] = grainsizeavg
                daily_grain_sizes_std[j] = grainsizestd
            else:
                daily_grain_sizes_avg[j] = np.nan
                daily_grain_sizes_std[j] = np.nan
            
        all_sigmac_dicts['date_'+dates[i]]['grain_sizes_avg_mm'] = daily_grain_sizes_avg
        all_sigmac_dicts['date_'+dates[i]]['grain_sizes_std_mm'] = daily_grain_sizes_std
        
    
#%%

all_sigmac_dicts = open_all_sigmac_dicts(sigmac_dict_dir=sigmac_dict_dir)

# now add to this dictionary the grain sizes information
load_all_available_grainsizes(all_sigmac_dicts=all_sigmac_dicts)


# %%

# save all_sigmac_dicts in a pkl file
output_path = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/all_sigmac_dicts.pkl"
with open(output_path, 'wb') as file:
    pickle.dump(all_sigmac_dicts, file)
