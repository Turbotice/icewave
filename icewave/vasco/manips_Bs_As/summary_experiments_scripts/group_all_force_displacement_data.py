#%%
import numpy as np
import os
import pickle
import csv

#%%
disk = 'D:'
force_disp_dict_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"
output_path = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/all_force_displacement_dicts.pkl"

def open_single_experiment_result_dict(date, acq, serie,force_disp_dict_dir=force_disp_dict_dir):
    if serie==None:
        dict_path = force_disp_dict_dir + f'results_dict_force_displacement_{date}_acq{acq}.pkl'
    else:
        dict_path = force_disp_dict_dir + f'results_dict_force_displacement_{date}_serie{serie}_acq{acq}.pkl'
    
    with open(dict_path, 'rb') as file:
        results_dict = pickle.load(file)

    return results_dict

def list_experiments_to_load(force_disp_dict_dir=force_disp_dict_dir):
    files = os.listdir(force_disp_dict_dir)
    dates = []
    acqs = []
    series = []
    for file in files:
        if file.startswith('results_dict_force_displacement_') and file.endswith('.pkl'):
            print(file)
            parts = file.split('_')
            print(parts)
            date = parts[4]
            if 'serie' in file:
                serie = int(parts[5][5:])
                acq = int(parts[6].split('.')[0].replace('acq', ''))
            else:
                serie = None
                acq = int(parts[5].split('.')[0].replace('acq', ''))
            print(date, acq)
            dates.append(date)
            acqs.append(acq)
            series.append(serie)
    return dates, acqs, series

def open_all_experiments_result_dicts(force_disp_dict_dir=force_disp_dict_dir):
    
    dates, acqs, series = list_experiments_to_load(force_disp_dict_dir)
    print(dates)
    print(acqs)
    print(series)
    
    all_results_dicts = []
    for date, acq, serie in zip(dates, acqs, series):
        results_dict = open_single_experiment_result_dict(date, acq, serie, force_disp_dict_dir)
        all_results_dicts.append(results_dict)
    
    return all_results_dicts

#%%

all_results_dicts = open_all_experiments_result_dicts(force_disp_dict_dir=force_disp_dict_dir)

# add the info about what are xdata2fit and ydata2fit
for results_dict in all_results_dicts:
    results_dict['info_fit'] = 'xdata2fit is the force data (newtons), ydata2fit is the displacement data (meters), and the fit is a linear fit to get the slope and the intercept'

# %%
# save all_results_dicts in a pkl file
# output_path defined at the beginning of the script
with open(output_path, 'wb') as file:
    pickle.dump(all_results_dicts, file)

