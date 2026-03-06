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

# load force displacement data
force_disp_dict_path = f"{dict_dir}all_force_displacement_dicts.pkl"
with open(force_disp_dict_path, 'rb') as f:
    all_force_disp_dicts = pickle.load(f)

#load sigma_c data
sigmac_dict_path = f"{dict_dir}all_sigmac_dicts.pkl"
with open(sigmac_dict_path, 'rb') as f:
    all_sigmac_dicts = pickle.load(f)



#%%
matched_data = {}
for i in range(len(all_force_disp_dicts)):
    date2match = all_force_disp_dicts[i]['date']
    prefix2match = all_force_disp_dicts[i]['prefix_filename_thickness']   

    matched_data[f'date_{date2match}_sample{prefix2match}'] = {}
    for key in all_force_disp_dicts[i]:
        matched_data[f'date_{date2match}_sample{prefix2match}'][key] = all_force_disp_dicts[i][key]

    for key in all_sigmac_dicts[f'date_{date2match}']:
        for j in range(len(all_sigmac_dicts[f'date_{date2match}']['thicknesses_avg'])):
            prefix = list(all_sigmac_dicts[f'date_{date2match}']['dict_epaisseurs'].keys())[j]
            print(f'comparing {prefix2match} and {prefix}')
            if prefix==prefix2match:
                for KEY in all_sigmac_dicts[f'date_{date2match}']:
                    print('match')
                    print(KEY)
                    print(all_sigmac_dicts[f'date_{date2match}'][KEY])
                    if KEY=='dict_epaisseurs':
                        matched_data[f'date_{date2match}_sample{prefix2match}'][KEY] = all_sigmac_dicts[f'date_{date2match}'][KEY][prefix2match]
                        print('dict_epaisseurs enregistré')
                    elif (len(all_sigmac_dicts[f'date_{date2match}'][KEY])==len(all_sigmac_dicts[f'date_{date2match}']['thicknesses_avg'])):
                        matched_data[f'date_{date2match}_sample{prefix2match}'][KEY] = all_sigmac_dicts[f'date_{date2match}'][KEY][j]
                    else:
                        print(f'Key {KEY} is not of the right size')

    
# %%
# save the matched data dictionary
matched_data_path = f"{dict_dir}all_sigmac_force_disp_matching_data.pkl"
with open(matched_data_path, 'wb') as f:
    pickle.dump(matched_data, f)



################################################
################################################
################################################
# %% 
