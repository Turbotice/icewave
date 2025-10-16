#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:07:32 2025

@author: moreaul
"""


import os
import matplotlib.pyplot as plt
import gpxpy
import pickle
import numpy as np

import icewave.tools.rw_data as rw




date = '0211'
year = '2025'
# path2data = f'/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/{year}_BICWIN/{date}/DAS'
path2data = f'U:/Data/{date}/DAS/'

save_path = os.path.join(path2data, f'disp_curv_{date}.pkl')


plt.close("all")
with open(save_path, 'rb') as f:
    data = pickle.load(f)

freq_dict = data['freq']      # now a dict
kQS_dict = data['k_QS']       # dict
xpos_dict = data['xposition'] # dict

plt.figure(figsize=(8,6))

# On normalise les positions pour les mapper à une colormap
xvals = np.array(list(xpos_dict.values()))
norm = plt.Normalize(vmin=xvals.min(), vmax=xvals.max())
cmap = plt.cm.viridis   # tu peux choisir 'plasma', 'inferno', 'turbo', etc.


plt.figure(figsize=(8,6))
for key in freq_dict.keys():
    f_mode = freq_dict[key]
    k_mode = kQS_dict[key]
    xpos = xpos_dict[key]
    color = cmap(norm(xpos))
    plt.plot(f_mode, k_mode, color=color, label=f"x={xpos:.1f} m")

plt.xlabel("Frequency [Hz]")
plt.ylabel("Wavenumber [1/m]")
plt.legend()  # affiche toutes les légendes
plt.show()


#%% Save data dictionnary in a h5 file 

new_data = {}
main_keys = ['freq','k_QS','xposition']
for key in main_keys :
    new_data[key] = {}
    for sub_key in freq_dict.keys():
        
        new_key = str(sub_key)
        new_data[key][new_key] = data[key][sub_key]

save_path = os.path.join(path2data, f'new_disp_curv_{date}.h5')
rw.save_dict_to_h5(new_data, save_path)
