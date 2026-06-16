#%% imports modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import pickle
import os
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage
from vasco.tools.clickonfigures import get_n_points

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken
from profiles import *

#%% definition des chemins des données
disk = 'L:'# disk is Elements on adour
date = '0211'

data_path = f'{disk}/Share_hublot/Data'
daily_drone_data_path = f'{data_path}/{date}/Drones'

velocity_field_path = f'{daily_drone_data_path}/exact_solution_real_field_stereo_0211_2024_rectangular_grid.h5'

# traitement stéphane : postitions des fractures (obtenue avec divergence du champ de vitesses)
fractures_positions_path = f'{daily_drone_data_path}/Results/fracture_positions.pkl'

#%% Cellule à exécuter une ceule fois (chargement champs vitesses from stereo piv)
#load dict from h5 file

dict_stereo_pivdata = rw.load_dict_from_h5(velocity_field_path)

#%% chargement des positions des fractures détectées
with open(fractures_positions_path, 'rb') as file:
    fractures_positions_data = pickle.load(file)

# %% definitions des variables utiles (champs d'élévation etc.)
dict_stereo_pivdata.keys()

vx = dict_stereo_pivdata['u'][0,:,:,:]
vy = dict_stereo_pivdata['u'][1,:,:,:]
vz = dict_stereo_pivdata['u'][2,:,:,:]

dt = dict_stereo_pivdata['t'][1] - dict_stereo_pivdata['t'][0]

ux = np.cumsum(vx, axis=2)*dt
uy = np.cumsum(vy, axis=2)*dt
uz = np.cumsum(vz, axis=2)*dt

facq_x = dict_stereo_pivdata['SCALE']['facq_x']

#%% chargement de dict_frac
%matplotlib inline

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)
image = vz[:,:,1000]
plt.imshow(image)
plt.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')

if 'dict_frac' in globals():
    for k in dict_frac:
        if 'single_frac' in k:
            plt.plot(dict_frac[k]['idcs_single_frac'][0][0], dict_frac[k]['idcs_single_frac'][0][1], '^r', markersize=4)
plt.show()

#%%
from scipy.ndimage import gaussian_filter1d
from find_wave_front import find_wave_fronts_on_image, find_lines, group_lines_to_dict, compute_kappa_n_profiles_along_wavecrest

k = 'dict_single_frac_yind10_xind39'

idcs_single_frac = dict_frac[k]['idcs_single_frac']
time_frac_approx = dict_frac[k]['times_frac_sec_approx_ref_noncassee'][0]
index_time_frac_approx = np.where(dict_stereo_pivdata['t']>=time_frac_approx)[0][0]


matrice2d = uz[:,:,index_time_frac_approx]

peaks2d, matrice2d_smoothed = find_wave_fronts_on_image(matrice2d, sigma_smooth=5, axis=0, plot=False)

plt.figure()
plt.imshow(matrice2d_smoothed)
plt.imshow(peaks2d)
plt.plot(idcs_single_frac[:,0], idcs_single_frac[:,1], '^k')
plt.show()

#%%


'''
# test compute local angle
plt.figure()
plt.plot(dict_lines['line_0']['x_ind'],gaussian_filter1d(dict_lines['line_0']['y_ind'], sigma=1))

x_test = dict_lines['line_0']['x_ind']
ysm_test = gaussian_filter1d(dict_lines['line_0']['y_ind'],sigma=2)

#tanalpha = np.diff(ysm_test)/np.diff(x_test)
alpha = np.arctan2(np.gradient(x_test), np.gradient(ysm_test))

plt.plot(x_test, alpha)'''
%matplotlib inline
for i in range(0,500,50):
    dict_lines, (peaks2d, matrice2d_smoothed, matrice2d) = compute_kappa_n_profiles_along_wavecrest(uz, index_time=500+i, facq_x=dict_stereo_pivdata['SCALE']['facq_x'])

