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
#disk = 'L:'# disk is Elements on adour
disk = 'C:'
date = '0211'

#data_path = f'{disk}/Share_hublot/Data'
data_path = f'C:/Users/Vasco Zanchi/Desktop/Saguenay2024'
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

# %% aller chercher le fichier contenant les infos sur les lignes de front d'onde
dict_all_wavefront_lines_path = f'{daily_drone_data_path}/Results/traitement_vasco/dict_all_wavefront_lines.pkl'

with open(dict_all_wavefront_lines_path, 'rb') as file:
    dict_all_wavefront_lines = pickle.load(file)

#%%
(yind_positions_fractures, xind_positions_fractures) = np.where(fractures_positions_data['binary']==1)

array_tfrac_approx = np.zeros(len(yind_positions_fractures))
for i in range(len(yind_positions_fractures)):
    array_tfrac_approx[i] = fractures_positions_data['tmaxs_kappa2_t'][yind_positions_fractures[i],xind_positions_fractures[i]]

def shortest_distance_to_any_wavefront(xind=np.ndarray, yind=np.ndarray, tfrac=np.ndarray, dict_all_wavefront_lines=dict_all_wavefront_lines, dict_stereo_pivdata=dict_stereo_pivdata):
    t_ind = np.where(dict_stereo_pivdata['t']>=tfrac)[0][0]
    dict_lines = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']

    min_distances = np.zeros(len(dict_lines))
    xshifts = np.zeros(len(dict_lines))
    yshifts = np.zeros(len(dict_lines))
    for i in range(len(dict_lines)):
        xidcs_line = dict_lines['line_'+str(i)]['x_ind']
        yidcs_line = dict_lines['line_'+str(i)]['y_ind']
        
        distances = np.hypot(xidcs_line - xind, yidcs_line - yind)
        agmn_dist = np.argmin(distances)
        min_dist = distances[agmn_dist]

        min_distances[i] = min_dist
        xshifts[i] = (xidcs_line - xind)[agmn_dist]
        yshifts[i] = (yidcs_line - yind)[agmn_dist]
    
    id_closest_line = np.argmin(min_distances)
    closest_distance = min_distances[id_closest_line]
    xshift = xshifts[id_closest_line]
    yshift = yshifts[id_closest_line]


    # ajouter le calcule de "où est" le front par rapport à la fracture    

    return closest_distance, id_closest_line, xshift, yshift

#%%
# test
array_closest_distance = np.zeros(len(array_tfrac_approx))
yshifts = np.zeros(len(array_tfrac_approx))
for i in range(len(array_tfrac_approx)):
    closest_distance, id_closest_line, xshift, yshift = shortest_distance_to_any_wavefront(xind=xind_positions_fractures[i], yind=yind_positions_fractures[i], tfrac=array_tfrac_approx[i])
    array_closest_distance[i] = closest_distance
    yshifts[i] = yshift
%matplotlib qt
plt.style.use('dark_background')
plt.figure()
plt.scatter(xind_positions_fractures, yind_positions_fractures,c=array_closest_distance*np.sign(yshifts), cmap='seismic', vmin=-20, vmax=20)
plt.colorbar()
plt.show()
# %%
# partie angles !!