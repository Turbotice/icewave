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

from find_wave_front import find_wave_fronts_on_image, find_lines, group_lines_to_dict

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


img = peaks2d

lines = find_lines(img, min_points_line=20, tolerance_px=6.0)

print(f"{len(lines)} ligne(s) trouvée(s) :\n")
for i, line in enumerate(lines):
    print(f"  Ligne {i} — {len(line)} point(s) : {line}")

dict_lines = group_lines_to_dict(lines)

plt.figure()
plt.imshow(matrice2d, cmap='gray')
ct=0
for k in dict_lines:
    xvals2plot = dict_lines[k]['x_ind']
    yvals2plot = dict_lines[k]['y_ind']
    a,b = np.polyfit(xvals2plot,yvals2plot, 1)
    plt.plot(xvals2plot, a*xvals2plot+b,'k')
    plt.plot(xvals2plot, yvals2plot, label='line '+str(ct))
    dict_lines[k]['linear_fit'] = {}
    dict_lines[k]['linear_fit']['slope'] = a
    dict_lines[k]['linear_fit']['intercept'] = b
    ct+=1

plt.legend()
plt.show()




# tracé de la courbure "transverse" le long d'une crête de vague

# tout d'abord, on peut calculer le champ de courbures :
# ATTENTION : CE N'EST PAS ENCORE DANS LES BONNES DIMENSIONS POUR L'ESPACE !!!
kappa_x = np.gradient(np.gradient(matrice2d_smoothed,axis=1),axis=1)
kappa_y = np.gradient(np.gradient(matrice2d_smoothed,axis=0),axis=0)

plt.figure()
for k in dict_lines:

    slope = dict_lines[k]['linear_fit']['slope']
    theta = np.arctan(slope) + np.pi/2
    x_ind = dict_lines[k]['x_ind']
    y_ind = dict_lines[k]['y_ind']
    
    kappa_n_prof = np.zeros(len(x_ind))
    for i in range(len(x_ind)):
        kappa_n_prof[i] = np.cos(theta) * kappa_x[y_ind[i],x_ind[i]] + np.sin(theta) * kappa_y[y_ind[i],x_ind[i]]


    dict_lines[k]['kappa_n_profile'] = kappa_n_prof

    plt.plot(np.abs(kappa_n_prof))

# %%
