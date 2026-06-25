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
from vasco.tools.tridimfits import fit_3dparabola, measure_bidimensional_curvature_around_central_point, convertir_base_tournee, change_order_from_polyfeatures_to_polycoefs


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
from find_wave_front import find_wave_fronts_on_image, find_lines, group_lines_to_dict, compute_kappa_n_profiles_along_wavecrest, bidim_curvature_oneline

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

#%% Find lines for all frames and save pkl file
dic_all_lines = {}
for i in range(vz.shape[2]):
    dic_all_lines['frame_'+str(i)] = {}
    dict_lines, (peaks2d, matrice2d_smoothed, matrice2d) = compute_kappa_n_profiles_along_wavecrest(uz, index_time=i, facq_x=dict_stereo_pivdata['SCALE']['facq_x'], plot=False)
    dic_all_lines['frame_'+str(i)]['dict_lines'] = dict_lines
    dic_all_lines['frame_'+str(i)]['other_infos'] = (peaks2d, matrice2d_smoothed, matrice2d)

pkl_filename = 'dict_all_wavefront_lines'
pkl_filepath = f'{daily_drone_data_path}/Results/traitement_vasco/{pkl_filename}.pkl'
save_input = input('save ? (y/n)')
if save_input=='y':
    pickle.dump(dic_all_lines, open(pkl_filepath, "wb"))
else:
    pass

# %% plot histogramme des temps de fractures
(yind_positions_fractures, xind_positions_fractures) = np.where(fractures_positions_data['binary']==1)

array_tfrac_approx = np.zeros(len(yind_positions_fractures))
for i in range(len(yind_positions_fractures)):
    array_tfrac_approx[i] = fractures_positions_data['tmaxs_kappa2_t'][yind_positions_fractures[i],xind_positions_fractures[i]]

%matplotlib inline
hh = plt.hist(array_tfrac_approx, bins=40)
xvals = hh[0]
yvals = hh[1]
plt.xlim(0, np.max(xvals)*1.05)
plt.ylim(0, np.max(yvals)*1.05)
plt.xlabel('Fracture time [sec]', fontsize=15)
plt.ylabel('# of pixels', fontsize=15)
plt.show()

# et plot t_frac_approx vs yind
%matplotlib inline
plt.figure()
plt.plot(yind_positions_fractures,array_tfrac_approx,'o')
plt.show()

#%%
# on va chercher kappa_n et kappa_t 
# pour les 3 premieres fractures d'intérêt

# on va faire "à moitié" à  la main :
# sélection de la ligne de front d'onde 
# qui nous intéresse, après avoir choisi la fracture 
# qui nous intéresse et donc le temps (frame) à afficher

k1 = 'dict_single_frac_yind23_xind28'
tfrac_approx = dict_frac[k1]['times_frac_sec_approx_ref_noncassee'][0]
ind_tfrac_approx = np.where(dict_stereo_pivdata['t']>=tfrac_approx)[0][0]

plt.figure()
plt.imshow(uz[:,:,ind_tfrac_approx])
plt.plot(dict_frac[k1]['idcs_single_frac'][:,0], dict_frac[k1]['idcs_single_frac'][:,1], 'k^')
for i in range(len(dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines'])):
    xvals = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(i)]['x_ind']
    yvals = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(i)]['y_ind']
    plt.plot(xvals, yvals, label='line '+str(i))

plt.legend()
plt.show()

# choix par exemple ici : line 0
line_index = 0
# on peut donc choisir dans le dic_all_lines la frame 
# qui nous intéresse et la ligne qui nous intéresse
kappa_n_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_n_profile']
kappa_t_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_t_profile']

plt.figure()
plt.plot(kappa_n_profile)
plt.plot(kappa_t_profile)



arr_coefficients_invariantbasis, arr_coefficients_movingbasis, alpha_arr = bidim_curvature_oneline(uz=uz, 
                                                                                                   dic_all_lines=dic_all_lines,
                                                                                                   ind_t=ind_tfrac_approx, 
                                                                                                   line_index=line_index,
                                                                                                   window_size=(11,11))
#####################
kappa_xx = 2 * arr_coefficients_invariantbasis[:,0]
kappa_yy = 2 * arr_coefficients_invariantbasis[:,1]


#####################
kappa_tt_profile = 2 * arr_coefficients_movingbasis[:,0]
kappa_nn_profile = 2 * arr_coefficients_movingbasis[:,1]
K_gauss_movingbasis = arr_coefficients_movingbasis[:,0]*arr_coefficients_movingbasis[:,1] - arr_coefficients_movingbasis[:,2]**2
