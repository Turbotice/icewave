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
from find_wave_front import find_wave_fronts_on_image, find_lines, group_lines_to_dict, compute_kappa_n_profiles_along_wavecrest, bidim_curvature_oneline, closest_points_between_2curves, Show_two_times_wavefront_and_frac

#k = 'dict_single_frac_yind10_xind39'
k = 'dict_single_frac_yind23_xind28'

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
dict_lines, (peaks2d, matrice2d_smoothed, matrice2d) = compute_kappa_n_profiles_along_wavecrest(uz, index_time=index_time_frac_approx, facq_x=dict_stereo_pivdata['SCALE']['facq_x'])

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

#%% étude "complète" pour 1 crack
# on va chercher kappa_n et kappa_t 
# pour les 3 premieres fractures d'intérêt

# on va faire "à moitié" à  la main :
# sélection de la ligne de front d'onde 
# qui nous intéresse, après avoir choisi la fracture 
# qui nous intéresse et donc le temps (frame) à afficher

# Les 3 fractures qui sont sur le front d'onde au moment de leur détection

#keyfrac = 'dict_single_frac_yind23_xind28'
#line_index = 0
#keyfrac = 'dict_single_frac_yind10_xind39'
#line_index = 0
keyfrac = 'dict_single_frac_yind34_xind25'
line_index = 1

# Les autres, détection en retard :

#keyfrac = 'dict_single_frac_yind56_xind14'
#line_index = 1
#keyfrac = 'dict_single_frac_yind50_xind19'
#line_index = 1
#keyfrac = 'dict_single_frac_yind80_xind9'
#line_index = 4


indt_distmin = Show_two_times_wavefront_and_frac(keyfrac, 
                                                 delta_indt=120,
                                                   uz=uz,
                                                   method='last',
                                                     dic_all_lines=dic_all_lines,
                                                       dict_frac=dict_frac,
                                                         dict_stereo_pivdata=dict_stereo_pivdata)

ind_tfrac_approx = indt_distmin

# on peut donc choisir dans le dic_all_lines la frame 
# qui nous intéresse et la ligne qui nous intéresse
kappa_x_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_x_profile']
kappa_y_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_y_profile']
kappa_n_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_n_profile']
kappa_t_profile = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['kappa_t_profile']

plt.figure()
plt.plot(kappa_n_profile)
plt.plot(kappa_t_profile)


#%%


arr_coefficients_invariantbasis, arr_coefficients_movingbasis, alpha_arr = bidim_curvature_oneline(uz=uz, 
                                                                                                   dic_all_lines=dic_all_lines,
                                                                                                   ind_t=ind_tfrac_approx, 
                                                                                                   line_index=line_index,
                                                                                                   window_size=(11,11))

#%%
#####################
# extraction des courbures (converties en unités physiques)
facq_x = dict_stereo_pivdata['SCALE']['facq_x']

kappa_xx = 2 * arr_coefficients_invariantbasis[:,0] * (facq_x**2)
kappa_yy = 2 * arr_coefficients_invariantbasis[:,1] * (facq_x**2)
K_gauss_invariantbasis = (facq_x**4) * (arr_coefficients_invariantbasis[:,0]*arr_coefficients_invariantbasis[:,1] - arr_coefficients_invariantbasis[:,2]**2)


#####################
kappa_tt_profile = 2 * arr_coefficients_movingbasis[:,0] * (facq_x**2)
kappa_nn_profile = 2 * arr_coefficients_movingbasis[:,1] * (facq_x**2)
K_gauss_movingbasis = (facq_x**4) * (arr_coefficients_movingbasis[:,0]*arr_coefficients_movingbasis[:,1] - arr_coefficients_movingbasis[:,2]**2)

"""plt.figure()
plt.plot(kappa_nn_profile)
plt.plot(kappa_tt_profile)
plt.plot(K_gauss_movingbasis)"""


#######################

# vérif que les 2 methodes de mesure de courbure sont cohérentes
"""plt.figure()
plt.plot(kappa_xx,label='kappaxx')
plt.plot(kappa_yy,label='kappayy')
plt.plot(kappa_x_profile,label='kappax')
plt.plot(kappa_y_profile,label='kappay')
plt.legend()
plt.figure()
plt.plot(kappa_tt_profile, label='kappa_tt')
plt.plot(kappa_nn_profile, label='kappa_nn')
plt.plot(kappa_t_profile, label='kappa_t_profile')
plt.plot(kappa_n_profile, label='kappa_n_profile')
plt.legend()"""





x1 = dict_frac[keyfrac]['idcs_single_frac'][:,0]
y1 = dict_frac[keyfrac]['idcs_single_frac'][:,1]

x2 = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['x_ind']
y2 = dic_all_lines['frame_'+str(ind_tfrac_approx)]['dict_lines']['line_'+str(line_index)]['y_ind']

ind_closests = closest_points_between_2curves(x1=x1, y1=y1, x2=x2, y2=y2)

s_ind = np.zeros(len(x2))
for i in range(len(x2)-1):
    s_ind[i+1] = s_ind[i] + np.sqrt( (x2[i+1]-x2[i])**2 + (y2[i+1]-y2[i])**2 )

#xvals2plot_meters = np.arange(len(kappa_nn_profile)) * (1/facq_x)
xvals2plot_meters = s_ind * (1/facq_x)

%matplotlib inline


plt.figure(figsize=(10,7))
plt.plot(xvals2plot_meters, kappa_nn_profile, label='curvature normal to wave front')
plt.plot(xvals2plot_meters, kappa_tt_profile, label='curvature tangent to wave front')

plt.plot(xvals2plot_meters, np.sign(K_gauss_movingbasis) * np.sqrt(np.abs(K_gauss_movingbasis)))


#plt.plot(xvals2plot_meters, kappa_n_profile,label='kappa_n_profile method 1')
#plt.plot(xvals2plot_meters, kappa_t_profile,label='kappa_t_profile method 1')

plt.vlines(np.min(xvals2plot_meters[ind_closests]),np.min(kappa_nn_profile), np.max(kappa_nn_profile),'r',linestyle='--')
plt.vlines(np.max(xvals2plot_meters[ind_closests]),np.min(kappa_nn_profile), np.max(kappa_nn_profile),'r',linestyle='--')
plt.plot([],[],'r--',label='limites inf et sup de la fracture')
plt.legend()
plt.ylim(np.min(kappa_nn_profile), np.max(kappa_nn_profile))
plt.xlim(0, np.max(xvals2plot_meters[ind_closests] * 1.5))
plt.grid()
plt.xlabel('Position along the wave front [m]', fontsize=15)
plt.ylabel('Curvature [m^-1]', fontsize=15)
# %%
