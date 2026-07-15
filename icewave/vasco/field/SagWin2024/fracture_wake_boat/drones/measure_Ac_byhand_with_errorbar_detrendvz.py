#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken, clic2extract_amplitude_withref_singlepixel, oneclick2extract_amplitude_singlepixel, extract_amplitude_singlepixel_automatic

#%%
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

# %%
dict_stereo_pivdata.keys()

vx = dict_stereo_pivdata['u'][0,:,:,:]
vy = dict_stereo_pivdata['u'][1,:,:,:]
vz = dict_stereo_pivdata['u'][2,:,:,:]

dt = dict_stereo_pivdata['t'][1] - dict_stereo_pivdata['t'][0]

def detrend_along_time_axis(v, window_size=120): # à peu près nombre de frames par periode d'onde):
    window = np.ones(window_size) / window_size
    v_detrended = np.zeros_like(v)

    # Supposons vz.shape = (N_x, N_y, N_t)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            v_detrended[i, j, :] = v[i, j, :] - np.convolve(v[i, j, :], window, mode='same')
    return v_detrended

vx_detrended = detrend_along_time_axis(vx)
vy_detrended = detrend_along_time_axis(vy)
vz_detrended = detrend_along_time_axis(vz)


ux = np.cumsum(vx_detrended, axis=2)*dt
uy = np.cumsum(vy_detrended, axis=2)*dt
uz = np.cumsum(vz_detrended, axis=2)*dt


#%% PLOTS
%matplotlib inline
for i in range(0,1000,100):
    plt.figure()
    plt.imshow(vz_detrended[:,:,1000+i])



#%%
%matplotlib inline

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac_detrendvz.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)

image = uz[:,:,1000]
plt.imshow(image)
plt.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')

if 'dict_frac' in globals():
    for k in dict_frac:
        if 'single_frac' in k:
            plt.plot(dict_frac[k]['idcs_single_frac'][0][0], dict_frac[k]['idcs_single_frac'][0][1], '^r', markersize=4)
plt.show()
#%%
# soit on refait en cliquant pour choisir un nouveau crack (code commenté) :
"""
n_points = 10
%matplotlib qt

matrix_temp_evol_uz, times_frac_sec_approx, idcs_single_frac = click_on_fracture_path_plot_time_evol(n_points, uz, dict_stereo_pivdata,fractures_positions_data)
# LE PREMIER CLIC DOIT ETRE LE PLUS PROCHE POSSIBLE DE LA LIMITE AVEC GLACE NON CASSEE !!
"""
# soit on fait en allant chercher les informations des cracks déjà sélectionnés auparavant
path_dict_frac_ref = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
with open(path_dict_frac_ref, 'rb') as file:
    dict_frac_ref = pickle.load(file)

# Now choose the 'keyfrac' :
keyfrac = 'dict_single_frac_yind23_xind28'
times_frac_sec_approx = dict_frac_ref[keyfrac]['times_frac_sec_approx']
idcs_single_frac = dict_frac_ref[keyfrac]['idcs_single_frac']

matrix_temp_evol_uz = np.zeros((idcs_single_frac.shape[0], len(dict_stereo_pivdata['t'])))
for i in range(idcs_single_frac.shape[0]):
    matrix_temp_evol_uz[i,:] = uz[idcs_single_frac[i,1], idcs_single_frac[i,0],:]

# %%
#array_amplitudes_frac, wave_period_sec, array_t1_sec, array_t2_sec = click2extract_amplitude(matrix_temp_evol_uz, times_frac_sec_approx, dict_stereo_pivdata)

array_amplitudes_frac = []
#array_amplitudes_err_frac = []
dict_all_results_atfracpixels = {}

for i in range(len(matrix_temp_evol_uz)):
    dict_results_atfracpixels = extract_amplitude_singlepixel_automatic(temp_evol_uz=matrix_temp_evol_uz[i], t_frac_sec_approx=times_frac_sec_approx[i], dict_stereo_pivdata=dict_stereo_pivdata, arr_t_frac_sec=times_frac_sec_approx)
    array_amplitudes_frac.append(dict_results_atfracpixels['amplitudemax'])
    #array_amplitudes_err_frac.append(dict_results_atfracpixels['amplitude_err'])
    dict_all_results_atfracpixels[str(i)] = dict_results_atfracpixels

#%%

%matplotlib inline

array_D = np.array([2,4,6,8,10,12,14]) # on définit une longueur (en pixel) où on veut se placer pour regarder 
#array_D = np.array([2,4,6,8])
# la ref en glace non cassée

# on met dans des matrices de la même dimensions que le cas précédent (donc si une seule valeur on la transforme en tableau)
matrix_temp_evol_uz_ref_noncassee = np.zeros((len(array_D), uz.shape[2]))


for i in range(len(array_D)):
    D=array_D[i]
    uz_ref_noncassee, x_ind_ref_noncasse_int, y_ind_ref_noncasse_int, time_frac_approx, (x_ind_ref_noncasse,y_ind_ref_noncasse,xdata2plot,ydata2plot,slope,intercept, deltaD) = plot_elevation_refnotbroken_and_broken(D, uz, idcs_single_frac, times_frac_sec_approx,
                                            dict_stereo_pivdata, matrix_temp_evol_uz,
                                            fractures_positions_data)
    
    matrix_temp_evol_uz_ref_noncassee[i,:] = uz_ref_noncassee

#times_frac_sec_approx_ref_noncassee = np.ones(len(array_D)) * time_frac_approx

# maintenant mesurer également l'amplitude pour la ref "glace non cassée"
%matplotlib qt

array_amplitudes_ref_noncassee = np.zeros(len(array_D))
#array_amplitudes_err_ref_noncassee = np.zeros(len(array_D))
dict_all_results_ref_noncassee = {}

for i in range(len(matrix_temp_evol_uz_ref_noncassee)):
    dict_results_ref_noncassee = oneclick2extract_amplitude_singlepixel(temp_evol_uz=matrix_temp_evol_uz_ref_noncassee[i], t_frac_sec_approx=np.mean(times_frac_sec_approx), dict_stereo_pivdata=dict_stereo_pivdata, arr_t_frac_sec=times_frac_sec_approx)
    #dict_results_ref_noncassee = extract_amplitude_singlepixel_automatic(temp_evol_uz=matrix_temp_evol_uz_ref_noncassee[i], t_frac_sec_approx=np.inf, dict_stereo_pivdata=dict_stereo_pivdata, arr_t_frac_sec=times_frac_sec_approx)
    array_amplitudes_ref_noncassee[i] = dict_results_ref_noncassee['amplitude']
    #array_amplitudes_err_ref_noncassee[i] = dict_results_ref_noncassee['amplitude_err']
    dict_all_results_ref_noncassee[str(i)] = dict_results_ref_noncassee
    

#%%

%matplotlib inline
# on veut calculer l' "abscisse" sur la droite fittée de fracture de tous ces points
alpha = np.arctan(slope)

coord_onfracline_noncassee = - array_D
x_ind_onfrac = idcs_single_frac[:,0]
y_ind_onfrac = idcs_single_frac[:,1]

def func_fit(x,a=slope,b=intercept):
    return a*x+b

coord_onfracline_cassee = (func_fit(x_ind_onfrac)-func_fit(x_ind_onfrac[0]))/np.sin(alpha) + ((func_fit(x_ind_onfrac[0]) - y_ind_onfrac[0]) - (func_fit(x_ind_onfrac) - y_ind_onfrac))*np.sin(alpha)

%matplotlib qt

array_amplitudes_err_frac = np.zeros(len(array_amplitudes_frac))
array_amplitudes_err_ref_noncassee = np.zeros(len(array_amplitudes_ref_noncassee))

plt.figure()
plt.errorbar(coord_onfracline_cassee, array_amplitudes_frac, array_amplitudes_err_frac, linestyle='', marker='o', color='k', ecolor='gray')
plt.errorbar(coord_onfracline_cassee, array_amplitudes_frac, array_amplitudes_err_frac, linestyle='', marker='o', color='k', ecolor='gray')
plt.errorbar(coord_onfracline_noncassee, array_amplitudes_ref_noncassee, array_amplitudes_err_ref_noncassee, linestyle='', marker='o', color='k', ecolor='gray')
plt.ylabel('Amplitude [m]')
plt.xlabel('Relative distance [px] to the first "fractured" pixel detected')
plt.show()

#%%

# save dico example

dict_example_timeevol = {}
dict_example_timeevol['matrix_temp_evol_uz'] = matrix_temp_evol_uz
dict_example_timeevol['tvals_sec'] = dict_stereo_pivdata['t']
dict_example_timeevol['times_frac_sec_approx'] = times_frac_sec_approx
dict_example_timeevol['keyfrac'] = keyfrac
dict_example_timeevol['idcs_single_frac'] = idcs_single_frac

dict_example_timeevol['matrix_temp_evol_uz_ref_noncassee'] = matrix_temp_evol_uz_ref_noncassee

dict_example_timeevol['coord_onfracline_cassee'] = coord_onfracline_cassee
dict_example_timeevol['coord_onfracline_noncassee'] = coord_onfracline_noncassee
dict_example_timeevol['array_amplitudes_frac'] = array_amplitudes_frac


path2dict_example_timeevol = "C:/Users/Vasco Zanchi/Desktop/presentations/reunions_hebdo/figures_data_article/dict_timeevol_example.pkl"
pickle.dump(dict_example_timeevol, open(path2dict_example_timeevol, "wb"))

#%% save data from one fracture line
frac_id = f'yind{idcs_single_frac[0,1]}_xind{idcs_single_frac[0,0]}'

dict_single_frac = {}
dict_single_frac['matrix_temp_evol_uz'] = matrix_temp_evol_uz
dict_single_frac['times_frac_sec_approx'] = times_frac_sec_approx
dict_single_frac['idcs_single_frac'] = idcs_single_frac
dict_single_frac['array_amplitudes_frac'] = array_amplitudes_frac
#dict_single_frac['wave_period_sec'] = array_wave_period_sec
dict_single_frac['frac_id'] = frac_id
dict_single_frac['coord_onfracline_cassee'] = coord_onfracline_cassee


dict_single_frac['array_amplitudes_ref_noncassee'] = array_amplitudes_ref_noncassee
#dict_single_frac['wave_period_sec_ref_noncassee'] = array_wave_period_sec_ref_noncassee
#dict_single_frac['times_frac_sec_approx_ref_noncassee'] = times_frac_sec_approx_ref_noncassee
dict_single_frac['matrix_temp_evol_uz_noncassee'] = matrix_temp_evol_uz_ref_noncassee
dict_single_frac['slope'] = slope
dict_single_frac['intercept'] = intercept
dict_single_frac['alpha'] = alpha
dict_single_frac['coord_onfracline_noncassee'] = coord_onfracline_noncassee
dict_single_frac['array_D'] = array_D



path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac_detrendvz.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)

else:
    dict_frac = {}

dict_frac['dict_single_frac_'+dict_single_frac['frac_id']] = dict_single_frac

save_input = input('save ? (y/n)')
if save_input=='y':
    pickle.dump(dict_frac, open(path_dict2save, "wb"))
else:
    pass



#%%

%matplotlib inline
plt.figure()

colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
ct = 0
for k in dict_frac:

    coord_onfracline_cassee = dict_frac[k]['coord_onfracline_cassee']
    coord_onfracline_noncassee = dict_frac[k]['coord_onfracline_noncassee']
    array_amplitudes_frac = dict_frac[k]['array_amplitudes_frac']
    array_amplitudes_ref_noncassee = dict_frac[k]['array_amplitudes_ref_noncassee']

    #wave_period_sec = dict_frac[k]['wave_period_sec']

    plt.plot(coord_onfracline_noncassee, array_amplitudes_ref_noncassee,'s',color=colors[ct])
    plt.plot(coord_onfracline_cassee, array_amplitudes_frac,'^',color=colors[ct])

    plt.ylabel('Amplitude [m]')
    plt.xlabel('Relative distance [px] to the first "fractured" pixel detected')

    ct+=1


# %%
