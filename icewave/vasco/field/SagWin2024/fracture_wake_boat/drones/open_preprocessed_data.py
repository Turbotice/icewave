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

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken

#%%
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

# %%
dict_stereo_pivdata.keys()

vx = dict_stereo_pivdata['u'][0,:,:,:]
vy = dict_stereo_pivdata['u'][1,:,:,:]
vz = dict_stereo_pivdata['u'][2,:,:,:]

dt = dict_stereo_pivdata['t'][1] - dict_stereo_pivdata['t'][0]

ux = np.cumsum(vx, axis=2)*dt
uy = np.cumsum(vy, axis=2)*dt
uz = np.cumsum(vz, axis=2)*dt


#%% PLOTS
%matplotlib inline
for i in range(0,1000,100):
    plt.figure()
    plt.imshow(vz[:,:,1000+i])

#%% FFT2 in space
%matplotlib qt
# mesure de l'angle de l'onde pour une image
image = vz[:,:,1000]

image_croped = image[:,:60] # car sur la gauche l'angle est un peu different

padding_factor = 5

absfft2shifted = np.fft.fftshift(np.abs(np.fft.fft2(image_croped, s=(vz.shape[0]*padding_factor,vz.shape[1]*padding_factor))))

# detection of maxima
maxval = np.max(absfft2shifted)
npwh = np.where(absfft2shifted==maxval)
kyindmax = npwh[0][0]
kxindmax = npwh[1][0]
print(npwh)

kyzeroind = int(vz.shape[0]*padding_factor/2)
kxzeroind = int(vz.shape[1]*padding_factor/2)

diff_kxind = kxindmax - kxzeroind
diff_kyind = kyindmax - kyzeroind 

dominating_angle_wave_deg = 180/np.pi * np.arctan(diff_kyind/diff_kxind)
print('angle wave =',dominating_angle_wave_deg,'deg')

plt.figure()
plt.imshow(absfft2shifted)
plt.plot(kxindmax, kyindmax,'r+')
plt.show()




#%%

#%matplotlib qt

# faire un tableau de profil dans la direction de propagation de l'onde
from profiles import *

angle_deg = dominating_angle_wave_deg

profiles = all_profiles(image, angle_deg)



plot_profilelines_onimage(image, profiles)


#%%
%matplotlib inline

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)

plt.imshow(image)
plt.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')

if 'dict_frac' in globals():
    for k in dict_frac:
        plt.plot(dict_frac[k]['idcs_single_frac'][0][0], dict_frac[k]['idcs_single_frac'][0][1], '^r', markersize=4)
plt.show()
#%%

n_points = 11
%matplotlib qt

matrix_temp_evol_uz, times_frac_sec_approx, idcs_single_frac = click_on_fracture_path_plot_time_evol(n_points, uz, dict_stereo_pivdata,fractures_positions_data)
# LE PREMIER CLIC DOIT ETRE LE PLUS PROCHE POSSIBLE DE LA LIMITE AVEC GLACE NON CASSEE !!

#frac_id = f'xind{}'

# %%
array_amplitudes_frac, wave_period_sec, array_t1_sec, array_t2_sec = click2extract_amplitude(matrix_temp_evol_uz, times_frac_sec_approx, dict_stereo_pivdata)

#%%

%matplotlib inline

array_D = np.array([2,4,6,8,10,12,14]) # on définit une longueur (en pixel) où on veut se placer pour regarder 
#la ref en glace non cassée

# on met dans des matrices de la même dimensions que le cas précédent (donc si une seule valeur on la transforme en tableau)
matrix_temp_evol_uz_ref_noncassee = np.zeros((len(array_D), uz.shape[2]))


for i in range(len(array_D)):
    D=array_D[i]
    uz_ref_noncassee, x_ind_ref_noncasse_int, y_ind_ref_noncasse_int, time_frac_approx, (x_ind_ref_noncasse,y_ind_ref_noncasse,xdata2plot,ydata2plot,slope,intercept, deltaD) = plot_elevation_refnotbroken_and_broken(D, uz, idcs_single_frac, times_frac_sec_approx,
                                            dict_stereo_pivdata, matrix_temp_evol_uz,
                                            fractures_positions_data)
    
    matrix_temp_evol_uz_ref_noncassee[i,:] = uz_ref_noncassee

times_frac_sec_approx_ref_noncassee = np.ones(len(array_D)) * time_frac_approx

# maintenant mesurer également l'amplitude pour la ref "glace non cassée"
%matplotlib qt
array_amplitudes_ref_noncassee, wave_period_sec_ref_noncassee, array_t1_sec_noncassee, array_t2_sec_noncassee = click2extract_amplitude(matrix_temp_evol_uz_ref_noncassee, times_frac_sec_approx_ref_noncassee, dict_stereo_pivdata)

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
plt.figure()
plt.plot(coord_onfracline_noncassee, array_amplitudes_ref_noncassee,'o')
plt.plot(coord_onfracline_cassee, array_amplitudes_frac,'o')

plt.ylabel('Amplitude [m]')
plt.xlabel('Relative distance [px] to the first "fractured" pixel detected')
plt.show()

#%% save data from one fracture line
frac_id = f'yind{idcs_single_frac[0,1]}_xind{idcs_single_frac[0,0]}'

dict_single_frac = {}
dict_single_frac['matrix_temp_evol_uz'] = matrix_temp_evol_uz
dict_single_frac['times_frac_sec_approx'] = times_frac_sec_approx
dict_single_frac['idcs_single_frac'] = idcs_single_frac
dict_single_frac['array_amplitudes_frac'] = array_amplitudes_frac
dict_single_frac['wave_period_sec'] = array_wave_period_sec
dict_single_frac['array_t1_sec'] = array_t1_sec
dict_single_frac['array_t2_sec'] = array_t2_sec
dict_single_frac['frac_id'] = frac_id
dict_single_frac['coord_onfracline_cassee'] = coord_onfracline_cassee


dict_single_frac['array_amplitudes_ref_noncassee'] = array_amplitudes_ref_noncassee
dict_single_frac['wave_period_sec_ref_noncassee'] = array_wave_period_sec_ref_noncassee
dict_single_frac['array_t1_sec_noncassee'] = array_t1_sec_noncassee
dict_single_frac['array_t2_sec_noncassee'] = array_t2_sec_noncassee
dict_single_frac['times_frac_sec_approx_ref_noncassee'] = times_frac_sec_approx_ref_noncassee
dict_single_frac['matrix_temp_evol_uz_noncassee'] = matrix_temp_evol_uz_ref_noncassee
dict_single_frac['slope'] = slope
dict_single_frac['intercept'] = intercept
dict_single_frac['alpha'] = alpha
dict_single_frac['coord_onfracline_noncassee'] = coord_onfracline_noncassee
dict_single_frac['array_D'] = array_D



path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
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

    wave_period_sec = dict_frac[k]['wave_period_sec']

    plt.plot(coord_onfracline_noncassee, array_amplitudes_ref_noncassee,'s',color=colors[ct])
    plt.plot(coord_onfracline_cassee, array_amplitudes_frac,'^',color=colors[ct])

    plt.ylabel('Amplitude [m]')
    plt.xlabel('Relative distance [px] to the first "fractured" pixel detected')

    ct+=1




#%%


"""
frac_id = dict_frac['dict_single_frac_yind23_xind28']['frac_id']
array_amplitudes_frac = dict_frac['dict_single_frac_yind23_xind28']['array_amplitudes_frac']
idcs_single_frac = dict_frac['dict_single_frac_yind23_xind28']['idcs_single_frac']
matrix_temp_evol_uz = dict_frac['dict_single_frac_yind23_xind28']['matrix_temp_evol_uz']
times_frac_sec_approx = dict_frac['dict_single_frac_yind23_xind28']['times_frac_sec_approx']
wave_period_sec = dict_frac['dict_single_frac_yind23_xind28']['wave_period_sec']
coord_onfracline_cassee = dict_frac['dict_single_frac_yind23_xind28']['coord_onfracline_cassee']

array_amplitudes_ref_noncassee = dict_frac['dict_single_frac_yind23_xind28']['array_amplitudes_ref_noncassee']
wave_period_sec_ref_noncassee = dict_frac['dict_single_frac_yind23_xind28']['wave_period_sec_ref_noncassee']
matrix_temp_evol_uz_ref_noncassee = dict_frac['dict_single_frac_yind23_xind28']['matrix_temp_evol_uz_noncassee']
times_frac_sec_approx_ref_noncassee = dict_frac['dict_single_frac_yind23_xind28']['times_frac_sec_approx_ref_noncassee']
slope = dict_frac['dict_single_frac_yind23_xind28']['slope']
intercept = dict_frac['dict_single_frac_yind23_xind28']['intercept']
alpha = dict_frac['dict_single_frac_yind23_xind28']['alpha']
coord_onfracline_noncassee = dict_frac['dict_single_frac_yind23_xind28']['coord_onfracline_noncassee']
array_D = dict_frac['dict_single_frac_yind23_xind28']['array_D']"""
# %%
