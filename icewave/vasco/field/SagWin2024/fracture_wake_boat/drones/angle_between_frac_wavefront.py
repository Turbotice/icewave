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

from measure_anglefrac import complete_routine_measure_angle_frac_lines, projection_coord_oninclinated_line
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
    """
    va chercher, pour un point (xind,yind) avec pour instant de fracture t_frac, le point sur le front d'onde détecté le plus proche
    
    Returns :
    closest_distance, id_closest_line, xshift, yshift, agmn_dist
    (avec agmn_dist l'indice, dans la ligne de front d'onde d'indice
    id_closest_line, du point trouvé le plus proche)
    """
    
    t_ind = np.where(dict_stereo_pivdata['t']>=tfrac)[0][0]
    dict_lines = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']

    min_distances = np.zeros(len(dict_lines))
    xshifts = np.zeros(len(dict_lines))
    yshifts = np.zeros(len(dict_lines))
    agmn_dists = np.zeros(len(dict_lines))
    for i in range(len(dict_lines)):
        xidcs_line = dict_lines['line_'+str(i)]['x_ind']
        yidcs_line = dict_lines['line_'+str(i)]['y_ind']
        
        distances = np.hypot(xidcs_line - xind, yidcs_line - yind)
        agmn_dist = np.argmin(distances)
        min_dist = distances[agmn_dist]

        min_distances[i] = min_dist
        xshifts[i] = (xidcs_line - xind)[agmn_dist]
        yshifts[i] = (yidcs_line - yind)[agmn_dist]
        agmn_dists[i] = agmn_dist
    
    id_closest_line = np.argmin(min_distances)
    closest_distance = min_distances[id_closest_line]
    xshift = xshifts[id_closest_line]
    yshift = yshifts[id_closest_line]
    ind_closest_onwavefront = int(agmn_dists[id_closest_line])


    # ajouter le calcule de "où est" le front par rapport à la fracture    

    return closest_distance, id_closest_line, xshift, yshift, ind_closest_onwavefront, t_ind

#%%
# test
array_closest_distance = np.zeros(len(array_tfrac_approx))
yshifts = np.zeros(len(array_tfrac_approx))
for i in range(len(array_tfrac_approx)):
    closest_distance, id_closest_line, xshift, yshift, _, _ = shortest_distance_to_any_wavefront(xind=xind_positions_fractures[i], yind=yind_positions_fractures[i], tfrac=array_tfrac_approx[i])
    array_closest_distance[i] = closest_distance
    yshifts[i] = yshift
%matplotlib qt
#plt.style.use('dark_background')
plt.figure()
plt.scatter(xind_positions_fractures, yind_positions_fractures,c=array_closest_distance*np.sign(yshifts), cmap='seismic', vmin=-20, vmax=20)
plt.colorbar()
plt.show()
# %%
# partie angles !!
dict_lines_frac, dict_results_anglesfrac = complete_routine_measure_angle_frac_lines(fractures_positions_data=fractures_positions_data, minsize_frac_px=6)

#%%

Xind_arr = dict_results_anglesfrac['Xind_arr']
Yind_arr = dict_results_anglesfrac['Yind_arr']
Alpha_linfit_arr = dict_results_anglesfrac['Alpha_linfit_arr']

Alpha_wavefront_closest_arr = np.zeros(len(Xind_arr))

for i in range(len(Xind_arr)):

    ind_same = np.where((Xind_arr[i]==xind_positions_fractures)&(Yind_arr[i]==yind_positions_fractures))[0][0]

    xind = xind_positions_fractures[ind_same]
    yind = yind_positions_fractures[ind_same]
    tfrac = array_tfrac_approx[ind_same]

    closest_distance, id_closest_line, xshift, yshift, ind_closest_onwavefront, t_ind = shortest_distance_to_any_wavefront(xind=xind, yind=yind, tfrac=tfrac, dict_all_wavefront_lines=dict_all_wavefront_lines)

    xind_wavefront = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']['line_'+str(id_closest_line)]['x_ind']
    yind_wavefront = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']['line_'+str(id_closest_line)]['y_ind']
    
    u_t = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']['line_'+str(id_closest_line)]['polynomial_fit']['u_t']
    v_t = dict_all_wavefront_lines['frame_'+str(t_ind)]['dict_lines']['line_'+str(id_closest_line)]['polynomial_fit']['v_t']

    u_t_closest = u_t[ind_closest_onwavefront]
    v_t_closest = v_t[ind_closest_onwavefront]
    alpha_wavefront_closest = np.arctan2(v_t_closest, u_t_closest)
    Alpha_wavefront_closest_arr[i] = alpha_wavefront_closest


# %%
%matplotlib inline

maskplot = (Yind_arr>=20)&(Yind_arr<70) & ((Alpha_linfit_arr - Alpha_wavefront_closest_arr)>0)

plt.figure()
plt.scatter(Xind_arr, Yind_arr, c=(Alpha_linfit_arr) * 180/np.pi, vmin=0, vmax=90, cmap='rainbow')
plt.colorbar()
plt.title('Angle of fracture lines [deg]')
plt.show()

plt.figure()
plt.plot(Yind_arr, (Alpha_linfit_arr - Alpha_wavefront_closest_arr) * 180/np.pi, '.')
plt.xlabel('y [px]', fontsize=13)
plt.ylabel('($alpha_{frac}$ - $alpha_{wavefront}$) [deg]', fontsize=13)
plt.show()

#%matplotlib qt

alpha_frac_deg = 40 # on prend un angle typique
alpha_frac_rad = np.radians(alpha_frac_deg)
slope = np.tan(alpha_frac_rad)
intercept = 0 # la valeur de intercept n'importe pas

Xprim_along_line = projection_coord_oninclinated_line(x_arr=Xind_arr, y_arr=Yind_arr, slope_line=slope, intercept_line=intercept)

plt.figure()
plt.plot(Xprim_along_line[maskplot], (Alpha_linfit_arr - Alpha_wavefront_closest_arr)[maskplot] * 180/np.pi, '.')
#plt.plot(Xprim_along_line[np.logical_not(maskplot)], (Alpha_linfit_arr - Alpha_wavefront_closest_arr)[np.logical_not(maskplot)] * 180/np.pi, '.', alpha=0.2)
plt.xlabel('x\' along average direction of wave front [px]', fontsize=13)
plt.ylabel('($alpha_{frac}$ - $alpha_{wavefront}$) [deg]', fontsize=13)
#plt.ylim(0,90)
#plt.xlim(0)

from scipy.ndimage import gaussian_filter1d
indsort = np.argsort(Xprim_along_line[maskplot])
xsm = Xprim_along_line[maskplot][indsort]
ysm = gaussian_filter1d((Alpha_linfit_arr - Alpha_wavefront_closest_arr)[maskplot][indsort] * 180/np.pi,sigma=10)
plt.plot(xsm, ysm)
plt.ylim(0, 75)
plt.xlim(0, np.max(xsm)*1.05)
plt.show()

# %%
