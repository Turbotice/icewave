#%%
import numpy as np
import matplotlib.pyplot as plt

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage, get_n_points_anyfigure

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken



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

ux = np.cumsum(vx, axis=2)*dt
uy = np.cumsum(vy, axis=2)*dt
uz = np.cumsum(vz, axis=2)*dt


#%%
%matplotlib inline

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)

fig, ax = plt.subplots()
ax.imshow(uz[:,:,1000])
ax.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')
plt.show()

#%%
dict_longfrac = {}

nb_longfrac = 30

fig, ax = plt.subplots()
ax.imshow(uz[:,:,1000])
ax.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')


%matplotlib qt

for i in range(nb_longfrac):
    n_points = int(input('number of pixels in the fracture'))

    coords = get_n_points_onimage(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan), n_points=n_points)

    dict_longfrac['longfrac_'+str(i)] = {}
    dict_longfrac['longfrac_'+str(i)]['coords'] = coords

#%%
Lines = find_lines(fractures_positions_data['binary'], tolerance_px=2.5)

import matplotlib.colors as mcolors
mcolors.CSS4_COLORS
colors = [list(mcolors.CSS4_COLORS.keys())[i] for i in range(len(list(mcolors.CSS4_COLORS)))]

for i in range(len(Lines)):
    for j in range(len(Lines[i])):
        plt.plot(Lines[i][j][0], Lines[i][j][1], 'o', color=colors[i])

# %%
