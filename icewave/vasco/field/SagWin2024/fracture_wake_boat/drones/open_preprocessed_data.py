#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage

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

%matplotlib qt

# faire un tableau de profil dans la direction de propagation de l'onde
from profile import *

angle_deg = dominating_angle_wave_deg

profiles = all_profiles(image, angle_deg)



plot_profilelines_onimage(image, profiles)



theta = np.deg2rad(angle_deg)
d = np.array([np.sin(theta), np.cos(theta)])
plt.figure()

for p in profiles:
    src = np.array(p['src'])

    # position du début du profil dans le repère de la direction d
    s0 = np.dot(src, d)

    x = s0 + np.arange(len(p['profile']))

    plt.plot(x, p['profile'], 'k', alpha=0.1)

plt.xlabel("Position le long de la direction du profil")
plt.show()

