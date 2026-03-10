#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import csv
import re
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator
from PIL import Image
from scipy.signal import correlate2d
from scipy.ndimage import gaussian_filter1d
#%%

%matplotlib qt

# %% load data
W = 32
#Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 300
N = 1200
refimg = 300

frame_frac = 1146
frame_force_application = 300 #début de l'application de la force

dcm = 5
dpx = 221 # à changer en fonction de la manip
L = 8e-2 # distance entre 2 points d'appui

date = '20250514'
time = '143000'

system_loc = 'windows_server'

disk = 'D:/Grenoble/'
general_folder = f'{disk}/Gre25/Data/PIV_results/3points_bending/{date}_{time}/'

"""
if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/PIV_results/3points_bending/{date}_{time}/'
"""

path2data = general_folder

matfile = f'{path2data}PIV_processed_displacement_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

#matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'


from scipy.io import loadmat

mat_dict = loadmat(matfile)


u_original = mat_dict['u_original'][:,0]
v_original = mat_dict['v_original'][:,0]

def reshape_array(arr):
    array_new = np.zeros((len(arr),arr[0].shape[0],arr[0].shape[1]))
    for i in range(len(arr)):
        array_new[i,:,:] = arr[i]
    return array_new

u = reshape_array(u_original)
v = reshape_array(v_original)


xpix = mat_dict['x'][0][0][0]
ypix = mat_dict['y'][0][0][:,0]



# %%

yind = 2
xind = 5

 # average neighbours and smooth vs time :
avg_smooth = False
if avg_smooth:
    v_field = gaussian_filter1d(np.mean(v[:,yind-2:yind+2,xind-2:xind+2],axis=(1,2)),5)
else:
    v_field = v[:,yind,xind]

plt.figure()
plt.plot(np.arange(v.shape[0])+i0,v_field)
plt.xlim(i0,frame_frac)
plt.ylim(-np.max(v_field)/4,np.max(v_field)/4)
plt.show()




plt.figure()
plt.title('vertical displacement (m) vs time (frames)')
plt.plot(np.arange(v.shape[0])+i0,1e-2 * v_field * dcm/dpx)
plt.xlim(i0,frame_frac)
plt.ylim(-1e-2 * (dcm/dpx)*np.max(v_field[:frame_frac-i0])/4,1e-2 * (dcm/dpx)*np.max(v_field[:frame_frac-i0]))
plt.show()



h_measurements = np.array([5.25, 5.38, 5.0, 4.8])
h_avg_mm = np.mean(h_measurements)
h_std_mm = np.std(h_measurements)

h_avg = 1e-3 * h_avg_mm
h_std = 1e-3 * h_std_mm
print('h='+str(h_avg)+'+-'+str(h_std))


def compute_kappa(delta,L):
    return 8*delta/(L**2)

def compute_epsilon(kappa,h):
    return kappa*h/2
# calcul approx pour odg :
delta_c = 6.9e-4 # fleche max

kappa_c = compute_kappa(delta_c,8e-2) # pour l'insant on va chercher la flèche à l'oeil
epsilon_c = compute_epsilon(kappa_c,h_avg)

epsilon_c_2 = (6*h_avg/((L)**2)) * delta_c
print(epsilon_c_2)

print('courbure critique : ',kappa_c,'m^-1')
print('déformation critique : ',epsilon_c)




# %%
