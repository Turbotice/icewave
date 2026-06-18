#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.io import loadmat
import pickle
from scipy.optimize import curve_fit

import sys

vasco_git_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/vasco'
sys.path.append(vasco_git_path)

from tools.front_tracking import find_peak_near_x0

from manips_Bs_As.tracking_force_displacement.functions.load_force_vs_frame import *

#%% inputs

disk = 'D:'

date = '1203'
serie = 2
acq = 2
camnum = 1

img_dir = f"{disk}/manips_BsAs/Basler_images/{date}/serie{serie}/acq{acq}_expo2000_facq50/cam{camnum}/"
prefix_imgname = 'Basler_a2A1920-160umBAS__40575305__20251203_154603379_'


i0 = 1800
ifrac = 1861

# plate properties and spatial scales
dcm = 4
dpx = 165

L = 16e-2 # length plate [m]
w = 4e-2 # width plate [m]
h_avg = 2.15e-3 # plate thickness [m]
h_std = 0.14e-3 # plate thickness error [m]

# force dir path and scale on image to measure force
forcedirpath = f"{disk}/manips_BsAs/Basler_images/{date}/serie{serie}/acq{acq}_expo2000_facq50"
dg_sur_dpx = 300/405


# %% charger les images
# open 100 images, named from Basler_acA2040-90um__23029848__20260127_182642777_0420
# to Basler_acA2040-90um__23029848__20260127_182642777_0429, and store their pixel values in numpy arrays
# (use a 3d numpy array to store the pixel values of all 10 images, with shape (10, height, width))



images = []
arr_frames = np.arange(i0, ifrac+10, 1)

for i in range(len(arr_frames)):
    frame = arr_frames[i]
    filename = f"{img_dir}{prefix_imgname}{str(frame).zfill(4)}.tiff"
    image = plt.imread(filename)
    images.append(image)
images_array = np.array(images)


# %% find edge ice plate to have profile

x_plate_lim_approx = 420

yinf = 0
ysup = images_array.shape[1]

images_array_diffx = np.gradient(images_array, axis=2)

Peaks_found = np.zeros((images_array.shape[0], images_array.shape[1]))

for i in range(Peaks_found.shape[0]):
    for j in range(Peaks_found.shape[1]):
        pk = find_peak_near_x0(images_array_diffx[i,j,:], x_plate_lim_approx, window_size=20)
        Peaks_found[i,j] = pk

#%% load force vs time data

dict_force = load_force_vs_frames_1balance_from_dirpath(forcedirpath, dg_sur_dpx= dg_sur_dpx,plot=True)

force_N = 9.81 * dict_force['force_grams'] * 1e-3
arr_frames_forces = dict_force['frames']

common, idx_a, idx_b = np.intersect1d(arr_frames, arr_frames_forces, return_indices=True)
force_N_common = force_N[idx_b]


#%%
def beam_profile_theory(x, F=1, E=1e9, L=8e-2, w=4e-2, h=4e-3):
    x = np.asarray(x)
    I = (1/12) * w * h**3
    x_sym = L - x  # mirror of x with respect to x=L/2
    y1 = (F * x) / (48 * E * I) * (3 * L**2 - 4 * x**2)
    y2 = (F * x_sym) / (48 * E * I) * (3 * L**2 - 4 * x_sym**2)
    return np.where(x <= L/2, y1, y2)


yvals_meters = np.arange(images_array.shape[1]) * (dcm/dpx) * 1e-2
initial_profile = Peaks_found[0,:]
plt.figure()
for i in range(0, images_array.shape[0]-10, 10):
    mask1 = (Peaks_found[i,:]>=x_plate_lim_approx-10)&(Peaks_found[i,:]<x_plate_lim_approx+10)
    mask2 = (Peaks_found[0,:]>=x_plate_lim_approx-10)&(Peaks_found[0,:]<x_plate_lim_approx+10)
    mask = mask1 & mask2


    plt.plot(yvals_meters[mask], 1e-2*(dcm/dpx) * (Peaks_found[i,:]-initial_profile)[mask])

plt.plot(yvals_meters, beam_profile_theory(yvals_meters+0.012, F=0.85, E=3e9, L=L, h=h_avg), 'k')
# %%

ind_times_common = np.arange(images_array.shape[0])[idx_a]

Peaks_found_common = Peaks_found[idx_a,:]
initial_profile = Peaks_found_common[0,:]

colors = plt.cm.viridis(np.linspace(0, 1, len(range(0, len(ind_times_common), 10))))
idx=0
for i in range(0, len(ind_times_common), 10):
    force = force_N_common[i]

    mask1 = (Peaks_found_common[i,:]>=x_plate_lim_approx-10)&(Peaks_found_common[i,:]<x_plate_lim_approx+10)
    mask2 = (Peaks_found_common[0,:]>=x_plate_lim_approx-10)&(Peaks_found_common[0,:]<x_plate_lim_approx+10)
    mask = mask1 & mask2

    color = colors[idx]
    idx+=1
    plt.plot(yvals_meters[mask], 1e-2*(dcm/dpx) * (Peaks_found_common[i,:]-initial_profile)[mask], color=color, label='F='+str(np.round(force,2))+' N')
    plt.plot(yvals_meters, beam_profile_theory(yvals_meters+0.012, F=force, E=4e9, L=L, h=h_avg), '--', color=color, alpha=0.5)
plt.legend()

plt.ylim(-1e-4,1e-3)