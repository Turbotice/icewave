#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 10:00:10 2025

@author: moreaul
"""

import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import gpxpy
import pandas as pd
from scipy.fft import fft2, fftshift, fftfreq,fft, ifft
from matplotlib.path import Path
from datetime import datetime
from scipy.signal import butter, sosfiltfilt
from scipy.signal import resample
from scipy.signal import decimate
from geopy.distance import distance

get_ipython().run_line_magic('matplotlib', 'qt')
plt.close('all')

date = '0211'
year = '2025'
path2data = f'/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/{year}_BICWIN/{date}/DAS'


xe_map = {
    '0210': 0.5,
    '0211': 2,
    '0212': 1,
}

# choix automatique
xe = xe_map.get(date, None)  # None si la date n’est pas dans la liste

if xe is None:
    raise ValueError(f"Aucune valeur définie pour la date {date}")




# Get sorted list of .h5 files in path2data
h5_files = sorted([f for f in os.listdir(path2data) if f.endswith('.h5')])
print(f"Found {len(h5_files)} .h5 files in {path2data}.")

# Map acquisition numbers to file names by index
file_dict = {i: fname for i, fname in enumerate(h5_files)}

# Set acquisition numbers to plot
acquisition_numbers_to_plot = [2]  # based on file order

# Loop through selected acquisition numbers and inspect each file
for acq_num in acquisition_numbers_to_plot:
    if acq_num in file_dict:
        filepath = os.path.join(path2data, file_dict[acq_num])
        print(f"\nReading acquisition {acq_num}: {file_dict[acq_num]}")
        with h5py.File(filepath, 'r') as h5file:
            print("Available datasets:")
            def print_structure(name, obj):
                if isinstance(obj, h5py.Dataset):
                    print(f" - {name}")
            h5file.visititems(print_structure)

            # Load time and strain datasets
            try:
                time = h5file['fa1-20050033/Source1/time'][:]
                strain = h5file['fa1-20050033/Source1/Zone1/Strain Rate [nStrain|s]'][:]

                print("Time shape:", time.shape)
                print("Strain shape:", strain.shape)

                # # Plot strain vs time for first 5 channels
                # plt.figure(figsize=(10, 6))
                # for ch in range(min(5, strain.shape[0])):
                #     plt.plot(strain[ch, :], label=f'Ch {ch}')
                # plt.xlabel("Time (s)")
                # plt.ylabel("Strain Rate [nStrain/s]")
                # plt.title(f"Acquisition {acq_num}: {file_dict[acq_num]}")
                # plt.legend()
                # plt.grid(True)
                # plt.tight_layout()
                # plt.show()

            except KeyError as e:
                print(f"Dataset missing: {e}")
    else:
        print(f"\nAcquisition number {acq_num} not found among .h5 files.")



Nseg,fe,Nx = strain.shape

data = np.zeros([Nx,Nseg*fe])
for ix in range(Nx):
    tmp = []
    for seg in range(Nseg):
        tmp.append(strain[seg, :, ix])
    tmp = np.concatenate(tmp)  # Now tmp is a 1D vector
    data[ix,:] = tmp
    
low, high = 0.05, 0.5
# design stable second-order-sections bandpass
sos = butter(4, [low, high], btype='band', fs=fe, output='sos')
# zero-phase filtering along time axis (axis=1)
data = sosfiltfilt(sos, data, axis=1)
    

with h5py.File(filepath, 'r') as h5file:
    time_data = h5file['fa1-20050033/Source1/time'][:]
    print("Time dataset shape:", time_data.shape)
    # print("First few values:", time_data[:10])
#time_utc = np.array([datetime.utcfromtimestamp(t) for t in time_data])



time = np.arange(0,Nseg,1/fe)

#values for BF on 0211
# xmin = 350
# xmax = 500#Nx*xe   # 700 for 0212
# tmin = 130
# tmax = 180

#values for BF on 0212
xmin = 400
xmax = 600#Nx*xe   # 700 for 0212
tmin = 100
tmax = 200




x = np.arange(0,xmax,xe)


# Find indices corresponding to limits
ix_min = np.searchsorted(x, xmin)
ix_max = np.searchsorted(x, xmax)
it_min = np.searchsorted(time, tmin)
it_max = np.searchsorted(time, tmax)

# Slice the data
data_sub = data[ix_min:ix_max, it_min:it_max]
x_sub = x[ix_min:ix_max]
time_sub = time[it_min:it_max]

T, X = np.meshgrid(time_sub, x_sub)  # shape = data_sub.shape

# Flatten the grids to list of coordinates
points = np.vstack((T.ravel(), X.ravel())).T



fe_new = 10
xe_new = 2

time_factor = int(fe/fe_new)
if xe>xe_new:
    space_factor = int(xe/xe_new)
else:
    space_factor = int(xe_new/xe)
# Décimation temporelle (antialias filter inclus)
data_dec_time = decimate(data_sub, time_factor, axis=1, ftype='iir', zero_phase=True)
time_dec = time_sub[::time_factor]

# Décimation spatiale 
data_dec = decimate(data_dec_time, space_factor, axis=0, ftype='iir', zero_phase=True)
x_dec = x_sub[::space_factor]
vmin = np.percentile(data_dec, 5)
vmax = np.percentile(data_dec, 95)



xmid = np.mean(x_dec)
distances = []
for ista in range(len(x_dec)):
    # Calculate the distance between the receiver and the center of the array
    dist = np.sqrt((x_dec[ista]-xmid)**2)
    distances.append(dist)

min_distance_index = np.argmin(distances)
min_distance = distances[min_distance_index]
# print(f"Index of closest station to array centre is {min_distance_index}, and its distance from centre is {min_distance:.2f} meters")







def calculate_bearing(center_coords, point_coords):
    lat1, lon1 = np.radians(center_coords)
    lat2, lon2 = np.radians(point_coords)

    delta_lon = lon2 - lon1

    y = np.sin(delta_lon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(delta_lon)

    # Calculate the bearing in radians
    bearing_rad = np.arctan2(x, y)

    # Convert bearing from radians to degrees
    bearing_deg = np.degrees(bearing_rad)

    # Ensure the angle is between 0 and 360 degrees
    bearing_deg = (bearing_deg + 360) % 360

    return bearing_deg







c_values = np.linspace(5, 20, 32)
theta_inc_values = np.linspace(-np.pi/12, np.pi/12, 60)







# Create a 3D array to store delays for different combinations of 'c', 'theta_inc', and stations
delays = np.empty((len(x_dec), len(c_values), len(theta_inc_values)))

# Iterate over receivers positions, 'c' values, and 'theta_inc' values to calculate delays
for idx in range(len(x_dec)):
    angle_rad = np.radians(calculate_bearing([xmid, 0], [x_dec[idx], 0]))-np.pi/2
    for i, c in enumerate(c_values):
        for j, theta_inc in enumerate(theta_inc_values):
            delay = distances[idx] * np.cos(theta_inc - angle_rad) / c
            delays[idx, i, j] = delay
            
            
            


# def max_by_index(idx, arr):
#     return (idx,) + np.unravel_index(np.argmax(arr[idx]), arr.shape[1:])
# x_by_index(0, delays))

n = data_dec.shape[1]
frequencies = np.fft.fftfreq(n, 1 / (fe/time_factor))
data_spectrum = fft(data_dec, axis=1)

plt.plot(frequencies, np.abs(np.mean(data_spectrum,axis=0)))



# Initialize an empty array to store the delayed data spectra
delayed_data_spectrum = np.empty((len(x_dec), len(c_values), len(theta_inc_values), data_spectrum.shape[1]), dtype=np.complex128)

delayed_data_spectrum.shape
for idx in range(len(x_dec)):
    angle_rad = np.radians(calculate_bearing([xmid, 0], [x_dec[idx], 0]))
    for i, c in enumerate(c_values):
        for j, theta_inc in enumerate(theta_inc_values):
            # Perform the operation and store the result in the array
            delayed_data_spectrum[idx, i, j,:] = data_spectrum[idx] * np.exp(-1j * 2 * np.pi * frequencies * delays[idx, i, j])

delayed_data = np.real(ifft(delayed_data_spectrum, axis=3))


# Initialize variables to store the maximum value and corresponding i and j
max_value = -np.inf
max_i = None
max_j = None

# Initialize the beamformer matrix
beamformer = np.zeros((len(c_values), len(theta_inc_values)))

# Iterate over all combinations of i and j
for i in range(len(c_values)):
    for j in range(len(theta_inc_values)):
        # Calculate the quantity
        quantity = np.sum(np.sum(delayed_data[:, i, j, :], axis=0) ** 2)

        # Update the beamformer matrix
        beamformer[i, j] = quantity

        # Check if the current quantity is the maximum
        if quantity > max_value:
            max_value = quantity
            max_i = i
            max_j = j


C, Theta = np.meshgrid(c_values, theta_inc_values)

# Calculate X and Y using element-wise multiplication
X = C * np.cos(Theta)
Y = C * np.sin(Theta)

# Create the meshgrid
C, Theta = np.meshgrid(c_values, theta_inc_values)

# Calculate X and Y using element-wise multiplication
X = C * np.cos(Theta)
Y = C * np.sin(Theta)

# Create a 2D plot of the beamformer
plt.figure(figsize=(8, 8))
contour = plt.contourf(X, Y, beamformer.T, cmap='viridis')
plt.colorbar(contour, label='Beamformer Value')

# Marker for maximum
X_max = c_values[max_i] * np.cos(theta_inc_values[max_j])
Y_max = c_values[max_i] * np.sin(theta_inc_values[max_j])
plt.plot(X_max, Y_max, 'ro', markersize=10, label='Max Beamformer')

# Annotate with values
plt.text(X_max, Y_max, 
         f" {c_values[max_i]:.1f} m/s\n {theta_inc_values[max_j]*180/np.pi:.1f}°", 
         color='red', fontsize=10, ha='left', va='bottom')

# Set labels and title
plt.xlabel('velocity projected on X axis')
plt.ylabel('velocity projected on Y axis')
plt.title('2D Plot of Beamformer ')
plt.legend()

plt.show()






fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # 1 row, 2 columns

# --- First subplot ---
im0 = axes[0].imshow(data_dec, aspect='auto', cmap='gray', origin='lower',
                     extent=[time_dec[0], time_dec[-1], x_dec[0], x_dec[-1]],
                     vmin=vmin, vmax=vmax)
fig.colorbar(im0, ax=axes[0], label='Strain Rate [nStrain/s]')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('Distance (m)')
axes[0].set_title('strain avant BF')

# --- Second subplot ---
delayed_dta_BF = delayed_data[:, max_i, max_j, :]
x_dec = x_sub[::space_factor]
vmin = np.percentile(delayed_dta_BF, 5)
vmax = np.percentile(delayed_dta_BF, 95)

im1 = axes[1].imshow(delayed_dta_BF, aspect='auto', cmap='gray', origin='lower',
                     extent=[time_dec[0], time_dec[-1], x_dec[0], x_dec[-1]],
                     vmin=vmin, vmax=vmax)
fig.colorbar(im1, ax=axes[1], label='Strain Rate [nStrain/s]')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Distance (m)')
axes[1].set_title('strain après BF')

plt.tight_layout()
plt.show()








