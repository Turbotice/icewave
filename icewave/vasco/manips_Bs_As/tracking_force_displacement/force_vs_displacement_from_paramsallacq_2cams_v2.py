# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 18:21:41 2025

@author: Vasco Zanchi
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.io import loadmat
import pickle
import sys
from scipy.optimize import curve_fit


sys.path.append('functions')
from utils import load_csv_params_allacq, str_to_float, correspond_samplenum_acqnum, reshape_array, clean_array, load_displacement_data_PIV
from load_force_vs_frame import *
from load_displacement_vs_frame import load_disp_vs_frame

from force_disp import load_force_displacement_curve

%matplotlib inline

#%% Load csv data "params_all_acq.csv"


disk = 'D:'
csv_file_path = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/params_acquisitions/params_allacq.csv"

data_csv = load_csv_params_allacq(csv_file_path=csv_file_path)


# Apply the conversion to every element
#data_csv = np.vectorize(str_to_float)(data_csv)

data_csv_converted = np.zeros(data_csv.shape, dtype=object)  # Create an empty array with the same shape
for i in range(data_csv.shape[0]):
    for j in range(data_csv.shape[1]):
        data_csv_converted[i, j] = str_to_float(data_csv[i, j])
data_csv = data_csv_converted
#print(data_csv_converted)

#%% PART 1 : displacement vs frames
idx = 68
method = data_csv[idx, 25]

use_ref = int(data_csv[idx, 20])

frame_frac = data_csv[idx, 12]
frame_force_application = data_csv[idx,11] #début de l'application de la force

frame_ref_cam1 = data_csv[idx, 16]
frame_ref_cam2 = data_csv[idx, 17]

date = str(int(data_csv[idx, 0])).zfill(4)
acq = int(data_csv[idx, 1])
serie = data_csv[idx, 2] # si pas de serie, mettre serie = None

# tracking using coordonnées des vis

xpx_vis = int(data_csv[idx, 3])
ypx_vis = int(data_csv[idx, 4])
xpx_vis_ref = int(data_csv[idx, 5])
ypx_vis_ref = int(data_csv[idx, 6])

if np.isnan(serie):
    serie = None
else:
    serie = int(serie)

if np.isnan(frame_ref_cam2):
    pass
else:
    frameref_diff = frame_ref_cam2 - frame_ref_cam1
    frame_frac += frameref_diff

i0 = int(data_csv[idx, 7])
if method=='PIV':
    refimg = int(data_csv[idx, 8])
    W = int(data_csv[idx, 9])
    n_passages = int(data_csv[idx, 10])
else:
    refimg = np.nan
    W = np.nan
    n_passages = np.nan

dmm_sur_dpx = float(data_csv[idx, 14])
dmm = dmm_sur_dpx
dpx = 1


path2dailydir = get_path2dailydir(disk,date,serie)
acqname_found = find_acqname(path2dailydir,acq)


dict_results = load_disp_vs_frame(path2dailydir, acqname_found, method=method,
                                  n_passages=n_passages,i0=i0, refimg=refimg, W=W,
                                    xpx_vis=xpx_vis, ypx_vis=ypx_vis, xpx_vis_ref=xpx_vis_ref,
                                      ypx_vis_ref=ypx_vis_ref, use_ref=use_ref,
                                        frame_force_application=frame_force_application, 
                                        frame_frac=frame_frac, dmm=dmm,dpx=dpx)

force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

if serie==None:
    with open(force_disp_track_dir + f'data_displacement_{date}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)
else:
    with open(force_disp_track_dir + f'data_displacement_{date}_serie{serie}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)



#%% PART 2 : force vs frames


from functions.load_force_vs_frame import *

#dict_results = load_force_vs_frames_1balance(idx=idx, date=date, serie=serie, acq=acq, data_csv=data_csv, plot=True)

if os.path.exists(f'{path2dailydir}/{acqname_found}/pixel_balance_vs_frame.csv'):
    # cas où il y avait besoin d'une seule balance
    dict_results = load_force_vs_frames_1balance(idx=idx, date=date, serie=serie, acq=acq, data_csv=data_csv)
else:
    dict_results = load_force_vs_frames_2balances(idx=idx, date=date, serie=serie, acq=acq, data_csv=data_csv)


force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

if serie==None:
    with open(force_disp_track_dir + f'data_force_{date}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)
else:
    with open(force_disp_track_dir + f'data_force_{date}_serie{serie}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)

#%% PART 3 : force vs displacement


force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

if serie==None:
    with open(force_disp_track_dir + f'data_force_{date}_acq{acq}.pkl', 'rb') as file:
        dict_force = pickle.load(file)
else:
    with open(force_disp_track_dir + f'data_force_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
        dict_force = pickle.load(file)

if serie==None:
    with open(force_disp_track_dir + f'data_displacement_{date}_acq{acq}.pkl', 'rb') as file:
        dict_displacement = pickle.load(file)
else:
    with open(force_disp_track_dir + f'data_displacement_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
        dict_displacement = pickle.load(file)
        
# plot force vs displacement
if np.isnan(frame_ref_cam1):
    common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames'], return_indices=True)
else:
    common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames']+frameref_diff, return_indices=True)

common_old = common
common = common_old[common>=i0]
idx_a = idx_a[common_old>=i0]
idx_b = idx_b[common_old>=i0]

print("Common values:", common)
print("Indices in a:", idx_a)
print("Indices in b:", idx_b)


u0_relative = dict_displacement['u_field_relative_mm'][idx_a][0]
frame_frac = dict_displacement['frame_frac']

idx_frac_common = np.where(common==frame_frac)[0][0]

plt.figure()
plt.plot(dict_force['force_grams'][idx_b], dict_displacement['u_field_relative_mm'][idx_a] - u0_relative,'.')
#plt.xlim(0, ) # limites à bien choisir pour ensuite bien fitter
plt.ylim(-1e-1,1)
plt.xlabel('force (grams)')
plt.ylabel('displacement (mm)')
plt.savefig(force_disp_track_dir + 'f_grams_vs_displacement_mm.pdf',dpi=600)
plt.show()

# convert force in newtons and displacement in meters
force_newtons_comidx = dict_force['force_grams'][idx_b] * 1e-3 * 9.81
displacement_meters_comidx = 1e-3 * (dict_displacement['u_field_relative_mm'][idx_a] - u0_relative)

xdata2fit = force_newtons_comidx[:idx_frac_common]
ydata2fit = displacement_meters_comidx[:idx_frac_common]

xdata2fit = np.array(xdata2fit, dtype=float)
ydata2fit = np.array(ydata2fit, dtype=float)


try:
    popt, pcov = curve_fit(lambda x,a,b: a*x+b , xdata2fit, ydata2fit)
except:
    popt, pcov = np.array([np.nan,np.nan]), np.ones((2,2))*np.nan
slope, intercept = popt[0], popt[1]
slope_err, intercept_err = np.diag(pcov)[0], np.diag(pcov)[1]

plt.figure()
plt.plot(xdata2fit, ydata2fit,'.')
#plt.xlim(0, ) # limites à bien choisir pour ensuite bien fitter
#plt.ylim(-1e-1,1)
plt.plot(xdata2fit, slope*xdata2fit+intercept)
plt.xlabel('Fz (N)')
plt.ylabel('displacement (m)')
plt.savefig(force_disp_track_dir + f'f_N_vs_displacement_m_{date}_serie{serie}_acq{acq}.pdf',dpi=600)
plt.show()

# linear fit of force vs displacement
# F/displacement = E * 4*w*(h/L)**3

#if serie!=None:
#    fd = ((serie-1) * 3+1) + (acq-1)//3
#else:
#    fd = (acq-1)//3 + 1
#sd=acq%3
#if sd==0:
#    sd=3


fd, sd = correspond_samplenum_acqnum(date=date, acq=acq, serie=serie)

print(f'acq{acq} and serie{serie} correspond to sample name : M{fd}{sd}')

prefix_filename_thickness = f'M{fd}{sd}'


thickness_path1 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}.txt"
thickness_path2 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}_camera.txt"
thickness_path3 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}_post_frac.txt"
if os.path.exists(thickness_path1):
    thickness_path = thickness_path1
elif os.path.exists(thickness_path2):
    thickness_path = thickness_path2
else:
    thickness_path = thickness_path3

data_thickness = np.loadtxt(thickness_path,skiprows=1)
print('thicknesses measured for this sample:',data_thickness)
h_avg_mm = np.mean(data_thickness)
h_std_mm = np.std(data_thickness)

h_avg = 1e-3 * h_avg_mm
h_std = 1e-3 * h_std_mm

# compute young's modulus with the formula for elastic beam

w = 4e-2
L = data_csv[idx, 15]

E_beam = (1/(4*w))*(L/h_avg)**3 * (1/slope)
#E_err = (slope_err/slope) * E_beam

E_beam_upp = (1/(4*w))*(L/(h_avg-h_std))**3 * (1/slope)

E_err = E_beam_upp - E_beam

print('E = ',E_beam*1e-9,' GPa')
print('E_err = ',E_err*1e-9,' GPa')
print('E_upp = ',E_beam_upp*1e-9,' GPa')

print('E_upp - E = ',(E_beam_upp-E_beam)*1e-9,' GPa')

# %%
%matplotlib inline
N = len(data_csv[:,0])
print(N)
arr_h_avg = np.zeros(N)
arr_h_std = np.zeros(N)
arr_E = np.zeros(N)
arr_Eerr = np.zeros(N)
arr_prefix_filename_thickness = np.empty(N,dtype=object)
arr_dates = np.empty(N,dtype=object)


import pandas as pd
df = pd.read_csv(csv_file_path, delimiter=';')
qualite_fit = df['qualite_fit'].values

for i in range(N):
    print('index',i)
    if (qualite_fit[i]!='ok') and (qualite_fit[i]!='moyen'):
        _, arr_h_avg[i], arr_h_std[i], arr_E[i], arr_Eerr[i],arr_prefix_filename_thickness[i],arr_dates[i] = np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
    else:
        print(qualite_fit[i])
        _, arr_h_avg[i], arr_h_std[i], arr_E[i], arr_Eerr[i],arr_prefix_filename_thickness[i],arr_dates[i] = load_force_displacement_curve(idx=i,data_csv=data_csv,plot=True)


#mask = (arr_Eerr/arr_E < 1/2) & (arr_E>0) # on garde que les points pour lesquels l'erreur sur E est inférieure à 50% et E positif
mask = ((qualite_fit=='ok') | (qualite_fit=='moyen'))&((arr_Eerr/arr_E < 1/2) & (arr_E>0))

plt.figure()
plt.errorbar(arr_h_avg[mask], 1e-9 * arr_E[mask],yerr=1e-9 * np.abs(arr_Eerr[mask]),linestyle='',ecolor='k',color='b',marker='o')
plt.ylim(0,15)
plt.xlim(0,8e-3)
plt.xlabel('sample thickness (m)')
plt.ylabel('E (GPa)')
plt.show()


# %%
