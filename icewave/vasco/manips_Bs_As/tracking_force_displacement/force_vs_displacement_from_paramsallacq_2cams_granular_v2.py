# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 18:21:41 2025

@author: Vasco Zanchi
"""
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import re
from scipy.io import loadmat
import pickle
from scipy.optimize import curve_fit
import sys

%matplotlib qt

#%%


disk = 'D:'
csv_file_path = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/params_acquisitions/params_allacq_granular.csv"

with open(csv_file_path, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=';')
    count=0
    data = []
    for lines in csvFile:
        if count==0:
            header = lines
        else:
            data.append(lines)
        count+=1
print(header)

data_csv = np.array(data)

# --- Convert all string entries to float ---
def str_to_float(value):
    try:
        # Replace comma with dot (European decimal format)
        value = value.replace(',', '.')
        # Remove spaces just in case
        value = value.strip()
        # Convert to float
        return float(value)
    except ValueError:
        # If conversion fails, return NaN
        return np.nan

# Apply the conversion to every element
data_csv = np.vectorize(str_to_float)(data_csv)

print(data_csv)

#%% PART 1 : displacement vs frames
idx = 5

use_ref = int(data_csv[idx, 20])

frame_frac = data_csv[idx, 12]
frame_force_application = data_csv[idx,11] #début de l'application de la force

frame_ref_cam1 = data_csv[idx, 16]
frame_ref_cam2 = data_csv[idx, 17]

date = str(int(data_csv[idx, 0])).zfill(4)
acq = int(data_csv[idx, 1])
serie = data_csv[idx, 2] # si pas de serie, mettre serie = None

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
refimg = int(data_csv[idx, 8])
W = int(data_csv[idx, 9])
n_passages = int(data_csv[idx, 10])


if type(serie)==int:
    path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}/serie{serie}'
else:
    path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}'

lisdir = os.listdir(path2dailydir)
for acqname in lisdir:
    if f'acq{acq}_' in acqname:
        acqname_found = acqname

matfile = f'{path2dailydir}/{acqname_found}/PIV_processed_{n_passages}passages_i0{i0}_refimg{refimg}_W{W}.mat'

mat_dict = loadmat(matfile)

# deux facons d'importer les données en fonction de la structure du fichier .mat
def reshape_array(arr):
    array_new = np.zeros((len(arr),arr[0].shape[0],arr[0].shape[1]))
    for i in range(len(arr)):
        array_new[i,:,:] = arr[i]
    return array_new

def clean_array(arr):
    # On va parcourir le tableau et garder seulement les sous-tableaux qui ont au moins un élément non nul ou non NaN
    cleaned = []

    for sublist in arr:
        for subarr in sublist:
            # Vérifie si le tableau contient au moins un élément significatif
            if np.any(~np.isnan(subarr) & (subarr != 0)):
                cleaned.append(subarr)
    # Convertir en tableau numpy propre
    cleaned_array = np.array(cleaned, dtype=object)
    return cleaned_array

if 'u_original' in mat_dict:
    u_original = mat_dict['u_original'][:,0]
    v_original = mat_dict['v_original'][:,0]
    u = reshape_array(u_original)
    v = reshape_array(v_original)
    xpix = mat_dict['x'][0][0][0]
    ypix = mat_dict['y'][0][0][:,0]

    concat_with_zero = True
else:
    #u = np.moveaxis(clean_array(mat_dict['u']), 2, 1)
    #v = np.moveaxis(clean_array(mat_dict['v']), 2, 1)
    u = clean_array(mat_dict['u'])
    v = clean_array(mat_dict['v'])

    xpix = mat_dict['xpix'].flatten()
    ypix = mat_dict['ypix'].flatten()
    concat_with_zero = False


# pour ajouter un zero au début (dans le cas imgref et piv faite à la mains, car  sinon ca décale par rapport aux mesures de forces)
if concat_with_zero:
    u = np.concatenate((np.zeros((1,u.shape[1],u.shape[2])), u))
    v = np.concatenate((np.zeros((1,v.shape[1],v.shape[2])), v))


# tracking using coordonnées des vis

xpx_vis = int(data_csv[idx, 3])
ypx_vis = int(data_csv[idx, 4])
xpx_vis_ref = int(data_csv[idx, 5])
ypx_vis_ref = int(data_csv[idx, 6])

xind = (np.abs(xpix - xpx_vis)).argmin()
yind = (np.abs(ypix - ypx_vis)).argmin()

xind_ref = (np.abs(xpix - xpx_vis_ref)).argmin()
yind_ref = (np.abs(ypix - ypx_vis_ref)).argmin()


dmm_sur_dpx = float(data_csv[idx, 14])
dmm = dmm_sur_dpx
dpx = 1

nb_neighbors2avg = 0

# average neighbours and smooth vs time :
#u_field = gaussian_filter1d(np.mean(u[:,yind-2:yind+2,xind-2:xind+2],axis=(1,2)),5)
if nb_neighbors2avg!=0:
    u_field = np.nanmean(u[:,yind-nb_neighbors2avg:yind+nb_neighbors2avg,xind-nb_neighbors2avg:xind+nb_neighbors2avg],axis=(1,2))
    u_field_ref = np.nanmean(u[:,yind_ref-nb_neighbors2avg:yind_ref+nb_neighbors2avg,xind_ref-nb_neighbors2avg:xind_ref+nb_neighbors2avg],axis=(1,2))
else:
    u_field = u[:,yind,xind]
    u_field_ref = u[:,yind_ref,xind_ref]


if use_ref:
    u_field_relative = u_field - u_field_ref
#elif use_ref_edgeplate:
#    u_field_relative = u_field - u_field_ref_edge
else:
    u_field_relative = u_field


plt.figure()
plt.plot(np.arange(u.shape[0])+i0,u_field, label='moving point')
plt.plot(np.arange(u.shape[0])+i0,u_field_ref, label='ref (point immobile)')
plt.plot(np.arange(u.shape[0])+i0,u_field_relative, label='relative motion')

plt.xlim(i0,frame_frac)
#plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
plt.legend()
plt.show()


plt.figure()
plt.plot(np.arange(u.shape[0])+i0,u_field * dmm/dpx, label='moving point')
plt.plot(np.arange(u.shape[0])+i0,u_field_ref* dmm/dpx, label='ref (point immobile)')
plt.plot(np.arange(u.shape[0])+i0,u_field_relative* dmm/dpx, label='relative motion')

plt.xlim(i0,frame_frac)
#plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
plt.xlabel('t (frames) ')
plt.ylabel('displacement (mm)')
plt.legend()
plt.show()

# save pkl for displacement vs frames



dict_results = {}
dict_results['u_field'] = u_field
dict_results['u_field_mm'] = u_field * dmm/dpx
dict_results['u_field_ref'] = u_field_ref
dict_results['u_field_ref_mm'] = u_field_ref * dmm/dpx
dict_results['u_field_relative'] = u_field_relative
dict_results['u_field_relative_mm'] = u_field_relative * dmm/dpx
dict_results['i0'] = i0
dict_results['refimg'] = refimg
dict_results['frame_frac'] = frame_frac
dict_results['frame_force_application'] = frame_force_application

dict_results['xind'] = xind
dict_results['yind'] = yind
dict_results['xind_ref'] = xind_ref
dict_results['yind_ref'] = yind_ref
dict_results['dmm_sur_dpx'] = dmm/dpx
dict_results['xpix'] = xpix
dict_results['ypix'] = ypix
dict_results['u'] = u
dict_results['v'] = v
dict_results['frames_piv'] = np.arange(u.shape[0])+i0

force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

if serie==None:
    with open(force_disp_track_dir + f'data_displacement_{date}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)
else:
    with open(force_disp_track_dir + f'data_displacement_{date}_serie{serie}_acq{acq}.pkl', 'wb') as file:
        pickle.dump(dict_results, file)



#%% PART 2 : force vs frames

#csv_filename_px_single_balance = 'pixel_balance_vs_frame.csv'
#csv_filename_px_ref = 'pixel_vs_frame_ref.csv' # this file doesn't necessarily exist

from functions.load_force_vs_frame import load_force_vs_frames_1balance



if os.path.exists(f'{path2dailydir}/{acqname_found}/pixel_balance_vs_frame.csv'):
    # cas où il y avait besoin d'une seule balance
    dict_results = load_force_vs_frames_1balance(idx=idx, date=date, serie=serie, acq=acq, data_csv=data_csv)
else:
    # cas où il y avait beosin de plusieurs balances
    list_balances_maxweights = []
    list_dg_sur_dpx = []
    for fname in os.listdir(f'{path2dailydir}/{acqname_found}'):
        if ('g.txt' in fname)&('dg_sur_dpx_' in fname):
            list_dg_sur_dpx.append(fname)
            list_balances_maxweights.append(int(re.search(r'(\d+)g', fname).group(1)))
    
    dict_results = {}
    for i in range(len(list_dg_sur_dpx)):
        dict_results[f'balance_{list_balances_maxweights[i]}g'] = load_force_vs_frames_1balance(idx=idx, date=date, serie=serie, acq=acq, data_csv=data_csv, csv_filename_px_balance=f'pixel_balance_vs_frame_{list_balances_maxweights[i]}g.csv',plot=True)

    if len(list_dg_sur_dpx)==2:
        common, idx_a, idx_b = np.intersect1d(dict_results[f'balance_{list_balances_maxweights[0]}g']['frames'],dict_results[f'balance_{list_balances_maxweights[1]}g']['frames'], return_indices=True)
        frames_common = common
        force_grams = dict_results[f'balance_{list_balances_maxweights[0]}g']['force_grams'][idx_a] + dict_results[f'balance_{list_balances_maxweights[1]}g']['force_grams'][idx_b]
        dict_results['frames'] = frames_common
        dict_results['force_grams'] = force_grams
    else:
        print("Code non adapté pour le cas où il y a + de 2 balances !")

# et on plot la force grams vs frames (au cas où on a fait la somme de 2 balances temps)
plt.figure(figsize=(12,9))
plt.plot(dict_results['frames'], dict_results['force_grams'], '.')
plt.xlabel('frames')
plt.ylabel('force applied (in grams)')
plt.show()


# partie pour enregistrer dict
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

#xdata2fit = force_newtons_comidx[:idx_frac_common-60]
#ydata2fit = displacement_meters_comidx[:idx_frac_common-60]

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

if serie!=None:
    fd = ((serie-1) * 3+1) + (acq-1)//3
else:
    fd = (acq-1)//3 + 1
sd=acq%3
if sd==0:
    sd=3
print(fd, sd)
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
# %% on met tout ca dans une fonction load_force_displacement_curve()


def load_force_displacement_curve(idx=0, disk='D:', data_csv=data_csv, plot=True):

    
    use_ref = int(data_csv[idx, 20])

    frame_frac = data_csv[idx, 12]
    frame_force_application = data_csv[idx,11] #début de l'application de la force

    frame_ref_cam1 = data_csv[idx, 16]
    frame_ref_cam2 = data_csv[idx, 17]

    date = str(int(data_csv[idx, 0])).zfill(4)
    print(date)
    acq = int(data_csv[idx, 1])
    serie = data_csv[idx, 2] # si pas de serie, mettre serie = None

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
    refimg = int(data_csv[idx, 8])
    W = int(data_csv[idx, 9])
    n_passages = int(data_csv[idx, 10])

    force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

    if serie==None:
        with open(force_disp_track_dir + f'data_force_{date}_acq{acq}.pkl', 'rb') as file:
            dict_force = pickle.load(file)
            print('dict force loaded')
    else:
        with open(force_disp_track_dir + f'data_force_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
            dict_force = pickle.load(file)
            print('dict force loaded')

    if serie==None:
        with open(force_disp_track_dir + f'data_displacement_{date}_acq{acq}.pkl', 'rb') as file:
            dict_displacement = pickle.load(file)
            print('dict displacement loaded')
    else:
        with open(force_disp_track_dir + f'data_displacement_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
            dict_displacement = pickle.load(file)
            print('dict displacement loaded')
            
    # plot force vs displacement
    if np.isnan(frame_ref_cam1):
        common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames'], return_indices=True)
    else:
        common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames']+frameref_diff, return_indices=True)

    print("Common values:", common)
    print("Indices in a:", idx_a)
    print("Indices in b:", idx_b)


    u0_relative = dict_displacement['u_field_relative_mm'][idx_a][0]
    frame_frac = dict_displacement['frame_frac']

    idx_frac_common = np.where(common==frame_frac)[0][0]
    if plot:
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
    if plot:
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

    if serie!=None:
        fd = ((serie-1) * 3 + 1) + acq//3
    else:
        fd = (acq)//3
    sd = acq%3
    if sd==0:
        sd=3
    if fd==0:
        fd = 1
    print(fd)
    print(sd)
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
    E_err_fit = (slope_err/slope) * E_beam # ca serait l'incertitude liée au fit, qui est faible par rapport à celle liée à l'épaisseur

    E_beam_upp = (1/(4*w))*(L/(h_avg-h_std))**3 * (1/slope)

    E_err_h = E_beam_upp - E_beam

    if E_err_fit>E_err_h:
        E_err = E_err_fit
    else:
        E_err = E_err_h

    print('E = ',E_beam*1e-9,' GPa')
    print('E_upp = ',E_beam_upp*1e-9,' GPa')

    print('E_upp - E = ',(E_beam_upp-E_beam)*1e-9,' GPa')

    print(idx)

    return (h_avg, h_std, E_beam, E_err,prefix_filename_thickness,date)



# %%
N = 5
arr_h_avg = np.zeros(N)
arr_h_std = np.zeros(N)
arr_E = np.zeros(N)
arr_Eerr = np.zeros(N)
arr_prefix_filename_thickness = np.empty(N,dtype=object)
arr_dates = np.empty(N,dtype=object)


for i in range(N):
    arr_h_avg[i], arr_h_std[i], arr_E[i], arr_Eerr[i],arr_prefix_filename_thickness[i],arr_dates[i] = load_force_displacement_curve(idx=i,plot=False)
# %%

mask = (arr_Eerr/arr_E < 1/3) & (arr_E>0)

plt.figure()
plt.errorbar(arr_h_avg[mask], 1e-9 * arr_E[mask],yerr=1e-9 * np.abs(arr_Eerr[mask]),linestyle='',ecolor='k',color='b',marker='o')
plt.ylim(0,11)
plt.xlabel('sample thickness (mm)')
plt.ylabel('E (GPa)')
plt.show()
# %%
def load_grain_sizes_byhand(date,sample_name,disk=disk):
    path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{date}/{sample_name}_grains_sizes_mm.txt"
    datagrainsizes = np.loadtxt(path2grainsizes) 

    return datagrainsizes, np.mean(datagrainsizes), np.std(datagrainsizes)


plt.figure()
for i in range(len(arr_E)):
    sample_name = arr_prefix_filename_thickness[i]
    date = arr_dates[i]
    path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{date}/{sample_name}_grains_sizes_mm.txt"
    if os.path.exists(path2grainsizes):
        print('exists')
        datagrainsizes,grainsizeavg,grainsizestd = load_grain_sizes_byhand(date, sample_name)
    if (mask[i] & os.path.exists(path2grainsizes)):
        plt.errorbar(grainsizeavg, arr_E[i]*1e-9,xerr=grainsizestd,yerr=1e-9*arr_Eerr[i],marker='o',linestyle='',ecolor='k',color='b')
plt.ylim(0,11)
plt.xlim(0,5)
plt.show()
# %%
