import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

def get_path2dailydir(disk, date, serie):
    if type(serie)==int:
        path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}/serie{serie}'
    else:
        path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}'
    
    return path2dailydir

def find_acqname(path2dailydir, acq):
    lisdir = os.listdir(path2dailydir)
    for acqname in lisdir:
        if f'acq{acq}_' in acqname:
            acqname_found = acqname
    return acqname_found

def load_force_vs_frames_1balance(idx, date, serie, acq, data_csv, csv_filename_px_balance='pixel_balance_vs_frame.csv',
                                  csv_filename_px_ref = 'pixel_vs_frame_ref.csv', disk='D:',plot=False):
    frame_frac = data_csv[idx, 12]  

    path2dailydir = get_path2dailydir(disk, date, serie)

    acqname_found = find_acqname(path2dailydir, acq)

    path2csv = f'{path2dailydir}/{acqname_found}/{csv_filename_px_balance}'

    pd_data_forcepx = pd.read_csv(path2csv)

    frames = pd_data_forcepx['Y'].to_numpy()
    frames_int = [np.round(float(f)).astype(int) for f in frames]
    frames_int = np.array(frames_int)

    x_px = pd_data_forcepx['X'].to_numpy()
    x_px = [np.round(float(j)).astype(int) for j in x_px]
    x_px = np.array(x_px)

    # si la ref bouge un peu, on soustrait la ref :
    if os.path.exists(f'{path2dailydir}/{acqname_found}/{csv_filename_px_ref}'):
        pd_data_pxref = pd.read_csv(f'{path2dailydir}/{acqname_found}/{csv_filename_px_ref}')
        x_px_ref = pd_data_pxref['X'].to_numpy()
        frames_ref = pd_data_pxref['Y'].to_numpy()
        print(x_px_ref)

        x_px_ref_interp = np.interp(frames, frames_ref, x_px_ref)

        x_px = x_px - x_px_ref_interp

    if csv_filename_px_balance=='pixel_balance_vs_frame.csv':
        # and load the dg/dpx also:
        dg_sur_dpx = float(data_csv[idx, 13])
    
    elif csv_filename_px_balance.endswith('g.csv'):
        # par exemple pixel_balance_vs_frame_2500g.csv           
        match = re.search(r'(\d+)g', csv_filename_px_balance)
        if not match:
            print('Dont find the dg_sur_dpx file !!')
        nombre = match.group(1)
        txt_file_dg_sur_dpx = f"{path2dailydir}/{acqname_found}/dg_sur_dpx_{nombre}g.txt"
        dg_sur_dpx = float(np.loadtxt(txt_file_dg_sur_dpx))
    else:
        print('csv filename px balance ??')
        dg_sur_dpx = np.nan
        
    dm_g = dg_sur_dpx
    dx_px = 1

    # assume the first frame is for m0 = 0 grams

    m0 = (x_px * (dm_g/dx_px))[0]
    
    if plot:
        plt.figure(figsize=(12,9))
        plt.plot(frames_int, x_px * (dm_g/dx_px) - m0, '.')
        plt.xlabel('frames')
        plt.ylabel('force applied (in grams)')
        plt.show()

    # save data (force vs frames):
    dict_results_1balance = {}
    dict_results_1balance['frames'] = frames_int
    dict_results_1balance['x_px'] = x_px
    dict_results_1balance['force_grams'] = x_px * (dm_g/dx_px) - m0
    dict_results_1balance['dg_sur_dpx'] = dm_g/dx_px
    dict_results_1balance['m0'] = m0 # ref par rapport à laquelle la masse est mesurée

    return dict_results_1balance


def load_force_vs_frames_2balances(idx, date, serie, acq, data_csv, disk='D:',plot=False):
        # cas où il y avait beosin de plusieurs balances
    list_balances_maxweights = []
    list_dg_sur_dpx = []

    path2dailydir = get_path2dailydir(disk, date, serie)
    acqname_found = find_acqname(path2dailydir, acq)

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

    return dict_results