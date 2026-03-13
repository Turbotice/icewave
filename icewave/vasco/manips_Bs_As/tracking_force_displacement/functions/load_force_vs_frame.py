import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os



def load_force_vs_frames_1balance(idx, date, serie, acq, data_csv, csv_filename_px_balance='pixel_balance_vs_frame.csv',
                                  csv_filename_px_ref = 'pixel_vs_frame_ref.csv', disk='D:',plot=False):
    frame_frac = data_csv[idx, 12]  

    if type(serie)==int:
        path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}/serie{serie}'
    else:
        path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}'

    lisdir = os.listdir(path2dailydir)
    for acqname in lisdir:
        if f'acq{acq}_' in acqname:
            acqname_found = acqname

    path2csv = f'{path2dailydir}/{acqname_found}/{csv_filename_px_balance}'
    """
    with open(path2csv, mode ='r')as file:
        csvFile = csv.reader(file, delimiter=',')
        count=0
        data_forcepx = []
        for lines in csvFile:
            if count==0:
                header = lines
            else:
                data_forcepx.append(lines)
            count+=1
    print(header)

    data_forcepx = np.array(data_forcepx)"""

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