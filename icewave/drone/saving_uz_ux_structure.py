# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 10:32:16 2025

@author: sebas

This script  enables to compute real velocity fields uz(x,t) and ux(x,y,t) from apparent velocity fields 
Vx(x,y,t) and Vy(x,y,t) observed on the camera sensor. 

This can be used only if we assume that the real velocity field has no horizontal component uy. 
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


import h5py
import os
import glob

import sys
sys.path.append('C:/Users/sebas/git')

import icewave.tools.datafolders as df
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.drone.drone_projection as dp
import icewave.tools.rw_data as rw

parula_map = matcmaps.parula()

def load_mat_structure(base,date,drone_ID,exp_ID, idx_file = 0):
    """ Load matlab structure saved under .mat file (-v7.3)
    Inputs : - base, string, basis directory where field data are saved
             - date, string, date of experiment 'mmdd' 
             - drone_ID, string, drone name
             - exp_ID, string, experiment name, for instance 01-waves_001
             - idx_file, int, optional argument, designs file index which has a suffixe '*scaled.mat'
    Outputs : - S, dictionnary containing apparent velocity fields Vx and Vy 
              - file2load, string, file name that has been loaded """
    
    path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'

    if not os.path.isdir(path2data):
        print(f'Folder {path2data} does not exist')
        
    filelist = glob.glob(path2data + '*scaled.mat')
    file2load = filelist[0]
    
    # load data 
    with h5py.File(file2load, 'r') as fmat:
        S = {}
        print('Top-level keys : ', list(fmat.keys()))
        S = mat2py.mat_to_dict(fmat['m'],fmat['m'])
        S = mat2py.transpose_PIVmat_fields(S)

    return S,file2load

def plot_Vx_Vy(Vx,Vy,S,frame,cmap = parula_map):
    """ Plot apparent velocity fields Vx and Vy
    Inputs : - Vx and Vy, array like [nx,ny,nt], apparent velocity along horizontal and vertical axis of camera sensor. """

    extents_meter = np.array([S['x'].min(),S['x'].max(),S['y'].min(),S['y'].max()])
    
    set_graphs.set_matplotlib_param('single')
    fig, axs = plt.subplots(ncols = 2,sharey = True, figsize = (12,9))
    imsh = axs[0].imshow(Vx[:,:,frame].T,cmap = cmap,origin = 'lower',aspect = 'equal',extent = extents_meter)
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)
    cbar.set_label(r'$V_x \; \mathrm{(m.s^{-1})}$')

    imsh = axs[1].imshow(Vy[:,:,frame].T,cmap = cmap,origin = 'lower',aspect = 'equal',extent = extents_meter)
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)
    cbar.set_label(r'$V_y \; \mathrm{(m.s^{-1})}$')

    axs[0].set_ylabel(r'$y \; \mathrm{(m)}$')
    for ax in axs:
        ax.set_xlabel(r'$x \; \mathrm{(m)}$')

    plt.tight_layout()
    

    return fig, axs

def process_file(base,date,drone_ID,exp_ID, flip_field = 0, idx_file = 0):
    """ Load structure obtained from PIVlab, supress quadratic noise, plot the apparent velocity field and compute 
    real velocity fields (uz,ux) vertical and along wave field propagation, and save resulting structure. 
    Inputs : - base, string, basis directory where field data are saved
             - date, string, date of experiment 'mmdd' 
             - drone_ID, string, drone name
             - exp_ID, string, experiment name, for instance 01-waves_001
             - flip_field, boolean, 1 if wave field propagates from right to left on the video 
             - idx_file, int, optional argument, designs file index which has a suffixe '*scaled.mat'
    """
    
    S, file2load = load_mat_structure(base, date, drone_ID, exp_ID, idx_file = idx_file)
    
    Vx = FT.supress_quadratic_noise(np.transpose(S['Vx'],(1,0,2)),S['x'],S['y'])
    Vy = FT.supress_quadratic_noise(np.transpose(S['Vy'],(1,0,2)),S['x'],S['y'])
    Vx = np.transpose(Vx,(1,0,2))
    Vy = np.transpose(Vy,(1,0,2))
    if flip_field:
        Vx = np.flip(Vx,(0,1))
        Vy = np.flip(Vy,(0,1))
    print('Quadratic field supressed')
    
    frame = 0
    fig, axs = plot_Vx_Vy(Vx, Vy, S, frame)
    
    uz,ux,err_uz = dp.get_uz_ux_from_structure(Vx,Vy,S)

    data = S.copy()
    data['uz'] = uz
    data['ux'] = ux
    data['err_uz'] = err_uz
    del data['Vx']
    del data['Vy']
    del data['PIV_param']['p_param']
    del data['PIV_param']['s_param']
    
    mat_name = file2load.split('\\')[-1]
    filename = mat_name.replace('.mat','.h5')
    filename = filename.replace('PIV_processed','uz_ux')
    path2data = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
    
    print('Saving data...')
    rw.save_dict_to_h5(data, f'{path2data}{filename}')
    print(f'uz and ux structure saved under name : {path2data}{filename}')
    
    return 


def main():
    base = 'E:/Rimouski_2024/Data/'

    date = '0226'
    drone_ID = 'mesange'
    exp_ID = '12-FRAC_001'
    flip_field = 1
    
    process_file(base, date, drone_ID, exp_ID, flip_field = flip_field)

if __name__ == '__main__':
    main()

