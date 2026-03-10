import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import glob
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import time 

import h5py
import os

#kerdir = os.getcwd()
#os.chdir('/home/vzanchi/Bureau/Turbotice_git/')

#sys.path.append('C:/Users/Vasco/OneDrive - Université de Paris/Documents/git/icewave')
#os.chdir(kerdir)

system_loc = 'linux_server'

#kerdir = os.getcwd()
#if ordi=='Leyre':
#    os.chdir('/run/user/1003/gvfs/smb-share:server=thiou.local,share=homes/vasco/Grenoble_Nov2024/post_traitement/piv_grenoble/fracture/python_functions/')
if system_loc=='linux_server':
    sys.path.append('/media/turbots/DATA/thiou/homes/vasco/Grenoble_Nov2024/post_traitement/piv_grenoble/fracture/python_functions/')
    sys.path.append('/media/vasco/OS/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')
elif system_loc=='windows_server':
    sys.path.append('Z:/vasco/Grenoble_Nov2024/post_traitement/piv_grenoble/fracture/python_functions/')
    sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')
print(os.listdir())


import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp

from organize import *



#if ordi=='Leyre':
#    base = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'


def reconstruct_Vz_field_from_pivdata(input_params = dict):

    date = input_params['date'] 
    acq_num = input_params['acq_num']
    f_exc = input_params['f_exc']
    freq_acq = input_params['freq_acq']
    cam_SN = input_params['cam_SN']
    delta_x = input_params['delta_x']
    delta_y = input_params['delta_y']
    imwidth = input_params['imwidth']
    imheight = input_params['imheight']
    x0 = input_params['x0']
    y0 = input_params['y0']
    pixel_size = input_params['pixel_size']
    focale_mm = input_params['focale_mm']
    focale = input_params['focale']
    h_camera = input_params['h_camera']
    alpha_0_deg = input_params['alpha_0_deg']
    alpha_0 = input_params['alpha_0']
    unit_input_velocity = input_params['unit_input_velocity']
    box_position_information = input_params['box_position_information']
    a = input_params['a']
    W = input_params['W']
    Dt = input_params['Dt']
    i0 = input_params['i0']
    N = input_params['N']

    ########################################################
    ############# importation des donnees piv ##############
    ########################################################
    #if system_loc=='linux_server':
    #    base =f'R:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    #elif system_loc=='windows_server':
    #    base =f'R:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    if system_loc=='linux_server':
         base =f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    elif system_loc=='windows_server':
        base =f'R:/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    
    path2data = f'{base}matData/'
    matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'

    with h5py.File(matfile, 'r') as fmat:
        mat_dict = {}
        
        print('Top-level keys : ', list(fmat.keys()))

        mat_dict = mat_to_dict(fmat['m'],fmat['m'])

    print("KEYS : ",mat_dict.keys())
    ########################################################

    # Convert Vy displacement field to vertical field Vz 
    Vy = mat_dict['Vy'] # time, y, x
    Vy_transposed = np.transpose(Vy,(2,1,0)) 
    
    if unit_input_velocity == 'pxpersec':
        fps = freq_acq
        t = mat_dict['t']
    elif unit_input_velocity == 'pxperframe':
        t = np.arange(len(mat_dict['t']))/fps  # temps en secondes
    
    if box_position_information==False:
        x = np.arange(len(mat_dict['x']))*(W/2) + W + 0.5 # pixel position of each PIV box
        y = np.arange(len(mat_dict['y']))*(W/2) + W + 0.5 # pixel position of each PIV box 
    elif box_position_information==True:
        #x = mat_dict['xpix']
        #y = mat_dict['ypix']
        x = mat_dict['x'] * (W/2)
        y = np.flip(mat_dict['y']) * (W/2)
        
        if a!=0:
            x = x[a:-a]
            y = y[a:-a]




    # compute interpolation of pixel displacement vertical velocity field 
    Fy = RegularGridInterpolator((x,y,t), Vy_transposed)

    # compute vertical velocity Vz
    if unit_input_velocity=='pxpersec':
        Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed/fps,x,y,t,y0,h_camera,
                                                alpha_0,focale,fps,Dt)
    elif unit_input_velocity=='pxperframe':
        Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed,x,y,t,y0,h_camera,
                                                alpha_0,focale,fps,Dt)

    # Plot a plt.colormesh 
    #W = mat_dict['PIV_param']['w']
    x_edges = np.resize(x,x.size + 1)
    y_edges = np.resize(y,y.size + 1)
    x_edges[-1] = x_edges[-2] + W
    y_edges[-1] = y_edges[-2] + W
    x_edges = x_edges - W/2
    y_edges = y_edges - W/2


    Xedges,Yedges = np.meshgrid(x_edges,y_edges, indexing = 'ij')

    Xreal,Yreal = dp.projection_real_space(Xedges, Yedges, x0, y0, h_camera
                                        ,alpha_0,focale)

    fig, ax = plt.subplots()
    pc = ax.pcolormesh(Xreal,Yreal,Vz[:,:,20],shading = 'auto',vmin=-0.01,vmax=0.01) # pcolormesh object !
    fig.colorbar(pc)
    #ax.set_ylim(-0.7,-0.2)
    #ax.set_xlim(-1.0,1.0)

    X,Y = np.meshgrid(x,y)

    Xreal_centers,Yreal_centers = dp.projection_real_space(X, Y, x0, y0, h_camera, alpha_0, focale)

    return mat_dict,t,Xreal,Yreal,Vz,pc,Xreal_centers,Yreal_centers