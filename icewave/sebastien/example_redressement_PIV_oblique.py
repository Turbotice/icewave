# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:29:21 2024

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt 
import glob
from scipy.interpolate import RegularGridInterpolator
import time 

import h5py

import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp
import icewave.tools.matlab2python as mat2py

   
#%% 
date = '0211'
year = '2024'
drone_ID = 'bernache'
flight_ID = 'structured_data'

base = df.find_path('Hublot24')
base = 'H:/Rimouski_2024/Data/'
path2data = f'{base}{year}/{date}/Drones/{drone_ID}/matData/{flight_ID}/'
filelist_mat = glob.glob(f'{path2data}*.mat')

structured_mat = glob.glob(f'{path2data}*scaled.mat')[0]


#%%
# Load the .mat file and convert to a dictionary

with h5py.File(structured_mat, 'r') as fmat:
    mat_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    mat_dict = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    mat_dict['PIV_param']['p_param'] = mat2py.matcell2dict_PIV(mat_dict['PIV_param']['p_param'])
    mat_dict['PIV_param']['s_param'] = mat2py.matcell2dict_PIV(mat_dict['PIV_param']['s_param'])

#%% Convert Vy displacement field to vertical field Vz 
key_pix = 'PIXEL'

Vy = mat_dict['Vy'] # time, y, x
Vy_transposed = np.transpose(Vy,(2,1,0)) 

fps = 1/mat_dict['SCALE']['ft']
Dt = mat_dict['PIV_param']['Dt'] # time step between two image that are compared using PIV
t = mat_dict['t']*fps
x = mat_dict[key_pix]['x_pix'] # pixel position of each PIV box
y = mat_dict[key_pix]['y_pix'] # pixel position of each PIV box 
# compute interpolation of pixel displacement vertical velocity field 
Fy = RegularGridInterpolator((x,y,t), Vy_transposed)


#%%


Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed,x,y,t,mat_dict[key_pix]['y0'],mat_dict['DRONE']['h_drone'],
                                               mat_dict['DRONE']['alpha_0'],mat_dict['DRONE']['focale'],fps,Dt)
    



#%% Plot a plt.colormesh 
W = mat_dict['PIV_param']['w']
x_edges = np.resize(x,x.size + 1)
y_edges = np.resize(y,y.size + 1)
x_edges[-1] = x_edges[-2] + W
y_edges[-1] = y_edges[-2] + W
x_edges = x_edges - W/2
y_edges = y_edges - W/2


Xedges,Yedges = np.meshgrid(x_edges,y_edges, indexing = 'ij')

Xreal,Yreal = dp.projection_real_space(Xedges, Yedges, mat_dict[key_pix]['x0'], mat_dict[key_pix]['y0'], mat_dict['DRONE']['h_drone']
                                       ,mat_dict['DRONE']['alpha_0'],mat_dict['DRONE']['focale'])

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,Vz[:,:,1500],shading = 'auto')
fig.colorbar(c,ax = ax)





