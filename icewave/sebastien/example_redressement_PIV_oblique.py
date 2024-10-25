# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:29:21 2024

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 
import glob
from scipy.interpolate import RegularGridInterpolator

import h5py

import icewave.sebastien.seb as seb 
import icewave.tools.datafolders as df

#%%
def mat_to_dict(mat_object):
    """
    Recursively convert a MATLAB structure (HDF5 group or dataset) to a Python dictionary.
    
    """
    if isinstance(mat_object, h5py.Dataset):  # If it's a dataset, return its value
        data = mat_object[()]
                  
        # Handle MATLAB strings (stored as bytes)
        if data.dtype == 'uint16':  # Check if it's a string
            # Convert uint16 array to Python string (decode as Unicode characters)
            return ''.join(chr(code_point[0]) for code_point in data)
        return data
    
        if isinstance(data, np.ndarray):
            if data.size == 1:  # If the array contains only one element, convert to float
                    return float(data)
    
    elif isinstance(mat_object, h5py.Group):  # If it's a group (structure), create a dictionary
        result_dict = {}
        for key, item in mat_object.items():
            result_dict[key] = mat_to_dict(item)  # Recursively call for each element
        return result_dict
    
    else:
        raise TypeError(f"Unsupported type {type(mat_object)}")

def backward_projection(fun,x,y,t,y0,h,alpha_0,focale,fps,Dt):

    
    Fp = h*focale*fun((x,y,t))/((y - y0)*np.cos(alpha_0) + 
                              focale*np.sin(alpha_0))/((y - y0 + fun((x,y,t)))*np.sin(alpha_0) - focale*np.cos(alpha_0))*fps/Dt
    
    return Fp 


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
    # Iterate over all top-level keys in the .mat file
    for key in fmat['m'].keys():
        print(key)
        mat_dict[key] = mat_to_dict(fmat['m'][key])

#%% Convert Vy displacement field to vertical field Vz 
key_pix = 'PIXEL'

Vy = mat_dict['Vy'] # time, y, x
Vy_transposed = np.transpose(Vy,(2,1,0)) 

fps = 1/mat_dict['SCALE']['ft']
Dt = 4 # time step between two image that are compared using PIV
t = mat_dict['t'][:,0]*fps
x = mat_dict[key_pix]['x_pix'][:,0]
y = mat_dict[key_pix]['y_pix'][:,0]
# compute interpolation of pixel displacement vertical velocity field 
Fy = RegularGridInterpolator((x,y,t), Vy_transposed)

# compute interpolator of vertical velocity field F
Fz = backward_projection(Fy,x,y,t,mat_dict[key_pix]['y0'],mat_dict['DRONE']['h_drone'],mat_dict['DRONE']['alpha_0'],
                         mat_dict['DRONE']['focale'],fps,Dt)

#%%
t = t - t[0]
x = x - x[0]
y = y- y[0]
# compute interpolation of pixel displacement vertical velocity field 
Fy = RegularGridInterpolator((x,y,t), Vy_transposed)

a = np.array([1,2,3,4])
b = np.array([5,6,7])
c = np.array([1,2,3,4,5])

# Create the meshgrid of new points
new_x, new_y, new_t = np.meshgrid(a, b, c, indexing='ij')

# Combine the new coordinates into a single list of points to interpolate
new_points = np.array([new_x.ravel(), new_y.ravel(), new_t.ravel()]).T

# Perform the interpolation
interpolated_field = Fy(new_points)

# Reshape the interpolated field back to the new grid shape
interpolated_field = interpolated_field.reshape(len(a), len(b), len(c))

#%%
Fy((a,b,c))



