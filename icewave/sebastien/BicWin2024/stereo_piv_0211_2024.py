# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:57:43 2025

@author: sebas

This code aims at computing the real velocity field (ux,uy,uz) from the apparent velocity field (Vx,Vy)^(j)
measured by two drones

"""


import numpy as np
import os
import glob
import pandas as pd 
import h5py
import pickle
import cv2 as cv

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.path import Path

import scipy.signal as signal
from scipy.interpolate import LinearNDInterpolator, griddata
from datetime import datetime,timedelta
import pytz

from shapely.geometry import Polygon, Point
from scipy.spatial import ConvexHull

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.tools.rw_data as rw 

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
#%% FUNCTION SECTION

def transpose_PIVmat_fields(m):
    """ Change dimensions of different fields computed using PIVlab """
    key_fields = ['Vx','Vy','Vz','X','Y']

    for key in key_fields:
        m[key] = np.transpose(m[key])
        
    return m

#--------------------------------------------------------------------------------------

def synchronize_with_all_drones(m,istart = 0, iend = -1):
    """ Keep only data that is common to all drones """
    
    fields_with_time = ['Vx','Vy','Vz']
    for key in fields_with_time:
        m[key] = m[key][:,:,istart:iend]
    
    m['t'] = m['t'][istart:iend]
    
    return m


#---------------------------------------------------------------------------------------

def points_to_polygon(points):
    """Create a convex hull polygon from a set of points."""
    points = np.array(points)
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    return Polygon(hull_points)

#---------------------------------------------------------------------------------------

def interpolate_field(known_points,field,grid_x,grid_y):
    """ Interpolate a 2D field, known on a set of coordinates known_points, onto coordinates grid_x,grid_y 
    Inputs : """

    gridded_field = griddata(known_points, field.ravel(), 
                              (grid_x,grid_y),method = 'linear')
    
    return gridded_field

#---------------------------------------------------------------------------------------


def get_zstar(H,alpha,Y):
    """ Compute distance z_star in the projection geometry (distance between the camera and the plane 
    orthogonal to optical axis and crossing the observed point on the ice, which coordinate is Y in the reference 
    frame work of the drone) 
    
    Inputs : - H, drone height in meter
             - alpha, angle with respect to horizontal (in rad)
             - Y array like or float """
             
    z_star = H/np.sin(alpha) + np.cos(alpha)*Y
    return z_star

#---------------------------------------------------------------------------------------------------------------------------
# comute system coefficients 
def compute_coeffs(H,alpha,theta,focale,X,Y):
    """ Compute matrix coefficients for velocity field inversion
    Inputs : - H, float, drone height in meter
             - alpha, float, camera pitch angle with respect to the horizontal, in rad 
             - theta, float, yaw angle between the coordinate system (X,Y) that is projected and the coordinate system 
             of reference chosen to be projected. theta is taken clockwise and expressed in rad
             - focale, float, camera focal length
             - X,Y, coordinates X,Y in the coordinates system which is projected """

    if np.logical_or(alpha <0, alpha > np.pi/2):
        raise Exception("alpha does not belongs to [0,pi/2]")
    
    z_star = get_zstar(H, alpha, Y)
    
    a = (np.cos(theta) + np.sin(theta)*np.cos(alpha)*X/z_star)* focale / z_star
    b = (np.sin(theta) -np.cos(theta)*np.cos(alpha)*X/z_star)* focale / z_star
    c = np.sin(alpha)*X*focale/z_star**2
    d = (-np.sin(theta)*np.sin(alpha) + np.sin(theta)*np.cos(alpha)*Y*np.sin(alpha)/z_star)* focale / z_star
    e = (np.cos(theta) * np.sin(alpha) - np.cos(theta) * np.cos(alpha)*np.sin(alpha)*Y/z_star)*focale/z_star
    f = (np.cos(alpha) + np.sin(alpha)*np.sin(alpha)*Y/z_star)*focale/z_star

    coeffs = np.array((a,b,c,d,e,f))
    return coeffs

#-----------------------------------------------------------------------------------------------------------------------------

def create_matrix(param_ref,param_proj,grid_x,grid_y,grid_x_proj,grid_y_proj):
    """ Create matrix 4 x 3 for each position. This matrix needs to be solved in least square method to compute the real field 
    components. 
    Inputs : - param_ref , dictionnary, drone parameter of reference drone
             - param_proj, dictionnary, drone parameter of projected drone
             - grid_x,grid_y, array like shape [ny,nx], X,Y meshgrid coordinates in the reference coordinate system
             - grid_x_proj,y_grid_proj, array like, shape [ny,nx], X_proj, Y_proj meshgrid coordinates. Correspond
             to grid_x,grid_y coordinates projected in the projected coordinate system. 
    Output : M, numpy array,shape [4,3,ny,nx], M[:,:,i,j] contains the matrix to be inverted at position grid_x[i,j],grid_y[i,j]
    in order to compute the field (ux,uy,uz) at this specific position. """

    theta_ref = 0
    theta_proj = (param_ref['azimuth'] - param_proj['azimuth'])*np.pi/180

    coeffs_ref = compute_coeffs(param_ref['h'],param_ref['alpha_0'],theta_ref,param_ref['focale'],
                             grid_x,grid_y)
    coeffs_proj = compute_coeffs(param_proj['h'],param_proj['alpha_0'],theta_proj,param_proj['focale'],
                             grid_x_proj,grid_y_proj)

    # define matrix
    M = np.array([
    [coeffs_ref[0,:,:], coeffs_ref[1,:,:], coeffs_ref[2,:,:]],
    [coeffs_ref[3,:,:], coeffs_ref[4,:,:], coeffs_ref[5,:,:]],
    [coeffs_proj[0,:,:], coeffs_proj[1,:,:], coeffs_proj[2,:,:]],
    [coeffs_proj[3,:,:], coeffs_proj[4,:,:], coeffs_proj[5,:,:]]   
    ])

    return M

#-------------------------------------------------------------------------------------------------------------------------

def solve_linear_system_all_pos(M,b):
    """ Solve linear system for all positions, using 4 equations and least square method. 
    Inputs : - M, numpy array [4,3,ny,nx], matrix containing sub-matrices (4x3) for each position
             - b, numpy array [4,ny,nx] containing apparent velocity field observed on each camera for each position 
             of the studied grid. 
             - mask, array like, boolean array, contains 0 if the corresponding position does not belong to the studied grid
    Output : 
             - solutions, numpy array, [3,ny,nx], contain real velocity field (ux,uy,uz) for all positions of the studied grid
    """

    solutions = np.zeros((3,b.shape[1],b.shape[2]))

    for i in range(M.shape[2]):
        for j in range(M.shape[3]):
            # if ~mask[i,j] :
            #     solutions[:,i,j] = None

            sol = np.linalg.lstsq(M[:,:,i,j],b[:,i,j],rcond = None)
                
            solutions[:,i,j] = sol[0]
    return solutions


#%% Load data

date = '0211'
path2data = f'E:/Rimouski_2024/Data/{date}/Drones/'
drones = ['bernache','mesange']
data = {}
for i in range(len(drones)):
    filelist = glob.glob(f'{path2data}*{drones[i]}/matData/*scaled.mat')    
    file2load = filelist[0]
    with h5py.File(file2load, 'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
    
        data[drones[i]] = mat2py.mat_to_dict(fmat['m'],fmat['m'])
        
#%% Load synchro file 

filename = f'{path2data}synchro_3drones.pkl'

with open(filename,'rb') as pf:
   synchro = pickle.load(pf)
   
# Synchronize data and structure them
drone_keys = ['bernache','mesange']
for key in drone_keys :
    data[key] = transpose_PIVmat_fields(data[key])
    istart = synchro[key]
    data[key] = synchronize_with_all_drones(data[key],istart = istart)
    
#%% Load all procruste operations 

procruste_file = f'{path2data}procruste_operation_i0_520_Nbframes_4000_any_drones.pkl'

with open(procruste_file,'rb') as pf:
    all_op = pickle.load(pf)

# keep the given set of procruste operations 
ref_drone = 'mesange' # drone chosen as reference
projected_drone = 'bernache'

# get drone GPS position and orientation 
file2parameters = f'{path2data}parameters_3drones.pkl'
with open(file2parameters,'rb') as pf :
    param = pickle.load(pf)
    
#%% Define grid in the coordinate system of mesange 

minx = -58
maxx = 48
miny = -25
maxy = 50

artificial_facq_x = 1/0.8 # artificial spatial frequency in box/meter
grid_x, grid_y = np.meshgrid(
    np.linspace(minx, maxx, int((maxx - minx)*artificial_facq_x)),  
    np.linspace(miny, maxy, int((maxy - miny)*artificial_facq_x))
)

# Mask grid points outside intersection
# mask = np.array([Point(x, y).within(intersection) 
#                  for x, y in zip(grid_x.ravel(), grid_y.ravel())])
# mask = mask.reshape(grid_x.shape)

key_procruste = f'{projected_drone}_2_{ref_drone}'
procruste_op = all_op[key_procruste]

X_ref = data[ref_drone]['X']
Y_ref = data[ref_drone]['Y']
ref_points = np.array([X_ref.ravel(),Y_ref.ravel()]).T

X = data[projected_drone]['X']
Y = data[projected_drone]['Y']

# initialize final velocity field
Nb_frames = data[ref_drone]['Vx'].shape[2]
u = np.zeros((3,grid_x.shape[0],grid_x.shape[1],Nb_frames))

for frame in range(Nb_frames):
    #### Interpolate reference velocity fields 
    gridded_Vx_ref = interpolate_field(ref_points,data[ref_drone]['Vx'][:,:,frame],grid_x,grid_y) 
    gridded_Vy_ref = interpolate_field(ref_points,data[ref_drone]['Vy'][:,:,frame],grid_x,grid_y)
    
    # get procruste operations
    R = procruste_op['rot'][:,:,frame]
    translation = procruste_op['translat'][:,frame]
    scaling = procruste_op['scaling'][frame]
    
    # compute X,Y in the new reference frame 
    X_proj,Y_proj = dp.change_XY_reference_system(X, Y, param[projected_drone], param[ref_drone], R, translation, scaling)
    projected_points = np.array([X_proj.ravel(),Y_proj.ravel()]).T
    
    #### Interpolate projected velocity fields 
    gridded_Vx_proj = interpolate_field(projected_points,data[projected_drone]['Vx'][:,:,frame],grid_x,grid_y) 
    gridded_Vy_proj = interpolate_field(projected_points,data[projected_drone]['Vy'][:,:,frame],grid_x,grid_y) 
    
    #### Compute grid coordinates in coordinate system of projected drone
    key_procruste_back = f'{ref_drone}_2_{projected_drone}'
    procruste_op_back = all_op[key_procruste_back]
    
    # get procruste operations
    R = procruste_op_back['rot'][:,:,frame]
    translation = procruste_op_back['translat'][:,frame]
    scaling = procruste_op_back['scaling'][frame]
    
    # compute grid coordinates back to projected drone frame work
    grid_proj_back_X,grid_proj_back_Y = dp.change_XY_reference_system(grid_x, grid_y, 
                                                            param[ref_drone], param[projected_drone], 
                                                            R, translation, scaling)
    
    # create matrix to invert for all positions
    M = create_matrix(param[ref_drone],param[projected_drone],grid_x,grid_y,grid_proj_back_X,grid_proj_back_Y)
    
    # member b 
    b = np.array([
        gridded_Vx_ref,
        gridded_Vy_ref,
        gridded_Vx_proj,
        gridded_Vy_proj
    ])
    
    # solve linear system using least square method 
    solutions = solve_linear_system_all_pos(M,b)
    
    u[:,:,:,frame] = solutions
    print(f'Frame {frame} processed !')
 
#%% Create a structure to be saved

u_struct = {}
u_struct['u'] = u
u_struct['ref_drone'] = ref_drone
u_struct['projected_drone'] = projected_drone
u_struct['X'] = grid_x
u_struct['Y'] = grid_y
u_struct['extents'] = np.array([minx,maxx,miny,maxy])
u_struct['t'] = data[ref_drone]['t']
u_struct['units'] = {'u':'m/s','X':'m','Y':'m','extents':'m','t':'s'}
u_struct['u_shape'] = {'dim_0':'u_componants','dim_1':'ny','dim_2':'nx','dim_3':'nt'}

# scales
u_struct['SCALE'] = {}
u_struct['SCALE']['facq_t'] = data[ref_drone]['SCALE']['facq_t']
u_struct['SCALE']['facq_x'] = artificial_facq_x
u_struct['SCALE']['units'] = {'facq_t':'Hz','facq_x':'box/meter'}
 
# drone parameters
u_struct['param'] = {'ref_drone':param[ref_drone],'projected_drone':param[projected_drone]}
u_struct['synchro'] = synchro


file2save = f'{path2data}real_field_stereo_{date}_2024_rectangular_grid.h5'
rw.save_dict_to_h5(u_struct, file2save)






















