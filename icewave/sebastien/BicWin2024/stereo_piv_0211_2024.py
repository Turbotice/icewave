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


def get_zstar(H,alpha,Y):
    """ Compute distance z_star in the projection geometry (distance between the camera and the plane 
    orthogonal to optical axis and crossing the observed point on the ice, which coordinate is Y in the reference 
    frame work of the drone) 
    
    Inputs : - H, drone height in meter
             - alpha, angle with respect to horizontal (in rad)
             - Y array like or float """
             
    z_star = H/np.sin(alpha) + np.cos(alpha)*Y
    return z_star

# comute system coefficients 

def compute_coeffs(H,alpha,focale,X,Y):
    """ Compute matrix coefficients """
    
    z_star = get_zstar(H, alpha, Y)
    
    a = focale / z_star
    b = -np.cos(alpha)*X*focale/z_star**2
    c = np.sin(alpha)*X*focale/z_star**2
    d = np.zeros(np.shape(X))
    e = (np.sin(alpha) - np.cos(alpha)*np.sin(alpha)*Y/z_star)*focale/z_star
    f = (np.cos(alpha) + np.sin(alpha)*np.sin(alpha)*Y/z_star)*focale/z_star

    coeffs = np.array((a,b,c,d,e,f))
    return coeffs

#%% Load PIV data from both drones 

path2data = 'Y:/Banquise/Baptiste/Resultats/Analyse_1102/Data/'
drones = ['bernache','mesange']
data = {}
for i in range(len(drones)):
    filelist = glob.glob(f'{path2data}*{drones[i]}.mat')    
    file2load = filelist[0]
    with h5py.File(file2load, 'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
    
        data[drones[i]] = mat2py.mat_to_dict(fmat['m'],fmat['m'])
        
#%% Load synchro file 

path_save = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path_save}synchro_3drones.pkl'

with open(filename,'rb') as pf:
   synchro = pickle.load(pf)
   
# Synchronize data and structure them
drone_keys = ['bernache','mesange']
for key in drone_keys :
    data[key] = transpose_PIVmat_fields(data[key])
    istart = synchro[key]
    data[key] = synchronize_with_all_drones(data[key],istart = istart)
    
#%% Load all procruste operations 
path2procruste = 'W:/SagWin2024/Data/0211/Drones/'
procruste_file = f'{path2procruste}procruste_operation_i0_520_Nbframes_4000_any_drones.pkl'

with open(procruste_file,'rb') as pf:
    all_op = pickle.load(pf)

# keep the given set of procruste operations 
ref_drone = 'mesange' # drone chosen as reference
projected_drone = 'bernache'


# get drone GPS position and orientation 
path2parameters = 'W:/SagWin2024/Data/0211/Drones/'
file2parameters = f'{path2parameters}parameters_3drones.pkl'
with open(file2parameters,'rb') as pf :
    param = pickle.load(pf)
    
    
#%% Apply procruste_operations to bernache vertical wave field for a given frame 

Vz = data[projected_drone]['Vz']
key_procruste = f'{projected_drone}_2_{ref_drone}'
procruste_op = all_op[key_procruste]

# select frame 
frame = 500

# get procruste operations
R = procruste_op['rot'][:,:,frame]
translation = procruste_op['translat'][:,frame]
scaling = procruste_op['scaling'][frame]

# get X,Y positions of the whole field
X = data[projected_drone]['X']
Y = data[projected_drone]['Y']

# compute X,Y in the new reference frame 
X_proj,Y_proj = dp.change_XY_reference_system(X, Y, param[projected_drone], param[ref_drone], R, translation, scaling)


#%% Plot wave field 
Vz_ref = data[ref_drone]['Vz']
X_ref = data[ref_drone]['X']
Y_ref = data[ref_drone]['Y']

fig, ax = plt.subplots()
c = ax.pcolormesh(X_ref,Y_ref,Vz_ref[:,:,frame],shading = 'gouraud', cmap = 'viridis',vmin = -4, vmax = 4)
c.set_rasterized(True)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_title(f'{ref_drone} in frame work of {ref_drone}')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])
fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}Vz_i{frame}_{ref_drone}_framework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


fig, ax = plt.subplots()
c = ax.pcolormesh(X_proj,Y_proj,Vz[:,:,frame],shading = 'gouraud', cmap = 'viridis',vmin = -4, vmax = 4)
c.set_rasterized(True)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_title(f'{projected_drone} in frame work of {ref_drone}')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])
fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}Vz_i{frame}_{projected_drone}_framework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Superpose velocity fields from both views 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.pcolormesh(X_ref,Y_ref,Vz_ref[:,:,frame],shading = 'gouraud', cmap = parula_map,vmin = -4, vmax = 4,
                  alpha = 0.6)
c.set_rasterized(True)

cproj = ax.pcolormesh(X_proj,Y_proj,Vz[:,:,frame],shading = 'gouraud', cmap = parula_map,vmin = -4, vmax = 4,
                      alpha = 0.3)
cproj.set_rasterized(True)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])
fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}superposition_Vz_i{frame}_{projected_drone}_inframework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Keep only coordinates common to both drones 

top_left = [np.min(X_ref),np.max(Y_ref)]
top_right = [np.max(X_ref),np.max(Y_ref)]

# compute bottom left and bottom right vertices coordinates
Y_min = np.min(Y)

mask = Y == Y_min
bottom_left = [np.min(X_ref[mask]),Y_min]
bottom_right = [np.max(X_ref[mask]),Y_min]

vertices_ref = np.array((bottom_left,top_left,top_right,bottom_right))

# keep only coordinates inside this polygon 
projected_points = []
for i in range(X_proj.shape[0]):
    for j in range(X_proj.shape[1]):
        projected_points.append((X_proj[i,j],Y_proj[i,j]))
        
projected_points = np.array(projected_points)
polygon = Path(vertices_ref)
inside = polygon.contains_points(projected_points)

inside = inside.reshape(X_proj.shape)

#%% Get coordinates of common system 

projected_points = []
for i in range(X_proj.shape[0]):
    for j in range(X_proj.shape[1]):
        projected_points.append((X_proj[i,j],Y_proj[i,j]))
        
ref_points = []
for i in range(X_ref.shape[0]):
    for j in range(X_ref.shape[1]):
        ref_points.append((X_ref[i,j],Y_ref[i,j]))

def points_to_polygon(points):
    """Create a convex hull polygon from a set of points."""
    points = np.array(points)
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    return Polygon(hull_points)

# Build polygons from convex hulls
poly_ref = points_to_polygon(ref_points)
poly_proj = points_to_polygon(projected_points)

# Compute intersection polygon
intersection = poly_ref.intersection(poly_proj)

# Keep only points inside intersection
points_in_common_ref = [pt for pt in ref_points if Point(pt).within(intersection)]
points_in_common_proj = [pt for pt in projected_points if Point(pt).within(intersection)]

points_in_common_ref = np.array(points_in_common_ref)
points_in_common_proj = np.array(points_in_common_proj)


#%% Plot field and common area

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.pcolormesh(X_ref,Y_ref,Vz_ref[:,:,i],shading = 'gouraud', cmap = parula_map,vmin = -4, vmax = 4,
                  alpha = 0.6)
c.set_rasterized(True)

cproj = ax.pcolormesh(X_proj,Y_proj,Vz[:,:,i],shading = 'gouraud', cmap = parula_map,vmin = -4, vmax = 4,
                      alpha = 0.3)
cproj.set_rasterized(True)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

ax.plot(points_in_common_ref[:,0],points_in_common_ref[:,1],'.',color = 'tab:blue')
ax.plot(points_in_common_proj[:,0],points_in_common_proj[:,1],'.',color = 'tab:red')

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])



########################################################################################################################
############################# Compute velocity fields on a regular grid ################################################


#%% Create a common coordinates system 

# ---------------------------
# Build grid over intersection
# ---------------------------

artificial_facq_x = 1/0.8 # artificial spatial frequency in box/meter
minx, miny, maxx, maxy = intersection.bounds
grid_x, grid_y = np.meshgrid(
    np.linspace(minx, maxx, int((maxx - minx)*artificial_facq_x)),   # resolution: increase/decrease 100
    np.linspace(miny, maxy, int((maxy - miny)*artificial_facq_x))
)

# Mask grid points outside intersection
mask = np.array([Point(x, y).within(intersection) 
                 for x, y in zip(grid_x.ravel(), grid_y.ravel())])
mask = mask.reshape(grid_x.shape)

# Interpolate set points of reference 

gridded_Vx_ref = griddata(np.array(ref_points), data[ref_drone]['Vx'][:,:,frame].ravel(), 
                          (grid_x,grid_y),method = 'linear')
gridded_Vy_ref = griddata(np.array(ref_points), data[ref_drone]['Vy'][:,:,frame].ravel(), 
                          (grid_x,grid_y),method = 'linear')

# keep only points within the common area
gridded_Vx_ref[~mask] = np.nan
gridded_Vy_ref[~mask] = np.nan

gridded_Vx_proj = griddata(np.array(projected_points),data[projected_drone]['Vx'][:,:,frame].ravel(),
                           (grid_x,grid_y),method = 'linear')
gridded_Vy_proj = griddata(np.array(projected_points), data[projected_drone]['Vy'][:,:,frame].ravel(), 
                          (grid_x,grid_y),method = 'linear')

gridded_Vx_proj[~mask] = np.nan
gridded_Vy_proj[~mask] = np.nan


#%% Plot the interpolated fields 

fig, axs = plt.subplots(nrows = 1,ncols = 2)
axs[0].imshow(gridded_Vx_ref,extent = (minx,maxx,miny,maxy))
axs[0].set_title('$V_x^{\mathrm{ref}}$')
axs[1].imshow(gridded_Vx_proj,extent = (minx,maxx,miny,maxy))
axs[1].set_title('$V_x^{\mathrm{proj}}$')

fig, axs = plt.subplots(nrows = 1,ncols = 2)
axs[0].imshow(gridded_Vy_ref,extent = (minx,maxx,miny,maxy))
axs[0].set_title('$V_y^{\mathrm{ref}}$')
axs[1].imshow(gridded_Vy_proj,extent = (minx,maxx,miny,maxy))
axs[1].set_title('$V_y^{\mathrm{proj}}$')

#%% Compute the coordinates of the new grid in the frame work of projected drone 

key_procruste = f'{ref_drone}_2_{projected_drone}'
print(key_procruste)
procruste_op = all_op[key_procruste]

# get procruste operations
R = procruste_op['rot'][:,:,frame]
translation = procruste_op['translat'][:,frame]
scaling = procruste_op['scaling'][frame]

# compute X,Y in the new reference frame 
grid_proj_back_X,grid_proj_back_Y = dp.change_XY_reference_system(grid_x, grid_y, 
                                                        param[ref_drone], param[projected_drone], 
                                                        R, translation, scaling)
#%% Superpose initial field and interpolated field 

fig, ax = plt.subplots()
c = ax.pcolormesh(grid_x,grid_y,gridded_Vx_ref, cmap = parula_map, alpha = 0.6)
c.set_rasterized(True)

c = ax.pcolormesh(data[ref_drone]['X'],data[ref_drone]['Y'],data[ref_drone]['Vx'][:,:,frame],
              shading = 'gouraud',cmap = parula_map,alpha = 0.3)
c.set_rasterized(True)
ax.set_title(f'{ref_drone}')

# bernache 
fig, ax = plt.subplots()
c = ax.pcolormesh(grid_proj_back_X,grid_proj_back_Y,gridded_Vx_proj, cmap = parula_map, alpha = 0.6)
c.set_rasterized(True)

c = ax.pcolormesh(data[projected_drone]['X'],data[projected_drone]['Y'],data[projected_drone]['Vx'][:,:,frame],
              shading = 'gouraud',cmap = parula_map,alpha = 0.3)
c.set_rasterized(True)

ax.set_title(f'{projected_drone}')




#%% Compute matrix coefficients 

#--------------------------
# ref_drone : mesange
#--------------------------

coeffs_ref = compute_coeffs(param[ref_drone]['h'],param[ref_drone]['alpha_0'],param[ref_drone]['focale'],
                             grid_x,grid_y)

#--------------------------
# projected_drone : bernache
#--------------------------

coeffs_proj = compute_coeffs(param[projected_drone]['h'],param[projected_drone]['alpha_0'],param[projected_drone]['focale'],
                             grid_proj_back_X,grid_proj_back_Y)

#%% Create matrix 

# M contains matrix 3x3 for all points of the grid 
M = np.array([
    [(coeffs_ref[0,:,:] + coeffs_proj[0,:,:])*0.5, (coeffs_ref[1,:,:] + coeffs_proj[1,:,:])*0.5, 
     (coeffs_ref[2,:,:] + coeffs_proj[2,:,:])*0.5 ], # coeffs of equation ((1) + (3))/2
    [coeffs_ref[3,:,:], coeffs_ref[4,:,:], coeffs_ref[5,:,:]], # coeffs of equation (2)
    [coeffs_proj[3,:,:], coeffs_proj[4,:,:],coeffs_proj[5,:,:]] # coeffs of equation (3)
    ])


# member b 
b = np.array([
    0.5*(gridded_Vx_ref + gridded_Vx_proj),
    gridded_Vy_ref,
    gridded_Vy_proj,
    ])

#%%
# solve linear system for all positions 

solutions = np.zeros(b.shape)

for i in range(M.shape[2]):
    for j in range(M.shape[3]):
        if ~mask[i,j] :
            solutions[:,i,j] = None

        else :
            sol = np.linalg.lstsq(M[:,:,i,j],b[:,i,j])
            
            solutions[:,i,j] = sol[0]

#%% Plot interpolated velocity fields 

# Plot ux
fig, ax = plt.subplots()
imsh = ax.imshow(solutions[0,:,:],extent = (minx,maxx,miny,maxy),origin = 'lower', cmap = parula_map)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_x$')

imsh.set_clim([-0.25,0.25])

# Plot uy
fig, ax = plt.subplots()
imsh = ax.imshow(solutions[1,:,:],extent = (minx,maxx,miny,maxy),origin = 'lower', cmap = parula_map)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_y$')

imsh.set_clim([-1,1])

# Plot uz
fig, ax = plt.subplots()
imsh = ax.imshow(solutions[2,:,:],extent = (minx,maxx,miny,maxy),origin = 'lower', cmap = parula_map)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$u_z$')

imsh.set_clim([-1,1])



#%% Compute coefficients for a given drone





fig, ax = plt.subplots()
c = ax.pcolormesh(X_ref,Y_ref,z_star,shading = 'gouraud', cmap = parula_map)
c.set_rasterized(True)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$z^* \; \mathrm{(m)}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')

# compute coeffs for a reference drone


# compute coeffs for projected_drone 

# build matrix 
    

########################################################################################################################
########################################################################################################################



#%% Interpolate the projected drone coordinates on the reference drone coordinates system 

mask_ref = np.array([Point(pt).within(intersection) for pt in ref_points])
mask_ref = np.reshape(mask_ref,X_ref.shape)

mask_proj = np.array([Point(pt).within(intersection) for pt in projected_points])
mask_proj = np.reshape(mask_proj,X_proj.shape)

interp_X = LinearNDInterpolator(points_in_common_proj, X_proj[mask_proj])
interp_Y = LinearNDInterpolator(points_in_common_proj, Y_proj[mask_proj])

new_X_proj = interp_X(X_ref[mask_ref],Y_ref[mask_ref])
new_Y_proj = interp_Y(X_ref[mask_ref],Y_ref[mask_ref])

#%% Plot new interpolated coordinates 

fig, ax = plt.subplots()

# ax.plot(X_ref[mask_ref],Y_ref[mask_ref],'d',color = 'tab:blue')
ax.plot(points_in_common_ref[:,0],points_in_common_ref[:,1],'d',color = 'tab:blue')
ax.plot(new_X_proj,new_Y_proj,'.',color = 'tab:red')


#%% Transform new interpolated coordinates (in mesange framework) back to X,Y coordinates of bernache 

key_procruste = f'{ref_drone}_2_{projected_drone}'
procruste_op = all_op[key_procruste]

# get procruste operations
R = procruste_op['rot'][:,:,frame]
translation = procruste_op['translat'][:,frame]
scaling = procruste_op['scaling'][frame]

# compute X,Y in the new reference frame 
X_proj_back,Y_proj_back = dp.change_XY_reference_system(new_X_proj, new_Y_proj, 
                                                        param[ref_drone], param[projected_drone], 
                                                        R, translation, scaling)


