# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 15:18:30 2025

@author: sebas
"""


import numpy as np
import os
import pandas as pd 
import h5py
import pickle
import cv2 as cv

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import matplotlib.cm as cm

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp


#%% FUNCTION SECTION

def outliers2nans(S,keys,threshold = 100):
    """ Supress outliers from buoys detection 
    Inputs : - S, python dictionnary
             - keys, list of keys of S, ordered in a convenient way
             - theshold, max distance between two consecutive points to keep it
             
    Outputs : S_out, python dictionnary, same size and same keys as S. But outliers have been 
            replace by nan values"""
            
    set_graphs.set_matplotlib_param('single')

    fig,ax = plt.subplots()
    for i in range(2,len(keys) - 1,2):
        print(keys[i])
        
        label_buoy = f'Buoy {(i-1)//2 + 1}'
        ax.plot(S[keys[i]],S[keys[i+1]],'.',label = label_buoy)

        # check if nans are contained in each trajectory
        Nb_nans = np.sum(np.isnan(S[keys[i]]))
        print(f'Nb nans = {Nb_nans}')

    ax.legend()

    # Compute distance in pixels from previous detected point

    def distance (x,y,x0,y0) : 
        dist = np.sqrt((x - x0)**2 + (y - y0)**2)
        return dist
    
    fig,ax = plt.subplots()
    for i in range(2,len(keys) - 1,2):
        x = S[keys[i]]
        y = S[keys[i+1]]
        
        dist = distance(x[1:],y[1:],x[:-1],y[:-1])
        # set outliers to nan
        S[keys[i]][np.where(dist > threshold)[0]] = float('nan')
        S[keys[i+1]][np.where(dist > threshold)[0]] = float('nan')
        

        label_buoy = f'Buoy {(i-1)//2 + 1}'
        ax.plot(dist,'.',label = label_buoy)

    ax.legend()
    
    return S
    
    
def procrustes_alignement(points_project,points_ref,output_rescaling = 0):
    """ Align two set of corresponding points, using a rotation and a translation. 
    Inputs : - point_project, numpy array, set of points to be projected
            Points coordinates should be the last dimension
             - points_ref, numpy array, set of points to used as reference
    Outputs : - R, rotation matrix to apply to points_project
              - t, translation vector to apply to points_project """
              
    # Compute centers
    c_proj = np.nanmean(points_project,axis = 0)
    c_ref = np.nanmean(points_ref,axis = 0)
    # reduce center position
    centered_proj = points_project - c_proj
    centered_ref = points_ref - c_ref
    
    # compute rescaling factor if output_rescaling = 1
    if output_rescaling :
        scale_factor = np.sqrt(np.sum(centered_ref**2)/np.sum(centered_proj**2))
    else :
        scale_factor = 1
    
    # Compute best rotation
    H = centered_proj.T @ centered_ref # Cross-covariance matrix
    U,S,Vt = np.linalg.svd(H) # SVD
    R = Vt.T @ U.T # Compute optimal rotation

    # Ensure a proper rotation (no reflection)
    if np.linalg.det(R) < 0:
        Vt[-1,:] *= -1
        R = Vt.T @ U.T  
        
    # Compute best translation 
    t = c_ref - scale_factor * R @ c_proj
        
    return scale_factor,R,t 


def XY2GPS(X,Y,Lat0,Long0,azimuth):
    """ Convert cartesian coordinates X,Y to GPS coordinates
    Inputs : - X,Y : numpy array, containing cartesian (X,Y) coordinates with respect to a reference point O 
             - Lat0,Long0 : GPS coordinates of reference point (in degrees)
             - azimuth : orientation of the cartesian system with respect to the geographic north (in degrees, between 
                                                                                                   0 and 360)
             This azimuth angle is measured between geographic north direction and axis Y of cartesian system
    Outputs : - Lat,Long : numpy array, GPS coordinates of the selected points"""
    
    # convert cartesian to polar coordinates
    rho,theta = dp.cart2pol(X,Y)
    theta = theta*180/np.pi
    # compute azimuth of all selected points
    local_azimuth = azimuth + 90 - theta
    Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                    local_azimuth,rho)
    return Lat,Long

def GPS2XY(lat,long,Lat0,Long0,azimuth,R_earth = 6371e3):
    """Compute cartesian coordinates (X,Y) from GPS coordinates
    Inputs : - lat,long : numpy array, GPS coordinates (in degrees)
             - Lat0,Long0 : GPS coordinates of the point chosen as a reference for the cartesian coordinates system (in degrees)
             - azimuth : angle with respect to geopgraphic north with which the cartesian system will be oriented (in degrees,
                                                                                                            between 0 and 360)
             This azimuth angle is measured between geographic north direction and axis Y of cartesian system
    Outputs : - X,Y : numpy array, cartesian coordinates (in meter) """
    
    rho = R_earth*np.sqrt(((lat - Lat0)*np.pi/180)**2 + np.cos(Lat0*np.pi/180)**2 * ((long - Long0)*np.pi/180)**2)

    psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(long - Long0)*np.pi/180,(lat - Lat0)*np.pi/180)*180/np.pi

    theta = azimuth - psi + 90
    # compute euclidean coordinates
    X,Y = dp.pol2cart(rho,theta*np.pi/180)
    
    return X,Y

def replace_nans(arr):
    """ Replace nan values of a 1D array """   
    valid_idx = np.where(~np.isnan(arr))[0]
    valid_values = arr[valid_idx]

    nan_idx = np.where(np.isnan(arr))[0]
    arr[nan_idx] = np.interp(nan_idx,valid_idx,valid_values)
    
    return arr

#%% Load all detections

root = 'W:/SagWin2024/Data/0211/Drones/'
drones = ['fulmar','bernache','mesange']

keys = ['idx','t','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','x6','y6']
Z = {}

for drone_ID in drones:
    
    base = f'{root}{drone_ID}/buoys_tracking_0211_V2/'
    file2load = f'{base}Buoys_tracking_pix_{drone_ID}.mat'
    
    with h5py.File(file2load, 'r') as fmat:
        data_wave = {}
        
        print('Top-level keys : ', list(fmat.keys()))
    
        S = mat2py.mat_to_dict(fmat['S_buoys'],fmat['S_buoys'])

    # Supress outliers + restructure data
    
    S = outliers2nans(S,keys)
    
    buoys = np.zeros((6,2,len(S['idx'])))
    for i in range(2,len(keys) - 1,2):
        print(keys[i])
        buoys[(i-1)//2,0,:] = S[keys[i]]
        buoys[(i-1)//2,1,:] = S[keys[i+1]]
        
    Z[drone_ID] = {'idx':S['idx'],'t':S['t'],'buoys':buoys}
    
Z['bernache']['buoys'] = np.flip(Z['bernache']['buoys'],axis = 0)

#%% Find Nan values and interpolate

# correct buoy 5 position for fulmar 
mask = np.where(Z['fulmar']['buoys'][4,0,:] < 500)[0]
for i in range(2):
    Z['fulmar']['buoys'][4,i,:][mask] = float('nan')

# correct buoy 4 position for fulmar
mask = np.where(np.logical_and(Z['fulmar']['buoys'][3,0,:] < 1500,Z['fulmar']['buoys'][3,1,:] > 1180))[0]
for i in range(2):
    Z['fulmar']['buoys'][3,i,:][mask] = float('nan')

for drone_key in Z.keys():
    for i in range(5):
        for j in range(2):
            arr = Z[drone_key]['buoys'][i,j,:]
            Z[drone_key]['buoys'][i,j,:] = replace_nans(arr)
    
#%%
key_drone = 'fulmar'
  
fig, ax = plt.subplots()
for i in range(5):
    ax.plot(Z[key_drone]['buoys'][i,0,:],Z[key_drone]['buoys'][i,1,:],'.')


#%% Define folder where figures are saved

fig_folder = f'{root}Resultats_f2700/'
if not os.path.isdir(fig_folder) :
    os.mkdir(fig_folder)


#%% Load image from a drone 
path2img = f'{root}Initial_images/mesange_im_11496.tiff'

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

#%% Set drone parameters

file2save = f'{root}parameters_3drones.pkl'
# if os.path.isfile(file2save):
#     print(f'{file2save} already exists')
#     with open(file2save,'rb') as pf:
#         param = pickle.load(pf)
        
# else :
    # Parameter Fulmar
param = {'fulmar':{},'bernache':{},'mesange':{}}

theta = 32.75 # AFOV fulmar
param['fulmar']['h'] = 140 
param['fulmar']['alpha_0'] = np.pi/2
param['fulmar']['focale'] = 3840/2/np.tan(theta*np.pi/180)
# param['fulmar']['focale'] = 2800
param['fulmar']['x0'] = (3840+1)/2
param['fulmar']['y0'] = (2160+1)/2
param['fulmar']['latitude'] = 48.253089
param['fulmar']['longitude'] = -70.091383
param['fulmar']['azimuth'] = 268


 
# Parameter bernache
param['bernache']['h'] = 89.9 
param['bernache']['alpha_0'] = 60.1*np.pi/180
param['bernache']['focale'] = 2700
param['bernache']['x0'] = (3840+1)/2
param['bernache']['y0'] = (2160+1)/2
param['bernache']['latitude'] = 48.25307
param['bernache']['longitude'] = -70.09227
param['bernache']['azimuth'] = 85.8

# Parameter mesange
param['mesange']['h'] = 90.3 
param['mesange']['alpha_0'] = 59.8*np.pi/180
param['mesange']['focale'] = 2700
param['mesange']['x0'] = (3840+1)/2
param['mesange']['y0'] = (2160+1)/2
param['mesange']['latitude'] = 48.25317
param['mesange']['longitude'] = -70.09056
param['mesange']['azimuth'] = 268.1

with open(file2save,'wb') as pf :
    pickle.dump(param,pf)
    print(f'{file2save} saved !')


# Create a colormap for each drone 

maps = {'fulmar':'Blues','bernache':'Oranges','mesange':'Greens'}
maps_drone = {}
for key_drone in maps.keys():
    maps_drone[key_drone] = {}
    full_map = mpl.colormaps[maps[key_drone]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    cnorm = mcolors.Normalize(vmin = 0, vmax = 5)
    cnorm_time = mcolors.Normalize(vmin = 0, vmax = Z[key_drone]['t'].max())
    maps_drone[key_drone]['new_map'] = new_map
    maps_drone[key_drone]['cnorm'] = cnorm
    maps_drone[key_drone]['cnorm_time'] = cnorm_time

    
#%% Compute (X,Y) coordinates of each buoy

fig, ax = plt.subplots()

for drone_key in Z.keys():
    
    Z[drone_key]['real'] = np.zeros(np.shape(Z[drone_key]['buoys']))
    for i in range(6):
        X,Y = dp.projection_real_space(Z[drone_key]['buoys'][i,0,:],Z[drone_key]['buoys'][i,1,:],
                                       param[drone_key]['x0'],param[drone_key]['y0'],
                                       param[drone_key]['h'],param[drone_key]['alpha_0'],param[drone_key]['focale'])
        
        Z[drone_key]['real'][i,0,:] = X
        Z[drone_key]['real'][i,1,:] = Y
        
        current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](i))
        ax.plot(X,Y,'.',color = current_color)
        
ax.set_xlabel(r'$X$')
ax.set_ylabel(r'$Y$')

figname = f'{fig_folder}Raw_XY_buoys_3drones'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


#%% Compute buoys position in Latitude / Longitude coordinates

for drone_key in ['fulmar','bernache','mesange']:
    
    Z[drone_key]['GPS'] = np.zeros(Z[drone_key]['real'].shape)
    
    # horizontal distance between center of metric coordinates and drone position 
    if param[drone_key]['alpha_0'] != np.pi/2 :
        
        dist2drone = param[drone_key]['h']/np.tan(param[drone_key]['alpha_0'])
        Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[drone_key]['latitude'],param[drone_key]['longitude'],
                                                    param[drone_key]['azimuth'],dist2drone)
    else :
        Lat0 = param[key_drone]['latitude']
        Long0 = param[key_drone]['longitude']
        
    #GPS coordinates of detected objects 
    POS_rho,POS_theta = dp.cart2pol(Z[drone_key]['real'][:,0,:],Z[drone_key]['real'][:,1,:])
    POS_theta = POS_theta*180/np.pi
    POS_azimuth = param[drone_key]['azimuth'] + 90 - POS_theta
    POS_Lat,POS_Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                    POS_azimuth,POS_rho)
    
    # shift_drone = shift[drone_key]
    # POS_Long = POS_Long + shift_drone['longitude']
    # POS_Lat = POS_Lat + shift_drone['latitude']
    Z[drone_key]['GPS'][:,0,:] = POS_Lat
    Z[drone_key]['GPS'][:,1,:] = POS_Long
    
print('Buoys GPS coordinates computed')

#%% Plot buoys GPS position
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()   

for drone_key in Z.keys():
    for i in range(6): 
        current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](i))
        ax.plot(Z[drone_key]['GPS'][i,1,:],Z[drone_key]['GPS'][i,0,:],'.',color = current_color) 
        
ax.set_ylabel(r'Latitude')
ax.set_xlabel(r'Longitude')

figname = f'{fig_folder}Raw_GPS_buoys_3drones'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Save buoys structure 
path2save = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path2save}structure_buoys_tracking_3drones.pkl'

with open(filename,'wb') as pf:
    pickle.dump(Z,pf)


###################################################################################################
#%% Compute optimal procruste operation to apply to one drone to an other - GENERAL
###################################################################################################

synchro = {'fulmar' : 0, 'mesange' : 520, 'bernache' : 520}
Nb_frames = 4000
Nb_buoys = 5
output_rescaling = 1
procruste_op = {}
for ref_drone in Z.keys():
    for projected_drone in Z.keys():
        if projected_drone != ref_drone :
            key_procruste = f'{projected_drone}_2_{ref_drone}'
            print(key_procruste)
    
            current_op = {}
            current_op['i0'] = max(synchro[ref_drone],synchro[projected_drone])
            current_op['Nb_frames'] = Nb_frames
            current_op['ref_drone'] = ref_drone
            current_op['projected_drone'] = projected_drone
            current_op['scaling'] = np.zeros((Nb_frames,1))
            current_op['rot'] = np.zeros((2,2,Nb_frames))
            current_op['translat'] = np.zeros((2,Nb_frames))
            
            # loop over all frames
            for i in range(Nb_frames):  
                points = {}
                for drone_key in [ref_drone,projected_drone]:
                    i0 = synchro[drone_key]
                    # Compute buoys GPS position at frame i
                    points[drone_key] = Z[drone_key]['GPS'][:Nb_buoys,:,i0 + i]
                
                scaling,R,t = procrustes_alignement(points[projected_drone], points[ref_drone],output_rescaling)
                current_op['rot'][:,:,i] = R
                current_op['translat'][:,i] = t
                current_op['scaling'][i] = scaling
            
            # Save current operations
            procruste_op[key_procruste] = current_op
            

file2save = f'{fig_folder}procruste_operation_i0_520_Nbframes_{Nb_frames}_any_drones.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(procruste_op,pf)
    print(f'{file2save} saved !')


#%% Compute (Lat,Long) coordinates using procruste operation

latlong_transfo = np.zeros(Z[ref_drone]['GPS'][:Nb_buoys,:,:Nb_frames].shape)
for i in range(Nb_frames):
    R = procruste_op['rot'][:,:,i]
    t = procruste_op['translat'][:,i]
    scaling = procruste_op['scaling'][i]
    latlong_transfo[:,:,i] = scaling * (R @ Z[projected_drone]['GPS'][:Nb_buoys,:,i].T).T + t


    







###################################################################################
################ FIND BEST TRANSFORMATION bernache / mesange ######################
###################################################################################

#%% Find best transformation for all time steps at once
Nb_frames = 4000
Nb_buoys = 5
ref_drone = 'mesange'
projected_drone = 'bernache'

procruste_op = {}
# procruste_op['scaling'] = np.zeros((Nb_frames,1))
# procruste_op['rot'] = np.zeros((2,2,Nb_frames))
# procruste_op['translat'] = np.zeros((2,Nb_frames))

points = {}
for drone_key in [ref_drone,projected_drone]:
    # Compute buoys GPS position at frame i
    current_points = Z[drone_key]['GPS'][:Nb_buoys,:,:Nb_frames]
    transposed = np.transpose(current_points,axes = (0,2,1))
    points[drone_key] = np.reshape(transposed,(transposed.shape[0]*transposed.shape[1],transposed.shape[2]))
        
scaling,R,t = procrustes_alignement(points[projected_drone], points[ref_drone],output_rescaling = 0)
procruste_op['rot'] = R
procruste_op['translat'] = t
procruste_op['scaling'] = scaling

#%% Compute (Lat,Long) coordinates using procruste operation

latlong_transfo = np.zeros(Z[ref_drone]['GPS'][:Nb_buoys,:,:Nb_frames].shape)
R = procruste_op['rot']
t = procruste_op['translat']
scaling = procruste_op['scaling']

for i in range(Nb_frames):
    latlong_transfo[:,:,i] = scaling * (R @ Z[projected_drone]['GPS'][:Nb_buoys,:,i].T).T + t
    

#%% Find best translation and rotation for each time step

Nb_frames = 4000 #align first Nb_frames
Nb_buoys = 5
ref_drone = 'mesange'
projected_drone = 'bernache'
output_rescaling = 1

procruste_op = {}
procruste_op['i0'] = 0
procruste_op['Nb_frames'] = Nb_frames
procruste_op['ref_drone'] = ref_drone
procruste_op['projected_drone'] = projected_drone
procruste_op['scaling'] = np.zeros((Nb_frames,1))
procruste_op['rot'] = np.zeros((2,2,Nb_frames))
procruste_op['translat'] = np.zeros((2,Nb_frames))

for i in range(Nb_frames):  
    points = {}
    for drone_key in [ref_drone,projected_drone]:
        # Compute buoys GPS position at frame i
        points[drone_key] = Z[drone_key]['GPS'][:Nb_buoys,:,i]
    
    scaling,R,t = procrustes_alignement(points[projected_drone], points[ref_drone],output_rescaling)
    procruste_op['rot'][:,:,i] = R
    procruste_op['translat'][:,i] = t
    procruste_op['scaling'][i] = scaling

file2save = f'{fig_folder}procruste_operation_Nbframes_{Nb_frames}_ref_drone_{ref_drone}_projected_drone_{projected_drone}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(procruste_op,pf)
    print(f'{file2save} saved !')


#%% Compute (Lat,Long) coordinates using procruste operation

# load procruste operations
file2load = f'{fig_folder}procruste_operation_Nbframes_{Nb_frames}_ref_drone_{ref_drone}_projected_drone_{projected_drone}.pkl'
with open(file2load,'rb') as pf :
    procruste_op = pickle.load(pf)
    print(f'{file2load} loaded !')


latlong_transfo = np.zeros(Z[ref_drone]['GPS'][:Nb_buoys,:,:Nb_frames].shape)
for i in range(Nb_frames):
    R = procruste_op['rot'][:,:,i]
    t = procruste_op['translat'][:,i]
    scaling = procruste_op['scaling'][i]
    latlong_transfo[:,:,i] = scaling * (R @ Z[projected_drone]['GPS'][:Nb_buoys,:,i].T).T + t


#%% Compare Lat,Long coordinates with image in background 
if ref_drone == 'mesange':
    path2img = f'{root}Initial_images/mesange_im_11496.tiff'
elif ref_drone == 'bernache':
    path2img = f'{root}Initial_images/bernache_im_11500.tiff'

set_graphs.set_matplotlib_param('single')

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
# create coordinates for image
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param[ref_drone]['h'],param[ref_drone]['alpha_0'],
                                       param[ref_drone]['focale'])

fig, ax = plt.subplots()
# Plot frame 
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)

Lat,Long = XY2GPS(Xreal,Yreal,Lat0,Long0,param[ref_drone]['azimuth'])
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm_time'](Z[ref_drone]['t'][:Nb_frames])) 
    ax.scatter(Z[ref_drone]['GPS'][i,1,:Nb_frames],Z[ref_drone]['GPS'][i,0,:Nb_frames],s = 20,marker = '.',color = current_color)
    
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm_time'](Z[projected_drone]['t'][:Nb_frames])) 
    ax.scatter(latlong_transfo[i,1,:],latlong_transfo[i,0,:],s = 20,marker = '.',color = current_color)
    
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}GPS_positions_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Superpose buoys position in X,Y mesange coordinates system
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)

# coordinates of central point in ref_drone coordinates system
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)
x,y = GPS2XY(latlong_transfo[:,0,:],latlong_transfo[:,1,:],Lat0,Long0,param[ref_drone]['azimuth'])

for i in range(5):
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm_time'](Z[ref_drone]['t'][:Nb_frames]))
    ax.scatter(Z[ref_drone]['real'][i,0,:Nb_frames],Z[ref_drone]['real'][i,1,:Nb_frames],s = 20,marker = '.',
               color = current_color)
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm_time'](Z[projected_drone]['t'][:Nb_frames]))
    ax.scatter(x[i,:],y[i,:],s = 20,marker = '.',color = current_color)

ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')

figname = f'{fig_folder}XY_positions_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Compute mean distance between two sets of points
Xref = Z[ref_drone]['real'][:Nb_buoys,0,:Nb_frames]
Yref = Z[ref_drone]['real'][:Nb_buoys,1,:Nb_frames]

error = np.sqrt((x - Xref)**2 + (y - Yref)**2)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
for i in range(5):    
    label_buoy = f'Buoy {i+1}'
    ax.plot(error[i,:],label = label_buoy)

ax.legend()
ax.set_ylabel(r'Error (m)')
ax.set_xlabel(r'Time step')

figname = f'{fig_folder}Error_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')
########################################################################################################################
########################################################################################################################






#%%
########################################################################################################################
################################ FIND BEST TRANSFORMATION fulmar -> mesange ############################################
########################################################################################################################

#%%

synchro = {'fulmar' : 0, 'mesange' : 520, 'bernache' : 520}
projected_drone = 'fulmar'
ref_drone = 'mesange'
Nb_buoys = 5
Nb_frames = 3000
output_rescaling = 1 #enable homothety

procruste_op = {}
procruste_op['i0'] = 0
procruste_op['Nb_frames'] = Nb_frames
procruste_op['ref_drone'] = ref_drone
procruste_op['projected_drone'] = projected_drone

procruste_op['scaling'] = np.zeros((Nb_frames,1))
procruste_op['rot'] = np.zeros((2,2,Nb_frames))
procruste_op['translat'] = np.zeros((2,Nb_frames))

for i in range(Nb_frames):  
    points = {}
    for drone_key in [ref_drone,projected_drone]:
        # Compute buoys GPS position at frame i
        i0 = synchro[drone_key]
        points[drone_key] = Z[drone_key]['GPS'][:Nb_buoys,:,i0 + i]
    
    scaling,R,t = procrustes_alignement(points[projected_drone], points[ref_drone],output_rescaling)
    procruste_op['rot'][:,:,i] = R
    procruste_op['translat'][:,i] = t
    procruste_op['scaling'][i] = scaling

file2save = f'{fig_folder}procruste_operation_Nbframes_{Nb_frames}_ref_drone_{ref_drone}_projected_drone_{projected_drone}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(procruste_op,pf)
    print(f'{file2save} saved !')

#%% Compute (Lat,Long) coordinates using procruste operation

latlong_transfo = np.zeros((Nb_buoys,2,Nb_frames))
i0 = synchro[projected_drone]
for i in range(Nb_frames):
    R = procruste_op['rot'][:,:,i]
    t = procruste_op['translat'][:,i]
    scaling = procruste_op['scaling'][i]
    latlong_transfo[:,:,i] = scaling * (R @ Z[projected_drone]['GPS'][:Nb_buoys,:,i0 + i].T).T + t

#%% Compare Lat,Long coordinates with image in background 
if ref_drone == 'mesange':
    path2img = f'{root}Initial_images/mesange_im_11496.tiff'
elif ref_drone == 'bernache':
    path2img = f'{root}Initial_images/bernache_im_11500.tiff'

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
# create coordinates for image
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param[ref_drone]['h'],param[ref_drone]['alpha_0'],
                                       param[ref_drone]['focale'])

fig, ax = plt.subplots()
# Plot frame 
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)

Lat,Long = XY2GPS(Xreal,Yreal,Lat0,Long0,param[ref_drone]['azimuth'])
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    i0 = synchro[ref_drone]
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](Z[ref_drone]['t'][:Nb_frames])) 
    ax.scatter(Z[ref_drone]['GPS'][i,1,i0:i0 + Nb_frames],Z[ref_drone]['GPS'][i,0,i0:i0 + Nb_frames],
               s = 20,marker = '.',color = current_color)
    
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](Z[projected_drone]['t'][:Nb_frames])) 
    ax.scatter(latlong_transfo[i,1,:],latlong_transfo[i,0,:],s = 20,marker = '.',color = current_color)
    
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}GPS_positions_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Superpose buoys position in X,Y mesange coordinates system

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)

# coordinates of central point in ref_drone coordinates system
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)
x,y = GPS2XY(latlong_transfo[:,0,:],latlong_transfo[:,1,:],Lat0,Long0,param[ref_drone]['azimuth'])

for i in range(5):
    i0 = synchro[ref_drone]
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](Z[ref_drone]['t'][:Nb_frames]))
    ax.scatter(Z[ref_drone]['real'][i,0,i0:i0 + Nb_frames],Z[ref_drone]['real'][i,1,i0:i0 + Nb_frames],s = 20,marker = '.',
               color = current_color)
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](Z[projected_drone]['t'][:Nb_frames]))
    ax.scatter(x[i,:],y[i,:],s = 20,marker = '.',color = current_color)

ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')

figname = f'{fig_folder}XY_positions_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Compute mean distance between two sets of points
i0 = synchro[ref_drone]
Xref = Z[ref_drone]['real'][:Nb_buoys,0,i0:i0+Nb_frames]
Yref = Z[ref_drone]['real'][:Nb_buoys,1,i0:i0+Nb_frames]

error = np.sqrt((x - Xref)**2 + (y - Yref)**2)

fig, ax = plt.subplots(figsize = (12,9))
for i in range(5):    
    label_buoy = f'Buoy {i}'
    ax.plot(error[i,:],label = label_buoy)

ax.legend()
ax.set_ylabel(r'Error (m)')
ax.set_xlabel(r'Time step')

figname = f'{fig_folder}Error_{projected_drone}_to_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%%
######################################################################################################################
######################################################################################################################



####################################################################################################################
####### Illustrate procruste transformation ########################################################################
####################################################################################################################
synchro = {'fulmar' : 0, 'mesange' : 520, 'bernache' : 520}
projected_drone = 'bernache'
ref_drone = 'mesange'
Nb_buoys = 5
output_rescaling = 1

# check procrustes alignement 
i = 0
i0 = synchro[ref_drone]
points = {}
for drone_key in [ref_drone,projected_drone]:
    # Compute buoys GPS position at frame i
    i0 = synchro[drone_key]
    points[drone_key] = Z[drone_key]['GPS'][:Nb_buoys,:,i0 + i]

points_project = points[projected_drone]
points_ref = points[ref_drone]

scale_factor,R,t =  procrustes_alignement(points_project,points_ref,output_rescaling = 1)

# apply procruste operations to Lat,Long
rotated = R @ points_project.T
scaled = scale_factor * rotated.T
latlong_transfo = scaled + t[None,...]

print(latlong_transfo.shape)
# # reshape matrix 
# latlong_transfo = np.reshape(latlong_transfo,(X.shape[0],X.shape[1],2))

#%% Plot raw positions 

if ref_drone == 'mesange':
    path2img = f'{root}Initial_images/mesange_im_11496.tiff'
elif ref_drone == 'bernache':
    path2img = f'{root}Initial_images/bernache_im_11500.tiff'

img = cv.imread(path2img)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
# create coordinates for image
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param[ref_drone]['h'],param[ref_drone]['alpha_0'],
                                       param[ref_drone]['focale'])

set_graphs.set_matplotlib_param('single')
mksize = 7
fig, ax = plt.subplots()
# Plot frame 
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)

Lat,Long = XY2GPS(Xreal,Yreal,Lat0,Long0,param[ref_drone]['azimuth'])
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    i0 = synchro[ref_drone]
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i)) 
    ax.plot(points_ref[i,1],points_ref[i,0],marker = 'o',color = current_color,
               markeredgecolor = 'k',markersize = mksize)
    
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i)) 
    ax.plot(points_project[i,1],points_project[i,0],marker = 'o',color = current_color,
               markeredgecolor = 'k',markersize = mksize)
    
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}illustration_procruste_operations_raw_positions_projected_{projected_drone}_ref_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Plot final positions 

fig, ax = plt.subplots()
# Plot frame 
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)

Lat,Long = XY2GPS(Xreal,Yreal,Lat0,Long0,param[ref_drone]['azimuth'])
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    i0 = synchro[ref_drone]
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i)) 
    ax.plot(points_ref[i,1],points_ref[i,0],marker = 'o',color = current_color,
               markeredgecolor = 'k',markersize = mksize)
    
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i)) 
    ax.plot(latlong_transfo[i,1],latlong_transfo[i,0],marker = 'o',color = current_color,
               markeredgecolor = 'k',markersize = mksize)
    
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}illustration_procruste_operations_finale_projected_{projected_drone}_ref_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


#%% Detailed computation of each operation 

# Compute centers
c_proj = np.nanmean(points_project,axis = 0)
c_ref = np.nanmean(points_ref,axis = 0)
# reduce center position
centered_proj = points_project - c_proj
centered_ref = points_ref - c_ref

# compute rescaling factor if output_rescaling = 1
if output_rescaling :
    scale_factor = np.sqrt(np.sum(centered_ref**2)/np.sum(centered_proj**2))
else :
    scale_factor = 1

# Compute best rotation
H = centered_proj.T @ centered_ref # Cross-covariance matrix
U,S,Vt = np.linalg.svd(H) # SVD
R = Vt.T @ U.T # Compute optimal rotation

# Ensure a proper rotation (no reflection)
if np.linalg.det(R) < 0:
    Vt[-1,:] *= -1
    R = Vt.T @ U.T  
    
# Compute best translation 
t = c_ref - scale_factor * R @ c_proj


#%%




#%% Plots to illustrate procruste operations

set_graphs.set_matplotlib_param('single')

mksize = 8
fig, ax = plt.subplots()
for i in range(5):
    ref_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i))
    ax.plot(points_ref[i,1],points_ref[i,0],'o',color = ref_color,markeredgecolor = 'k',markersize = mksize)
    project_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i))
    ax.plot(points_project[i,1],points_project[i,0],'o',color = project_color,markeredgecolor = 'k',
            markersize = mksize)
    
ax.plot(c_proj[1],c_ref[0],'rd',markeredgecolor = 'k',markersize = mksize)
ax.plot(c_ref[1],c_ref[0],'gd',markeredgecolor = 'k',markersize = mksize)

ax.set_ylabel(r'Latitude')
ax.set_xlabel(r'Longitude')

figname = f'{fig_folder}Illustration_procruste_mass_center_ref_{ref_drone}_projected_{projected_drone}_frame_520'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%%

set_graphs.set_matplotlib_param('single')

mksize = 8
fig, ax = plt.subplots()
for i in range(5):
    ref_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i))
    ax.plot(centered_ref[i,1],centered_ref[i,0],'o',color = ref_color,markeredgecolor = 'k',markersize = mksize)
    project_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i))
    ax.plot(centered_proj[i,1],centered_proj[i,0],'o',color = project_color,markeredgecolor = 'k',
            markersize = mksize)

ax.set_ylabel(r'Latitude')
ax.set_xlabel(r'Longitude')

figname = f'{fig_folder}Illustration_procruste_centered_points_ref_{ref_drone}_projected_{projected_drone}_frame_520'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')



#%%
scaled = scale_factor * centered_proj
set_graphs.set_matplotlib_param('single')

mksize = 8
fig, ax = plt.subplots()
for i in range(5):
    ref_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i))
    ax.plot(centered_ref[i,1],centered_ref[i,0],'o',color = ref_color,markeredgecolor = 'k',markersize = mksize)
    project_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i))
    ax.plot(scaled[i,1],scaled[i,0],'o',color = project_color,markeredgecolor = 'k',
            markersize = mksize)

ax.set_ylabel(r'Latitude')
ax.set_xlabel(r'Longitude')

figname = f'{fig_folder}Illustration_procruste_scaled_points_ref_{ref_drone}_projected_{projected_drone}_frame_520'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


#%%
rotated = R @ scaled.T
rotated = rotated.T

set_graphs.set_matplotlib_param('single')

mksize = 8
fig, ax = plt.subplots()
for i in range(5):
    ref_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](i))
    ax.plot(centered_ref[i,1],centered_ref[i,0],'o',color = ref_color,markeredgecolor = 'k',markersize = mksize)
    project_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](i))
    ax.plot(rotated[i,1],rotated[i,0],'o',color = project_color,markeredgecolor = 'k',
            markersize = mksize)

ax.set_ylabel(r'Latitude')
ax.set_xlabel(r'Longitude')

figname = f'{fig_folder}Illusatration_procruste_rotated_points_ref_{ref_drone}_projected_{projected_drone}_frame_520'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')




########################################################################################################
############################ TRANSFORMATION BY HAND ####################################################
########################################################################################################
#%% Apply transformation to each drone to have a common framework to all drones
folder = f'{root}transformations/'
filetransfo = f'{folder}bernache_fulmar2mesange_transfo.pkl'

# dictionnary for transformations to apply to each drone
transfo = {'fulmar':{},'bernache':{},'mesange':{}}

# Time synchronization (reference : Fulmar)
ilag = 520 # or 521
transfo['fulmar']['ilag'] = 0
transfo['bernache']['ilag'] = 520
transfo['mesange']['ilag'] = 520

# Homothetie factor for each drone 
transfo['fulmar']['fh'] = 0.9
transfo['bernache']['fh'] = -1
transfo['mesange']['fh'] = 1

# Shift in (X,Y) coordinates
transfo['fulmar']['delta_x'] = -6.8
transfo['bernache']['delta_x'] = -0.5
transfo['mesange']['delta_x'] = 0

transfo['fulmar']['delta_y'] = 11.7
transfo['bernache']['delta_y'] = 26.5
transfo['mesange']['delta_y'] = 0

# Rotation in (X,Y) coordinates 
transfo['fulmar']['rot'] = -0.02
transfo['bernache']['rot'] = -0.06
transfo['mesange']['rot'] = 0

with open(filetransfo,'wb') as pf :
    pickle.dump(transfo,pf)
    print(f'{filetransfo} saved !')

# for drone_key in Z.keys():
#     ilag = transfo[drone_key]['ilag']
#     Z[drone_key]['buoys'] = Z[drone_key]['buoys'][:,:,ilag:]
#     Z[drone_key]['real'] = Z[drone_key]['real'][:,:,ilag:]
#     Z[drone_key]['t'] = Z[drone_key]['t'][ilag:]
#     Z[drone_key]['idx'] = Z[drone_key]['idx'][ilag:]

# print('Drone synchronized')


#%% Time synchronization with Fulmar

fps = 30

fig,ax = plt.subplots()
for drone_key in Z.keys():
    for i in range(4):
        ilag = transfo[drone_key]['ilag']
        print(ilag)
        X_buoy = Z[drone_key]['real'][i,0,:]*transfo[drone_key]['fh']
        t = Z[drone_key]['t'][:] - Z[drone_key]['t'][0] 
        
        current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](i))
        ax.plot(t,X_buoy,'.',color = current_color)

ax.set_xlabel('$t$')
ax.set_ylabel('$X$')

figname = f'{fig_folder}Drones_synchro'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%%

# Plot (X,Y) coordinates after transformations
fig, ax = plt.subplots()
for drone_key in Z.keys():
    
    # homothetie
    fh = transfo[drone_key]['fh'] 
    X = Z[drone_key]['real'][:,0,:]*fh
    Y = Z[drone_key]['real'][:,1,:]*fh
    
    # Rotation
    theta = transfo[drone_key]['rot']
    X = X*np.cos(theta) + Y*np.sin(theta)
    Y = -X*np.sin(theta) + Y*np.cos(theta)
    
    # translation
    delta_x = transfo[drone_key]['delta_x']
    delta_y = transfo[drone_key]['delta_y']
    X = X + delta_x
    Y = Y + delta_y
    
    # time synchro with Fulmar
    ilag = transfo[drone_key]['ilag']
        
    for i in range(6):

        current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](i))
        label_buoy = f'Buoy {i+1}'
        
        # Take buoy 5 as reference
        # X_buoy = X[i,:] - X[4,:]
        # Y_buoy = Y[i,:] - Y[4,:]
        
        X_buoy = X[i,:]
        Y_buoy = Y[i,:]
        ax.plot(X_buoy,Y_buoy,'.',color = current_color,label = label_buoy)

ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_xlim([-90,70])
ax.set_ylim([-40,20])

# ax.set_xlim([-30,120])
# ax.set_ylim([-20,10])
#%% Superpose detection with image

# Create a colormap based on time 
maps = {'fulmar':'Blues','bernache':'Oranges','mesange':'Greens'}
maps_drone = {}
for key_drone in maps.keys():
    maps_drone[key_drone] = {}
    full_map = mpl.colormaps[maps[key_drone]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    cnorm = mcolors.Normalize(vmin = Z[key_drone]['t'][0], vmax = Z[key_drone]['t'][-1])
    
    maps_drone[key_drone]['new_map'] = new_map
    maps_drone[key_drone]['cnorm'] = cnorm


# create coordinates for image
[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param['mesange']['h'],param['mesange']['alpha_0'],
                                       param['mesange']['focale'])

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)

for drone_key in Z.keys():
    
    # homothetie
    fh = transfo[drone_key]['fh'] 
    X = Z[drone_key]['real'][:,0,:]*fh
    Y = Z[drone_key]['real'][:,1,:]*fh
    
    # Rotation
    theta = transfo[drone_key]['rot']
    X = X*np.cos(theta) + Y*np.sin(theta)
    Y = -X*np.sin(theta) + Y*np.cos(theta)
    
    # translation
    delta_x = transfo[drone_key]['delta_x']
    delta_y = transfo[drone_key]['delta_y']
    X = X + delta_x
    Y = Y + delta_y
    
    # time synchro with Fulmar
    ilag = transfo[drone_key]['ilag']
        
    for i in range(6):

        current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](Z[drone_key]['t']))
        label_buoy = f'Buoy {i+1}'
       
        X_buoy = X[i,:]
        Y_buoy = Y[i,:]
        ax.scatter(X_buoy,Y_buoy,s = 20,marker = '.',color = current_color,label = label_buoy)

ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_xlim([-90,90])
ax.set_ylim([-40,60])


figname = f'{fig_folder}Transformation_bernache_fulmar2mesange'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')






























#%% Compute for each time event, center of all buoys coordinates

center = {}
Nb_frames = 1000 #align first Nb_frames
Nb_buoys = 5
for drone_key in ['bernache','mesange']:
    center[drone_key] = np.zeros((Z[drone_key]['GPS'].shape[1],Nb_frames))
    for i in range(Nb_frames):
        center[drone_key][0,i] = np.nanmean(Z[drone_key]['GPS'][:Nb_buoys,0,i],axis = 0)
        center[drone_key][1,i] = np.nanmean(Z[drone_key]['GPS'][:Nb_buoys,1,i],axis = 0)
    
# find best translation
best_t = center['mesange'] - center['bernache']

#%% 
fig, ax = plt.subplots()
ax.plot(best_t[0,:],'.')

fig, ax = plt.subplots()
ax.plot(best_t[1,:],'.')

#%%
# Create a colormap based on time 
maps = {'fulmar':'Blues','bernache':'Oranges','mesange':'Greens'}
maps_drone = {}
for key_drone in maps.keys():
    maps_drone[key_drone] = {}
    full_map = mpl.colormaps[maps[key_drone]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    ilag = transfo[key_drone]['ilag']
    cnorm = mcolors.Normalize(vmin = Z[key_drone]['t'][0], vmax = Z[key_drone]['t'][-1])
    
    maps_drone[key_drone]['new_map'] = new_map
    maps_drone[key_drone]['cnorm'] = cnorm
    
fig, ax = plt.subplots()
# Plot frame 
dist2drone = param['mesange']['h']/np.tan(param['mesange']['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param['mesange']['latitude'],param['mesange']['longitude'],
                                                param['mesange']['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param['mesange']['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)


for i in range(Nb_frames):
    drone_key ='bernache'
    current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](Z[drone_key]['t'][i]))
    ax.scatter(Z[drone_key]['GPS'][:Nb_buoys,1,i]+best_t[1,i],Z[drone_key]['GPS'][:Nb_buoys,0,i]+best_t[0,i],
               s = 20,marker = '.',color = current_color)
    drone_key = 'mesange'
    current_color = maps_drone[drone_key]['new_map'](maps_drone[drone_key]['cnorm'](Z[drone_key]['t'][i]))
    ax.scatter(Z[drone_key]['GPS'][:Nb_buoys,1,i],Z[drone_key]['GPS'][:Nb_buoys,0,i],
               s = 20,marker = '.',color = current_color)

ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

ax.plot(Long0,Lat0,'r.')
#%% Convert bernache GPS coordinates to polar coordinates in mesange coordinates system 

ref_drone = 'mesange'
projected_drone = 'bernache'
R_earth = 6371e3
Nb_buoys = 5
Nb_frames = 1000

# coordinates of central point in ref_drone coordinates system
dist2drone = param[ref_drone]['h']/np.tan(param[ref_drone]['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param[ref_drone]['latitude'],param[ref_drone]['longitude'],
                                                param[ref_drone]['azimuth'],dist2drone)

# Compute polar coordinates 
lat = Z[projected_drone]['GPS'][:Nb_buoys,0,:Nb_frames] + best_t[0,:]
long = Z[projected_drone]['GPS'][:Nb_buoys,1,:Nb_frames] + best_t[1,:]
rho = R_earth*np.sqrt(((lat - Lat0)*np.pi/180)**2 + np.cos(Lat0*np.pi/180)**2 * ((long - Long0)*np.pi/180)**2)

psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(long - Long0)*np.pi/180,(lat - Lat0)*np.pi/180)*180/np.pi

POS_Lat,POS_Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                psi,rho)

theta = param[ref_drone]['azimuth'] - psi + 90
# compute euclidean coordinates
x,y = dp.pol2cart(rho,theta*np.pi/180)



    

#%% Check GPS positions

fig, ax = plt.subplots()
dist2drone = param['mesange']['h']/np.tan(param['mesange']['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param['mesange']['latitude'],param['mesange']['longitude'],
                                                param['mesange']['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param['mesange']['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](Z[ref_drone]['t'][:Nb_frames])) 
    ax.scatter(Z[ref_drone]['GPS'][i,1,:Nb_frames],Z[ref_drone]['GPS'][i,0,:Nb_frames],s = 20,marker = '.',color = current_color)
    
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](Z[projected_drone]['t'][:Nb_frames])) 
    ax.scatter(POS_Long[i,:],POS_Lat[i,:],s = 20,marker = '.',color = current_color)
    
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x
#%% Check X,Y position

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)

for i in range(5):
    current_color = maps_drone[ref_drone]['new_map'](maps_drone[ref_drone]['cnorm'](Z[ref_drone]['t'][:Nb_frames]))
    ax.scatter(Z[ref_drone]['real'][i,0,:Nb_frames],Z[ref_drone]['real'][i,1,:Nb_frames],s = 20,marker = '.',
               color = current_color)
    current_color = maps_drone[projected_drone]['new_map'](maps_drone[projected_drone]['cnorm'](Z[projected_drone]['t'][:Nb_frames]))
    ax.scatter(x[i,:],y[i,:],s = 20,marker = '.',color = current_color)





#%% Plot buoys position in Latitude / Longitude coordinates 

# Create a colormap based on time 
maps = {'fulmar':'Blues','bernache':'Oranges','mesange':'Greens'}
maps_drone = {}
for key_drone in maps.keys():
    maps_drone[key_drone] = {}
    full_map = mpl.colormaps[maps[key_drone]].resampled(256)
    new_map = mcolors.ListedColormap(full_map(np.linspace(0.2,1,256)))
    ilag = transfo[key_drone]['ilag']
    cnorm = mcolors.Normalize(vmin = Z[key_drone]['t'][ilag], vmax = Z[key_drone]['t'][-1])
    
    maps_drone[key_drone]['new_map'] = new_map
    maps_drone[key_drone]['cnorm'] = cnorm

fig, ax = plt.subplots()

# Plot image with longitude,latitude coordinates
# horizontal distance between center of metric coordinates and drone position 
dist2drone = param['mesange']['h']/np.tan(param['mesange']['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param['mesange']['latitude'],param['mesange']['longitude'],
                                                param['mesange']['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param['mesange']['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)

c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)



Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param['mesange']['latitude'],param['mesange']['longitude'],
                                                param['mesange']['azimuth'],dist2drone)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}transfo_LatLong_fulmarbernache2mesange'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


#%% Find best rotation for first 1000 frames

# nb frames to match
Nb_frames = 1000
rot = np.arange(-0.2,0.2,0.001)
Z0 = Z['mesange']['real']
Z1 = Z['bernache']['real']*(-1)
buoy_ref = 1

d = best_rotation(rot,Z0[:,:,:Nb_frames],Z1[:,:,:Nb_frames],buoy_ref)
min_index = np.argmin(d)
# best rotation 
best_rot = rot[min_index]
print(f'Best rotation : {best_rot}')











