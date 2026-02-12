# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 12:56:30 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
from scipy.io import loadmat
from scipy.interpolate import griddata
import csv 

import icewave.tools.matlab2python as mat2py
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.matlab_colormaps as matcmaps

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# PARULA COLORMAP 
parula_map = matcmaps.parula()

#%%
date = '0206'
drone_ID = 'bernache'
exp_ID = '06-situation_map_geophones_002'

path2data = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/movie_1/'
path2img = f'U:/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/movie_1/'

fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(f'{path2img}*.tiff')
file2read = filelist[0]

#%% Set parameters of drone 

param_dict = {}
param_dict['h_drone'] = 46.4 # drone height
param_dict['alpha_0'] = 37.4*np.pi/180 # rad
param_dict['focale'] = 2700 # focal length 

param_dict['latitude'] = 48.34689
param_dict['longitude'] = -68.81791
param_dict['azimuth'] = 355.8

# create UTC_t0
local_time = datetime.strptime('2025-02-06 14:48:05.440' + '000','%Y-%m-%d %H:%M:%S.%f')
tz_quebec = pytz.timezone('America/Montreal')
tz_utc = pytz.utc
local_time = tz_quebec.localize(local_time)
UTC_time = local_time.astimezone(tz_utc)

param_dict['UTC_t0'] = UTC_time

file2save = f'{path2data}Param_drone_{date}_{drone_ID}_{exp_ID}_movie_1.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(param_dict,pf)

#%% 

img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

fig,ax = plt.subplots()
ax.imshow(img)


#%% Set pixel coordinates of geophones 

POS = {}

POS['pix'] = np.array([
    [582,1572], 
    [839,1034],
    [1015.5,679],
    [1136,432],
    [1501,1556.6],
    [1582,1028],
    [1619,690],
    [1643,436],
    [2387,1542],
    [2276,1022], #10
    [2203,695], #11    
    [2149,442],
    [3259,1532],
    [2978,1009],
    [2783,676],
    [2642,439]
    ])

#%% Superpose frame and pixel positions of geophones

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.imshow(img)
ax.plot(POS['pix'][:,0],POS['pix'][:,1],'o',color = 'r', markeredgecolor = 'k')

# hide labels
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

figname = f'{fig_folder}situation_picture_tomo_0206_2025_coord_pixel'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Georectify image and positions 

[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['h_drone'],param_dict['alpha_0'],
                                       param_dict['focale'])

# compute object position
POS_realx,POS_realy = dp.projection_real_space(POS['pix'][:,0],POS['pix'][:,1],x0,y0,
                                               param_dict['h_drone'],param_dict['alpha_0'],param_dict['focale'])

POS['real'] = np.zeros(POS['pix'].shape)
POS['real'][:,0] = POS_realx
POS['real'][:,1] = POS_realy


#%% Show georectified image 

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
ax.plot(POS['real'][:,0],POS['real'][:,1],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 7)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

c.set_rasterized(True)

figname = f'{fig_folder}georectified_situation_picture_tomo_0206_2025_coord_meter'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Show image using GPS coordinates 

# horizontal distance between center of metric coordinates and drone position 
dist2drone = param_dict['h_drone']/np.tan(param_dict['alpha_0'])
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param_dict['latitude'],param_dict['longitude'],
                                                param_dict['azimuth'],dist2drone)

Lat,Long = dp.XY2GPS(Xreal,Yreal,Lat0,Long0,param_dict['azimuth'])

POS_Lat,POS_Long = dp.XY2GPS(POS['real'][:,0],POS['real'][:,1],Lat0,Long0,param_dict['azimuth'])
POS['GPS'] = np.zeros(POS['pix'].shape)
POS['GPS'][:,0] = POS_Lat
POS['GPS'][:,1] = POS_Long

#%% Show image in GPS coordinates 

fig,ax = plt.subplots()
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS['GPS'][:,1],POS['GPS'][:,0],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 5)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 5)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

c.set_rasterized(True)

figname = f'{fig_folder}georeferenced_situation_picture_tomo_0206_2025_coord_GPS'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Load ice thickness prelevment

maddie_gps = {}
maddie_gps['GPS'] = np.array([ 
    [48.347354 ,-68.818249],
    [48.347492, -68.818232],
    [48.347638, -68.818203],
    [48.347664, -68.817955],
    [48.347546, -68.817963],
    [48.347346, -68.817941],
    [48.34735, -68.817657],
    [48.347539, -68.817639],
    [48.347697, -68.817693],
 ]) 

maddie_gps['image_gps'] = np.zeros((3,3,2))
for i in range(3):
    for j in range(3):
        idx = j+i*4
        print(idx)
        maddie_gps['image_gps'][j,i,:] = np.mean(np.array([
            POS['GPS'][idx,:], POS['GPS'][idx + 1,:], POS['GPS'][idx + 4, : ], POS['GPS'][idx + 5, :]]),axis = 0)

maddie_gps['image_gps'] = maddie_gps['image_gps'].reshape(maddie_gps['GPS'].shape)

maddie_gps['h'] = np.array([57,57,51,62,59,68,66,71,66])
maddie_gps['free_board'] = np.array([5,6,4,7,6,10,5,6,6])


filename = f'{path2data}maddie_gps_tomography_data_{date}_2025.pkl'
with open(filename,'wb') as pf:
    pickle.dump(maddie_gps,pf)
    
# additional thickness measurements (geophone line)
additional_gps = {}
additional_gps['GPS'] = np.array([
    [48.347233,-68.818358],
    [48.347215,-68.818419],
    [48.347211,-68.818538],
    [48.347202,-68.818596],
    [48.347192,-68.818711],
    [48.347169,-68.81873],
    [48.347193,-68.818825],
    [48.347188,-68.818888],
    [48.347178,-68.818951],])
            
additional_gps['h'] = np.array([52,47,61,55,54,50,54,59,56])
additional_gps['free_board'] = np.array([3,5,4,6,4,5,5,6,5])


filename = f'{path2data}maddie_gps_tomography_data_{date}_2025_additional_positions.pkl'
with open(filename,'wb') as pf:
    pickle.dump(additional_gps,pf)
    
#%% Show geophone array and thickness measurements

fig,ax = plt.subplots()
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS['GPS'][:,1],POS['GPS'][:,0],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.plot(maddie_gps['artificial_gps'][:,1],maddie_gps['artificial_gps'][:,0],'.',markerfacecolor = 'tab:blue',
        markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 5)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 5)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

c.set_rasterized(True)


#%% interpolate

key_gps = 'artificial_gps'

min_long = np.min(maddie_gps[key_gps] [:,1])
max_long = np.max(maddie_gps[key_gps] [:,1])
min_lat = np.min(maddie_gps[key_gps] [:,0])
max_lat = np.max(maddie_gps[key_gps] [:,0])


long_array = np.linspace(min_long,max_long,20)
lat_array = np.linspace(min_lat,max_lat,20)

grid_long, grid_lat = np.meshgrid(long_array,lat_array)
grid_h = griddata(maddie_gps['artificial_gps'],maddie_gps['h'],(grid_lat,grid_long),method = 'linear')


fig, ax = plt.subplots()
c = ax.imshow(grid_h.T, origin = 'lower', extent = [min_long,max_long,min_lat,max_lat], cmap = parula_map)
ax.plot(POS['GPS'][:,1],POS['GPS'][:,0],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 5)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 5)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.01)
cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$h \; \mathrm{(cm)}$')

figname = f'{fig_folder}geophone_thickness_interpolation_tomo_0206_2025_coord_GPS'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot all thickness measurements 

fig,ax = plt.subplots()
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
ax.plot(POS['GPS'][:,1],POS['GPS'][:,0],'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 8)
ax.plot(additional_gps['GPS'][:,1],additional_gps['GPS'][:,0],'.',markerfacecolor = 'tab:blue',
        markeredgecolor = 'k',markersize = 8)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 5)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 5)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

c.set_rasterized(True)


#%%

filename = f'{path2data}maddie_gps_tomography_data_{date}_2025_additional_positions.pkl'
with open(filename,'rb') as pf:
    test = pickle.load(pf)

#%% Load pickles file 

filename = 'C:/Users/sebas/Downloads/maddie_gps_tomography_data_0206_2025_additional_positions.pkl'
with open(filename,'rb') as pf:
    test = pickle.load(pf)


