# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 16:30:15 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
import csv

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Set fig_folder path 

fig_folder = 'F:/Rimouski_2025/Data/Summary/DAS/Situation_tomo_0210/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load DAS GPS position 

base = 'F:/Rimouski_2025/Data/'

file_water_height = f'{base}0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

#%% Load Stephane and Maddie's GPS measurements

# with pickle
filename = 'F:/Rimouski_2025/Data/Summary/DAS/Stephane_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    gps_results = pickle.load(pf)
    
filename = 'F:/Rimouski_2025/Data/Summary/DAS/Maddie_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    MS_gps = pickle.load(pf)
    
#%% Load geophones gps positions
file2load = 'F:/Rimouski_2025/Data/Summary/DAS/Geophones_positions_0210_0212.pkl'
with open(file2load,'rb') as pf:
    GPS_geoph = pickle.load(pf)


#%% Load image

base = 'F:/Rimouski_2025/PIV_images/'
date = '0210'
drone_ID = 'bernache'
exp_ID = '11-situation_map_geophone_003'
movie_ID = 'movie_1'

path2tiff = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/{movie_ID}/'
filelist = glob.glob(f'{path2tiff}*.tiff')
file2read = filelist[0]
print(file2read)

img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

#%% Set drone parameters 

# Set drone parameters
param_dict = {}
param_dict['H'] = 91.8 #110.8
param_dict['alpha_0'] = 90
param_dict['focale'] = 2700 

param_dict['latitude'] =  48.34755# 48.3475505902513 
param_dict['longitude'] =  -68.81794# -68.8179351905923 
param_dict['azimuth'] = 351.5

#%% Georectify picture 

[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['H'],param_dict['alpha_0']*np.pi/180,
                                       param_dict['focale'])

#%% Compute image GPS position 

# horizontal distance between center of metric coordinates and drone position 
dist2drone = param_dict['H']/np.tan(param_dict['alpha_0']*np.pi/180)
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param_dict['latitude'],param_dict['longitude'],
                                                param_dict['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param_dict['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)

#%% Check GPS positioning of drone image and GPS_geophones

set_graphs.set_matplotlib_param('single')
fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

# plot geophones
norm_acq = colors.Normalize(vmin = 0,vmax = 3)
full_cmap = mpl.colormaps['Greens'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,0.8,256)))

for j in range(3,4):
    current_color = cmap(norm_acq(j))
    if j == 0:
        ax.plot(GPS_geoph[:,j,0],GPS_geoph[:,j,1],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3,label = 'Geophone')
    else :
        ax.plot(GPS_geoph[:,j,0],GPS_geoph[:,j,1],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3)

ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

c.set_rasterized(True)

figname = f'{fig_folder}Tomography_situation_GPS_{date}'
plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)

#%% Show positioning in (X,Y) coordinates

# convert Geophones GPS coord into (X,Y) coordinates
X_geoph, Y_geoph = dp.GPS2XY(GPS_geoph[:,:,1],GPS_geoph[:,:,0],Lat0,Long0,param_dict['azimuth'])

set_graphs.set_matplotlib_param('single')
fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)

# plot geophones
norm_acq = colors.Normalize(vmin = 0,vmax = 3)
full_cmap = mpl.colormaps['Greens'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,0.8,256)))

for j in range(3,4):
    current_color = cmap(norm_acq(j))
    if j == 0:
        ax.plot(X_geoph[:,j],Y_geoph[:,j],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3,label = 'Geophone')
    else :
        ax.plot(X_geoph[:,j],Y_geoph[:,j],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3)

ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 15)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 15)
ax.set_aspect(1) # scaling y/x

figname = f'{fig_folder}Tomography_situation_XY_{date}'
plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)


#%% Position fiber on image in (X,Y) coordinates 

# # convert DAS GPS coordinates into (X,Y) coordinates 
# X_DAS,Y_DAS = dp.GPS2XY(DAS_water_height['lat'], DAS_water_height['long'], Lat0, Long0, param_dict['azimuth'])

# # convert Geophones GPS coord into (X,Y) coordinates
# X_geoph, Y_geoph = dp.GPS2XY(GPS_geoph[:,:,1],GPS_geoph[:,:,0],Lat0,Long0,param_dict['azimuth'])

# # convert thickness GPS into (X,Y) coordinates
# X_gps = {}
# Y_gps = {}
# for key_date in gps_results.keys():
#     X_gps[key_date],Y_gps[key_date] = dp.GPS2XY(gps_results[key_date]['latitude'],gps_results[key_date]['longitude'],
#                                                 Lat0, Long0, param_dict['azimuth'])
# MS_date = '0210' 
# X_MS, Y_MS = dp.GPS2XY(MS_gps[MS_date]['latitude'],MS_gps[MS_date]['longitude'],
#                                             Lat0, Long0, param_dict['azimuth'])

# # interpolate new positions along the fiber 
# s  = np.linspace(0,1,600)
# X_interp = (1-s)*X_DAS[0] + s*X_DAS[-1] - X_DAS[0]
# Y_interp = (1-s)*Y_DAS[0] + s*Y_DAS[-1] - Y_DAS[0]

# # create line segments (pairs of consecutive points)
# points = np.array([Y_interp,X_interp]).T.reshape(-1,1,2)
# segments = np.concatenate([points[:-1], points[1:]], axis=1)

# # distance from fiber beginning
# distance = np.sqrt((X_interp)**2 + (Y_interp)**2)

# # create a norm and a line collection 
# norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
# lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 4)
# lc.set_array(distance) # values used for colormap


# set_graphs.set_matplotlib_param('single')
# markersize = 6
# fig, ax = plt.subplots(figsize = (12,6))
# c = ax.pcolormesh(Yreal - Y_DAS[0],Xreal - X_DAS[0],img[:,:,0],shading = 'auto', cmap = 'gray')
# c.set_rasterized(True)
# # ax.plot(X_DAS,Y_DAS,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 5)
# ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
# ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
# ax.set_aspect(1) # set aspect ratio to 1 

# # Plot DAS in color
# # ax.add_collection(lc)

# # Plot DAS without colors
# ax.plot(Y_interp ,X_interp,'k',lw = 2.5,label = 'DAS')

# # plot geophones
# norm_acq = colors.Normalize(vmin = 0,vmax = 3)
# full_cmap = mpl.colormaps['Greens'].resampled(256)
# cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,0.8,256)))

# for j in range(4):
#     current_color = cmap(norm_acq(j))
#     if j == 0:
#         ax.plot(Y_geoph[:,j]-Y_DAS[0],X_geoph[:,j]-X_DAS[0],'^',color = current_color,mec = 'k',
#             ms = 6,zorder = 3,label = 'Geophone')
#     else :
#         ax.plot(Y_geoph[:,j]-Y_DAS[0],X_geoph[:,j]-X_DAS[0],'^',color = current_color,mec = 'k',
#             ms = 6,zorder = 3)
        

# # plot drillings
# norm_h = colors.Normalize(vmin = 0.3,vmax = 0.9)# create colormap for drilling thickness
# cmap = mpl.colormaps['coolwarm'].resampled(256)

# for i,key_date in enumerate(gps_results.keys()):
#     for j in range(len(gps_results[key_date]['h'])):
#         current_color = cmap(norm_h(gps_results[key_date]['h'][j]))
#         if j == 0 and i == 0:
#             scatter = ax.plot(Y_gps[key_date][j]-Y_DAS[0],X_gps[key_date][j]-X_DAS[0],'p',color = current_color,mec = 'k',
#                     ms = 8,zorder = 2,label = 'Drilling')
#         else:
#             scatter = ax.plot(Y_gps[key_date][j]-Y_DAS[0],X_gps[key_date][j]-X_DAS[0],'p',color = current_color,mec = 'k',
#                     ms = 8,zorder = 2)

# for i in range(len(MS_gps[MS_date]['h'])):
#     current_color = cmap(norm_h(MS_gps[MS_date]['h'][i]))
#     ax.plot(Y_MS[i] - Y_DAS[0],X_MS[i] - X_DAS[0],'p',color = current_color,mec = 'k',
#             ms = 8,zorder = 2)
    
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="2%", pad=0.1)

# sm = cm.ScalarMappable(cmap = cmap, norm = norm_h)
# sm.set_array([])  # Only needed for the colorbar
# cbar = fig.colorbar(sm, cax=cax)
# cbar.set_label(r'$h \; \mathrm{(m)}$')

# # cbar = plt.colorbar(scatter,cax = cax)

# # ax.set_xlim(np.array([-225,390]) - Y_DAS[0]) # [-225,380]
# # ax.set_ylim(np.array([-120,120]) - X_DAS[0]) # [-120,120]

# ax.legend(ncols = 3, bbox_to_anchor = (0.11, 1),
#               loc='lower left')

# figname = fig_folder + 'Camera_coord_YX_situation_picture_0205_18_doc_010__with_MS_core_colorbar_h' 


