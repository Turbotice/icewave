# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 15:57:24 2025

@author: sebas

Display detailed bathymetry in Ha!Ha! Bay
"""


import os
import numpy as np 
import matplotlib.pyplot as plt 
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

import cmocean.cm as cmo
import alphashape
from shapely.geometry import MultiPoint
from shapely.geometry import Point, Polygon
from scipy.interpolate import griddata
import cartopy.io.shapereader as shpreader


import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps


import icewave.analysis.bathy as bathy
import stephane.display.graphes as graphes
import icewave.gps.gps as gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Set fig_folder

fig_folder = 'U:/Data/0211/DAS/Figures_article/Bathymetry/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Load DAS GPS position 

file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

#%% Load Geophone GPS position 

dates = ['0210','0211','0212']

GPS = {'Geophones':{},'thickness':{}}
date = '0210'
year = '2025'

path2logfile = f'U:/Data/{date}/Geophones'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
UTC_timearea = pytz.timezone('UTC')

# load geophone GPS coordinates 
new_matrix = geophone_gps.get_GPS_coordinates(path2logfile, geophones_table_path, date, year)
    
#%% Compute GPS_logs for each deployment 

GPS_geoph = np.zeros((len(new_matrix.keys()),len(new_matrix['G01'].keys()),2))

for i,geo_key in enumerate(new_matrix.keys()):
    for j,acqu_key in enumerate(new_matrix[geo_key].keys()):
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
            
            GPS_geoph[i,j,0] = GPS_logs['longitude']
            GPS_geoph[i,j,1] = GPS_logs['latitude']
        else:
            
            GPS_geoph[i,j,0] = GPS_logs['longitude'][0]
            GPS_geoph[i,j,1] = GPS_logs['latitude'][0]
                            
        print(GPS_logs)


        
#%% Collect Stephane GPS waypoints 

records = {}
for date in dates :
    record = field_gps.get_records(date,year = '2025')
    records[date] = record['gps']['garminSP']
  
# collect gps thickness measurements 
gps_results = {}
for date in records.keys():
    
    gps_results[date] = {'latitude':[],'longitude':[],'h':[],'time':[],'name':[],'elevation':[]}
    for key in records[date].keys():
        if 'H' in key :
            for key_property in records[date][key].keys():
                gps_results[date][key_property].append(records[date][key][key_property][0])
                
            gps_results[date]['h'].append(int(key[2:])/100)
    
    for key in gps_results[date].keys():
        gps_results[date][key] = np.array(gps_results[date][key])    
        
#%% Maddie SMith GPS waypoints
MS_gps = {}
MS_gps['longitude'] = [-68.814822, -68.819197]
MS_gps['latitude'] = [48.347624, 48.347229]
MS_gps['h'] = [0.31,0.56]

#%% Get HaHa Bay bathymetry

bathy_data = bathy.get_bathy()
# ax,fig = bathy.display_bathy(bathy_data,"haha")

site = 'haha'

if site=='capelans':
    BBox = gps.boxes('capelans')
elif site=='haha':
    BBox = gps.boxes('haha')

b = gps.check_box(bathy_data['Long'],bathy_data['Lat'],BBox)
    
#%% Show bathymetry as scatter plot

Lat0 = np.mean(bathy_data['Lat'][b])

set_graphs.set_matplotlib_param('single')
fig,ax = plt.subplots(figsize = (12,6))

cm = cmo.deep # set ocean deep colormap
sc = ax.scatter(bathy_data['Long'][b],bathy_data['Lat'][b],c=bathy_data['Depth'][b],cmap=cm,marker='o',vmin=-5,vmax=4)

# plot DAS
ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k',lw = 2.5)

# plot geophones
# norm_acq = colors.Normalize(vmin = 0,vmax = 3)
# full_cmap = mpl.colormaps['Greens'].resampled(256)
# cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,0.8,256)))

# for j in range(4):
#     current_color = cmap(norm_acq(j))
#     if j == 0:
#         ax.plot(GPS_geoph[:,j,0],GPS_geoph[:,j,1],'^',color = current_color,mec = 'k',
#             ms = 6,zorder = 3,label = 'Geophone')
#     else :
#         ax.plot(GPS_geoph[:,j,0],GPS_geoph[:,j,1],'^',color = current_color,mec = 'k',
#             ms = 6,zorder = 3)
        

# # plot drillings

# for i,key_date in enumerate(gps_results.keys()):
#     current_color = 'tab:blue'
#     if i == 0:
#         scatter = ax.plot(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],'p',color = current_color,
#                           mec = 'k',ms = 8,zorder = 2,label = 'Drilling')
#     else:
#         scatter = ax.plot(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],'p',color = current_color,mec = 'k',
#                 ms = 8,zorder = 2)

# ax.plot(MS_gps['longitude'],MS_gps['latitude'],'p',color = 'tab:blue',
#                   mec = 'k',ms = 8,zorder = 2)

ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(sc,cax = cax)#,ticks=[0, 15, 30, 45, 60])
cbar.set_label(r'Depth (m)')

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

figname = f'{fig_folder}Haha_bathymetry_DAS_GPS_superposition'
plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)



# =============================================================================
#%% TEST
# =============================================================================

#%% Interpolate depth

# grid resolution (degrees)
resolution = 0.0001  

lon_grid = np.arange(min(bathy_data['Long'][b]), max(bathy_data['Long'][b]), resolution)
lat_grid = np.arange(min(bathy_data['Lat'][b]), max(bathy_data['Lat'][b]), resolution)

grid_lon, grid_lat = np.meshgrid(lon_grid, lat_grid)

# compute closest concave polygon
coords = np.column_stack((bathy_data['Long'][b], bathy_data['Lat'][b]))
alpha = 0.001   # smaller = more concave
poly = alphashape.alphashape(coords, alpha)
polygon_coords = list(poly.exterior.coords)
poly = Polygon(polygon_coords)  # list of (lon,lat)

# compute mask 
mask = np.array([
    poly.contains(Point(lon, lat))
    for lon, lat in zip(grid_lon.ravel(), grid_lat.ravel())
]).reshape(grid_lon.shape)

# interpolation
# Points where you want the interpolation
grid_points = np.column_stack((grid_lon.ravel(), grid_lat.ravel()))

# Interpolate (“linear”, “cubic”, “nearest”)
grid_depth = griddata(coords, bathy_data['Depth'][b], grid_points, method='linear')
grid_depth = grid_depth.reshape(grid_lon.shape)

# Apply mask
grid_depth_masked = np.where(mask, grid_depth, np.nan)

#%% Show pcolormesh

set_graphs.set_matplotlib_param('single')
fig,ax = plt.subplots(figsize=(12,6))
cm = cmo.deep # set ocean deep colormap

imsh = ax.pcolormesh(grid_lon,grid_lat,grid_depth_masked,shading = 'auto', cmap = cm,
                     vmin=-5,vmax=4)
imsh.set_rasterized(True)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(imsh,cax = cax,)#,ticks=[0, 15, 30, 45, 60])
cbar.set_label('Depth (m)')
# cbar.set_clim()

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')


