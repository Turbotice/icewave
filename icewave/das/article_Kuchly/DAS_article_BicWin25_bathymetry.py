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

# module for map plotting
import cmocean.cm as cmo
import alphashape
from shapely.geometry import MultiPoint
from shapely.geometry import Point, Polygon
from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader

# icewave modules
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps


import icewave.analysis.bathy as bathy
# import stephane.display.graphes as graphes
import icewave.gps.gps as gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Set fig_folder

fig_folder = 'U:/Data/Summary/DAS/Bathymetry/'
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

# gps.display_haha(ax)
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
# plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
# plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
# plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)

#%% Try Chatgpt method

# load data 
lon = bathy_data['Long'][b]
lat = bathy_data['Lat'][b]
depth = bathy_data['Depth'][b]

# create a regular grid for colormap
xi = np.linspace(lon.min(), lon.max(), 300)
yi = np.linspace(lat.min(), lat.max(), 300)
Xi, Yi = np.meshgrid(xi, yi)

Zi = griddata((lon, lat), depth, (Xi, Yi), method='cubic')


proj = ccrs.PlateCarree()
# get shorelines
land = cfeature.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='black',
    facecolor='lightyellow',
)

# coastline = cfeature.NaturalEarthFeature('physical', 'coastline', '10m',
#                                          edgecolor = 'black')

shp = shpreader.natural_earth(resolution='10m',
                              category='physical',
                              name='coastline')

#%% Create figure 

fig, ax = plt.subplots(
    figsize=(10, 7),
    subplot_kw=dict(projection=proj)
)

# Set extent (lon_min, lon_max, lat_min, lat_max)
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=proj)

# Plot interpolated bathymetry
im = ax.pcolormesh(
    Xi, Yi, Zi,
    cmap='viridis',
    shading='auto',
    transform=proj
)

# Add contours
cs = ax.contour(
    Xi, Yi, Zi,
    levels=[-1,0,1,2,5],
    colors='black',
    linewidths=0.7,
    transform=proj
)
ax.clabel(cs, inline=True, fontsize=8, fmt='%dm')

# Add land/shorelines

# Add GSHHG coastlines (high res)
ax.add_feature(
    cfeature.GSHHSFeature(scale='full'),
    facecolor='lightyellow',
    edgecolor='black',
    linewidth=1.0
)

# ax.add_feature(land,zorder = 0)
# ax.add_feature(coastline,zorder = 0)
# ax.coastlines

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')


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


#%%


# 1. Configuration de la zone (Region du Bic / Baie du Ha! Ha!)
# Coordonnées approximatives lues sur votre carte
lat_min, lat_max = 48.30, 48.38  # ~ 48°18' à 48°23'
lon_min, lon_max = -68.90, -68.80 # ~ 68°54' à 68°48'

# 2. Création de la figure avec une projection
# La projection Mercator est souvent utilisée pour les cartes marines, 
# mais Lambert Conformal est mieux pour le Canada à grande échelle.
# Ici, PlateCarree suffit pour une petite zone, ou Mercator pour le look "carte marine".
projection = ccrs.Mercator()
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': projection})

# Définir l'étendue de la carte
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# --- SIMULATION DE VOS DONNÉES DE BATHYMÉTRIE ---
# (Remplacez ceci par votre code d'interpolation / pcolormesh)
# Je crée un gradient simple pour l'exemple
x = np.linspace(lon_min, lon_max, 100)
y = np.linspace(lat_min, lat_max, 100)
X, Y = np.meshgrid(x, y)
Z = -1 * (np.sin((X+68.9)*20) * np.cos((Y-48.3)*20) * 30) # Fausse bathymétrie
# Affichage de la bathymétrie
cmap = plt.cm.viridis_r  # Palette similaire à votre image
bathy = ax.pcolormesh(X, Y, Z, transform=ccrs.PlateCarree(), cmap=cmap, shading='auto', vmin=0, vmax=25)
# ------------------------------------------------

# 3. Ajout de la LIGNE DE CÔTE HAUTE PRÉCISION
# 'scale' peut être: 'c' (crude), 'l' (low), 'i' (intermediate), 'h' (high), 'f' (full)
# Pour votre échelle (quelques km), utilisez 'f' (full) si installé, sinon 'h'.
high_res_land = cfeature.GSHHSFeature(scale='f', levels=[1], facecolor='lightgray', edgecolor='black')

# Ajout des terres. IMPORTANT : zorder élevé pour cacher l'interpolation sous la terre.
ax.add_feature(high_res_land, zorder=10) 

# # 4. Ajout de la grille (Gridlines)
# gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linestyle='--')
# gl.top_labels = False
# gl.right_labels = False
# gl.xlabel_style = {'size': 10}
# gl.ylabel_style = {'size': 10}

# 5. Barre de couleur (Colorbar)
cbar = plt.colorbar(bathy, ax=ax, orientation='horizontal', pad=0.05, shrink=0.5)
cbar.set_label('Profondeur (m)')

# Titre
plt.title("Bathymétrie - Parc national du Bic (Baie du Ha! Ha!)")

plt.show()
