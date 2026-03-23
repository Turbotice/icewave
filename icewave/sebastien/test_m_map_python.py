# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:56:24 2026

@author: sebas
"""

import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np 

import h5py
import icewave.tools.rw_data as rw
import icewave.tools.matlab2python as mat2py

#%% Load batymetry and coastlines data

base = 'C:/Users/sebas/Matlab_tools/dany_map/'

path2bathy = f'{base}NONNA10_4830N06890W.tiff'
path2coastline = f'{base}estuaire_line.shp'

with rasterio.open(path2bathy) as src:
    bathymetry_data = src.read(1)
    transform = src.transform
    crs = src.crs
    
coastline = gpd.read_file(path2coastline)

# Check if CRS is the same for .shp file and .tiff file 
if coastline.crs != crs:
        print(f"Reprojection des côtes de {coastline.crs} vers {crs}...")
        coastline = coastline.to_crs(crs)

#%% Modify bathymetry data

bathymetry_data[bathymetry_data > 10] = None
depth = -1*bathymetry_data

long = np.arange(src.bounds.left,src.bounds.right,step = 1e-4)
lat = np.arange(src.bounds.bottom + 1e-4,src.bounds.top,step = 1e-4)
lat = np.flip(lat)

longmat,latmat = np.meshgrid(long,lat)

#%% Modify coastline data 

long_min = -68.85
long_max = -68.80
lat_min = 48.33
lat_max = 48.36

inside = coastline.clip_by_rect(long_min,lat_min,long_max,lat_max)

#%% Plot data

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

# Afficher la bathymétrie
imsh = ax.contourf(longmat,latmat,depth,transform = ccrs.PlateCarree())
ax.contour(longmat,latmat,depth,colors = 'k',transform = ccrs.PlateCarree())
# imsh = ax.imshow(depth, extent=[src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top],
               # cmap='viridis', transform=ccrs.PlateCarree())

inside.plot(ax=ax, facecolor = 'None', color='k', linewidth=1.5, transform=ccrs.PlateCarree())

cbar = plt.colorbar(imsh, ax=ax, label='Profondeur (m)')

ax.set_xlim([long_min,long_max])
ax.set_ylim([lat_min,lat_max])
imsh.set_clim([-5,15])

ax.set_xlabel('Longitude (°)')
ax.set_ylabel('Latitude (°)')

