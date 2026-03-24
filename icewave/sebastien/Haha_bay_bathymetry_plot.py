# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:56:24 2026

@author: sebas
"""

import rasterio
from rasterio.fill import fillnodata
from rasterio.plot import show 
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np 
import shapely 
import glob
import cv2 as cv

import icewave.tools.rw_data as rw
import icewave.drone.drone_projection as dp
import icewave.sebastien.set_graphs as set_graphs

full_ocean = mpl.colormaps['ocean_r'].resampled(256)
subset_ocean = mcolors.ListedColormap(full_ocean(np.linspace(0,0.7,256)))

# set path to fig_folder
fig_folder = 'K:/Share_hublot/Data/Summary/SWIIFT_article/'

#%% Load batymetry and coastlines data

base = 'C:/Users/sebas/Matlab_tools/dany_map/'

path2bathy = f'{base}NONNA10_4830N06890W.tiff'
path2coastline = f'{base}estuaire_line.shp'

with rasterio.open(path2bathy) as src:
    bathymetry_data = src.read(1)
    transform = src.transform
    crs = src.crs
    nodata_val = src.nodata   # Get official NoData value
    
coastline = gpd.read_file(path2coastline)

# Check if CRS is the same for .shp file and .tiff file 
if coastline.crs != crs:
        print(f"Reprojection des côtes de {coastline.crs} vers {crs}...")
        coastline = coastline.to_crs(crs)
        
#%% Load UAV data

date = '0226'
base = 'K:/Share_hublot/Data/'

drone_ID = 'mesange'
exp_ID = '12-FRAC_001'
path2drone = f'{base}{date}/Drones/{drone_ID}/matData/{exp_ID}/'
file2load = glob.glob(f'{path2drone}*scaled.h5')[0]

S = rw.load_dict_from_h5(file2load)

#%% Import image 

# Load image 
path2img = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/images/*0684*.tiff'
filelist_img = glob.glob(path2img)
file2load = filelist_img[0]

img = cv.imread(file2load)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

#%% Modify bathymetry data

# build a mask based 
if nodata_val is not None:
        mask_valid = bathymetry_data != nodata_val
else:
    # If NoData is not defined, 
    # we manually set the limit value (here all values below 500)
    mask_valid = bathymetry_data < 500
    
# Interpolate missing values 
print("Interpolation of missing values...")
bathy_interp = fillnodata(bathymetry_data, mask=mask_valid, max_search_distance=100,
                          smoothing_iterations = 1)

#bathymetry_data[bathymetry_data > 10] = None
depth = -1*bathy_interp
# depth = np.flipud(depth)

long = np.arange(src.bounds.left,src.bounds.right,step = 1e-4)
lat = np.arange(src.bounds.bottom + 1e-4,src.bounds.top,step = 1e-4)
lat = np.flip(lat)

longmat,latmat = np.meshgrid(long,lat)

#%% Modify coastline data 

long_min = -68.85
long_max = -68.80
lat_min = 48.33
lat_max = 48.36

# keep coastlines which are within the desired rectangle 
inside = coastline.clip_by_rect(long_min,lat_min,long_max,lat_max)
 
# Create a polygon (rectangle) based on the limits of the map 
map_limits = shapely.geometry.box(long_min,lat_min,long_max,lat_max)

# merge all coastline into a single continuous line
merged_line = inside.geometry.union_all()

# cut map rectangle by the merged line
polygon_cut = shapely.ops.split(map_limits,merged_line)
# create a new geopanda dataframe with the new created polygons
zones_gdf = gpd.GeoDataFrame(geometry=list(polygon_cut.geoms), crs=coastline.crs)

# Sort polygons by area 
zones_gdf['area'] = zones_gdf.geometry.area
zones_gdf = zones_gdf.sort_values(by = 'area',ascending = False).reset_index(drop = True)

coast = zones_gdf.drop(index = 1) # keep only land

#%% Georeference picture
GPS_D = (S['GPS']['latitude'],S['GPS']['longitude'],S['GPS']['azimuth'])
img_lat,img_long = dp.georeference_from_param(img,S['DRONE']['h_drone'],S['DRONE']['alpha_0'],
                                              S['DRONE']['focale'],GPS_D)

img_extent = [img_long.min(), img_long.max(), img_lat.min(),img_lat.max()]
#%% Plot data

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
# show(depth, transform=transform, ax=ax, cmap='ocean_r')
# imsh = ax.images[0]

# cbar = fig.colorbar(imsh, ax=ax, fraction=0.036, pad=0.04)
# cbar.set_label("Profondeur (m)")

# Afficher la bathymétrie
levels_cf = np.arange(-5,25,1)
levels_cl = np.arange(0,25,5)
cf = ax.contourf(longmat,latmat,depth, levels = levels_cf, cmap = subset_ocean,extend = 'both')
cl = ax.contour(longmat,latmat,depth,levels = levels_cl,colors = 'k', linewidths = 0.8)
ax.clabel(cl, inline=True, fontsize = 10, fmt='%1.0f') # add labels on contour lines

# create colorbar
axins = inset_axes(
    ax,
    width="50%",  # width: 50% of parent_bbox width
    height="5%",  # height: 5%
    loc="lower center",
    bbox_to_anchor=(0.1, 0.15, 1, 1),
    bbox_transform = ax.transAxes
)
axins.xaxis.set_ticks_position("bottom")
cbar = fig.colorbar(cf, cax=axins,ticks = levels_cl, orientation="horizontal")

# cbar = fig.colorbar(cf, ax=ax, shrink=0.8)
cbar.set_label('Depth (m)')

# B. Ajouter l'aplat de couleur pour la terre ferme
# On met une couleur (ex: 'whitesmoke', 'lightgrey' ou 'tan')
# zorder=2 force l'affichage de cette couche par-dessus le raster
coast.plot(ax=ax, facecolor='whitesmoke', edgecolor='black', linewidth=1.8, zorder=2)

# Add UAV picture 
#ax.pcolormesh(img_long,img_lat,img[:,:,2],shading = 'auto',cmap = 'gray',zorder = 1)
ax.plot(img_long.mean(),img_lat.mean(),'s',color = 'k')

ax.set_xlim([long_min,long_max])
ax.set_ylim([lat_min,lat_max])
ax.set_aspect(1/np.cos(np.mean(latmat)*np.pi/180)) # scaling y/x
ax.set_xlabel('Longitude (°)')
ax.set_ylabel('Latitude (°)')


figname = f'{fig_folder}Haha_bathymetry_SWIIFT_BicWin24'

plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')





