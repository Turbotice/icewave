# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:54:07 2026

@author: sebas
"""
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import os 

import icewave.sebastien.set_graphs as set_graphs
plt.rcParams.update({
    "text.usetex": True}) # use latex

#%%
fig_folder = 'F:/Rimouski_2025/DAS_article/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%%
# 1. Initialize figure and use a suitable projection for Eastern Canada
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots(figsize=(10, 10), 
                       subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-62, 
                                                                       central_latitude=46)})

# 2. Set map limits [xmin, xmax, ymin, ymax] using PlateCarree for the coordinate system
# Ensure longitude and latitude are ordered from minimum to maximum
ax.set_extent([-80, -50.0, 35, 60], crs=ccrs.PlateCarree())

# 3. Add background elements 
ax.add_feature(cfeature.LAND, facecolor='lightgrey', edgecolor='black', zorder=1)
ax.add_feature(cfeature.OCEAN, facecolor='#f0f8ff', zorder=0)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=2)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5, zorder=2) # Added for US/Canada border

# 4. Add gridlines tailored to Eastern Canada
lon_lines = np.arange(-80, -50, 10)
lat_lines = np.arange(40, 60, 10)
# gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--',
#                   xlocs=lon_lines, ylocs=lat_lines)

# 5. Plot Bic National Park (approx: 48.35° N, -68.80° W)
bic_lon, bic_lat = -68.80, 48.35
ax.plot(bic_lon, bic_lat, marker='o', color='red', markersize=8, 
        transform=ccrs.PlateCarree(), zorder=5, label='Bic National Park')

montreal_lon,montreal_lat = -73.56 , 45.50
ax.plot(montreal_lon, montreal_lat, marker='o', color='red', markersize=8, 
        transform=ccrs.PlateCarree(), zorder=5)

ny_lon,ny_lat = -74.00 , 40.71
ax.plot(ny_lon, ny_lat, marker='o', color='red', markersize=8, 
        transform=ccrs.PlateCarree(), zorder=5)


# Add a label next to the dot (optional)
# ax.text(bic_lon + 0.5, bic_lat, 'Bic National Park', 
        # transform=ccrs.PlateCarree(), fontsize=10, weight='bold', zorder=5)

# Save the figure
figname = f'{fig_folder}map_Eastern_Canada'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

plt.show()