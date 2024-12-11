# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:35:16 2024

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt 
import os 
import glob

import pandas as pd

import folium
import geopandas as gpd
from shapely.geometry import Point 
import contextily as ctx
import pyproj
#%%

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Présentations/Environnement_Arctic_Refuge/'
img_quality = 800 #image quality in dpi 
#%% 

# set of latitudes longitudes points that I would like to emphasize on the map 
my_dict = {}
key_sites = ['Haha','Saguenay','AR20240909','AR20240914']
for key in key_sites:
    my_dict[key] = {}
    
my_dict['Haha']['lat'] = 48.350538
my_dict['Haha']['long'] = -68.808609

my_dict['Saguenay']['lat'] = 48.251353
my_dict['Saguenay']['long'] = -70.124327

my_dict['AR20240909']['lat'] = 80.8260668
my_dict['AR20240909']['long'] = -66.7653705

my_dict['AR20240914']['lat'] = 79.8690422
my_dict['AR20240914']['long'] = -69.9443962


latitudes = np.zeros(len(my_dict))
longitudes = np.zeros(len(my_dict))
coordinates = []
for i,key_site in enumerate(my_dict.keys()):

    latitudes[i] = my_dict[key_site]['lat']
    longitudes[i] = my_dict[key_site]['long']
    coordinates.append((my_dict[key_site]['lat'],my_dict[key_site]['long']))
#%% Create a map using folium 

# compute center of the map 
map_center = [sum(x[0] for x in coordinates) / len(coordinates), 
              sum(x[1] for x in coordinates) / len(coordinates)]

m = folium.Map(location=map_center, tiles = 'Esri.WorldImagery',zoom_start=3)

# add a marker for each GPS coordinate 
for coord in coordinates:
    folium.Marker(location=coord).add_to(m)


figfile = fig_folder + 'map_JMC.html'
m.save(figfile)

#%% Create a map using geopandas 
# Define a custom latitude/longitude range
lat_min, lat_max = 30, 85  # Custom latitude range (e.g., between 30°N and 55°N)
lon_min, lon_max = -100, +20  # Custom longitude range (e.g., between 130°W and 10°E)

# Convert the list of coordinates to a GeoPandas GeoDataFrame
geometry = [Point(long,lat) for lat, long in coordinates]
gdf = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")

# Reproject to web mercator (EPSG:3857) for contextily
gdf_web_mercator = gdf.to_crs(epsg=3857)


# Create a projection object for conversion
proj_4326_to_3857 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
proj_3857_to_4326 = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True).transform

# Convert custom lat/lon bounds to EPSG:3857 for map display
x_min, y_min = proj_4326_to_3857(lon_min, lat_min)
x_max, y_max = proj_4326_to_3857(lon_max, lat_max)


# Create a plot
fig, ax = plt.subplots(figsize=(10, 10))
gdf_web_mercator.plot(ax=ax, marker='o', color='red', markersize=60, label='Locations', edgecolor = 'k')

# Set custom x and y axis limits (bounding box in EPSG:3857)
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Add a satellite base map using Contextily
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery)


# Create custom ticks for latitude and longitude
x_ticks = np.linspace(x_min, x_max, 5)  # X-axis (Longitude)
y_ticks = np.linspace(y_min, y_max, 5)  # Y-axis (Latitude)

# Convert the ticks from EPSG:3857 back to EPSG:4326 for labeling
x_tick_labels = [proj_3857_to_4326(x, ax.get_ylim()[0])[0] for x in x_ticks]  # Longitude labels
y_tick_labels = [proj_3857_to_4326(ax.get_xlim()[0], y)[1] for y in y_ticks]  # Latitude labels

# Set the new ticks and labels for x and y axes
ax.set_xticks(x_ticks)
ax.set_xticklabels([f"{lon:.2f}°" for lon in x_tick_labels])
ax.set_yticks(y_ticks)
ax.set_yticklabels([f"{lat:.2f}°" for lat in y_tick_labels])

# Save figure 
figfile = fig_folder + 'map_JMC_geophone_deployments'
print('Saving figure..')
plt.savefig(figfile + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figfile + '.png', dpi = img_quality, bbox_inches = 'tight')
print('Figure saved !')


#%% Create a map locally close to Saguenay river 

my_dict = {}
key_sites = {'Haha','Saguenay0210','Saguenay0211','StMathieu'}

for key in key_sites:
    my_dict[key] = {}
    
my_dict['Haha']['lat'] = 48.350538
my_dict['Haha']['long'] = -68.808609

my_dict['Saguenay0211']['lat'] = 48.251353
my_dict['Saguenay0211']['long'] = -70.124327

my_dict['Saguenay0210']['lat'] = 48.364657
my_dict['Saguenay0210']['long'] = -70.720788

my_dict['StMathieu']['lat'] = 48.157260
my_dict['StMathieu']['long'] = -69.017456

latitudes = np.zeros(len(my_dict))
longitudes = np.zeros(len(my_dict))
coordinates = []
for i,key_site in enumerate(my_dict.keys()):

    latitudes[i] = my_dict[key_site]['lat']
    longitudes[i] = my_dict[key_site]['long']
    coordinates.append((my_dict[key_site]['lat'],my_dict[key_site]['long']))
    
    
#%%


# Define a custom latitude/longitude range
lat_min, lat_max = 46, 50  # Custom latitude range (e.g., between 30°N and 55°N)
lon_min, lon_max = -73, -66  # Custom longitude range (e.g., between 130°W and 10°E)

# Convert the list of coordinates to a GeoPandas GeoDataFrame
geometry = [Point(long,lat) for lat, long in coordinates]
gdf = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")

# Reproject to web mercator (EPSG:3857) for contextily
gdf_web_mercator = gdf.to_crs(epsg=3857)


# Create a projection object for conversion
proj_4326_to_3857 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
proj_3857_to_4326 = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True).transform

# Convert custom lat/lon bounds to EPSG:3857 for map display
x_min, y_min = proj_4326_to_3857(lon_min, lat_min)
x_max, y_max = proj_4326_to_3857(lon_max, lat_max)


# Create a plot
fig, ax = plt.subplots(figsize=(10, 10))
gdf_web_mercator.plot(ax=ax, marker='o', color='red', markersize=60, label='Locations', edgecolor = 'k')

# Set custom x and y axis limits (bounding box in EPSG:3857)
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Add a satellite base map using Contextily
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery)


# Create custom ticks for latitude and longitude
x_ticks = np.linspace(x_min, x_max, 5)  # X-axis (Longitude)
y_ticks = np.linspace(y_min, y_max, 5)  # Y-axis (Latitude)

# Convert the ticks from EPSG:3857 back to EPSG:4326 for labeling
x_tick_labels = [proj_3857_to_4326(x, ax.get_ylim()[0])[0] for x in x_ticks]  # Longitude labels
y_tick_labels = [proj_3857_to_4326(ax.get_xlim()[0], y)[1] for y in y_ticks]  # Latitude labels

# Set the new ticks and labels for x and y axes
ax.set_xticks(x_ticks)
ax.set_xticklabels([f"{lon:.2f}°" for lon in x_tick_labels])
ax.set_yticks(y_ticks)
ax.set_yticklabels([f"{lat:.2f}°" for lat in y_tick_labels])


figfile = fig_folder + 'map_JMC_geophone_deployments_BicWin2024'
print('Saving figure..')
plt.savefig(figfile + '.pdf', dpi = img_quality, bbox_inches = 'tight')
plt.savefig(figfile + '.png', dpi = img_quality, bbox_inches = 'tight')
print('Figure saved !')