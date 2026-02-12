# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 09:31:41 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob 
import pickle 
import gpxpy
import cartopy.crs as ccrs
import re 
import csv
from datetime import datetime
import pytz

import scipy

import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.weather as weather

parula_map = matcmaps.parula()
plt.rcParams.update({
    "text.usetex": True}) # use latex


#%% Function section 

def dms_to_decimal(dms_str):
    """ Convert DMS string to decimal degrees."""
    parts = dms_str[:-1].split('-')
    degrees = float(parts[0])
    minutes = float(parts[1])
    seconds = float(parts[2])
    decimal = degrees + minutes / 60 + seconds / 3600
    if dms_str[-1] in ['S', 'W']:
        decimal = -decimal
    return decimal

#---------------------------------------------------------------------------

def get_gpx(gpx_file):
    """ Collect gpx file and show its content """
    
    with open(gpx_file,'r') as file:

        gpx = gpxpy.parse(file) 
        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    print('Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation))

        for waypoint in gpx.waypoints:
            print('waypoint {0} -> ({1},{2})'.format(waypoint.name, waypoint.latitude, waypoint.longitude))

        for route in gpx.routes:
            print('Route:')
            
    return gpx


#---------------------------------------------------------------------------

def waypoints_loc(gpx):
    """ Extract longitude and latitude of waypoints stored in a GPX file :
        - Inputs : gpx, a gpxpy.gpx.GPX object
        - Outputs : Long, array of longitude and Lat, array of latitude """
        
    Long,Lat = [],[]
    for waypoint in gpx.waypoints: # loop over all waypoints 
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)
            
        
    return Long, Lat 

#---------------------------------------------------------------------------

def distance_GPS(lat,long,Lat0,Long0,R_earth = 6371e3):
    """ Computes distance between GPS coordinates. 
        Inputs : - lat, numpy array of latitude, in degree
                 - long, numpy array of longitude, in degree
                 - Lat0, reference latitude, in degree
                 - Long0, reference longitude, in degree
        Output : - rho, array of distance between all points of (lat,long) coordinate compared to
        reference point of coordinate (Lat0,Long0) """
        
    rho = R_earth*np.sqrt(((lat - Lat0)*np.pi/180)**2 + np.cos(Lat0*np.pi/180)**2 * ((long - Long0)*np.pi/180)**2)
    
    return rho

# -------------------------------------------------------------------------

def closest_bathymetry(GPS_coords,data_bath):
    """ Get bathymetry for a given GPS coordinate. 
    Inputs: - GPS_coords, tuple, or array, [lat,long] of point for which we look for bathymetry 
            - data_bath, numpy array, dimensions [nb_samples,3], data_bath[:,0] = latitude, data_bath[:,1] = longitude
            data_bath[:,2] = water height, choosing low tide as a reference
    Output : - (lat,long,H), latitude, longitude and water height of the closest point saved in the data base of bathymetry"""
    
    dist = distance_GPS(data_bath[:,0], data_bath[:,1], GPS_coords[0], GPS_coords[1])
    idx_min = np.argmin(dist)
    
    lat,long,H = data_bath[idx_min,:]
    return (lat,long,H)

#--------------------------------------------------------------------------

def get_bathymetry_GPS(GPS_coords,path2interpolator):
    """ Get bathymetry at a given GPS coordinate from linear interpolation 
    Input : - GPS_coords, tuple or array, [lat,long]
            - path2interpolator, path to bathymetry interpolator
    Output : - H, bathymetry at the given GPS position, gives water level at lowest tide"""
    
    with open(path2interpolator,'rb') as pf:
        interp_H = pickle.load(pf)
        
    H = interp_H(GPS_coords[0],GPS_coords[1])
    return H

#%% Import bathymetry data 

path2bath = 'U:/General/Bathymetrie/NONNA_CHS/Bathymetry/'
file2load = f'{path2bath}NONNA10_4830N06890W.txt'

coordinates = []

with open(file2load, 'r') as file:
    next(file)  # Skip header

    for line in file:
        lat_dms, lon_dms, depth_str, *_ = line.strip().split('\t')
        lat = dms_to_decimal(lat_dms)
        lon = dms_to_decimal(lon_dms)
        depth = float(depth_str)
        coordinates.append((lat, lon, depth))

# create numpy array 
data_bath = np.array(coordinates)

#%% Save array in pickle file 
file2save = f'{path2bath}NONNA10_4830N06890W_array.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(data_bath,pf)


#%% Plot bathymetry data

set_graphs.set_matplotlib_param('single')

fig, ax = plt.subplots()
scatter = ax.scatter(data_bath[:,1],data_bath[:,0],s = 10,c = data_bath[:,2],cmap = parula_map)

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(scatter,cax = cax)
cbar.set_label(r'$H \; \mathrm{(m)}$')

#%% Zoom on Haha Bay 

Haha_edges = [-68.8395,-68.8017,48.3355,48.3575]
#keep only points within this area 
mask = np.logical_and(data_bath[:,1] > Haha_edges[0],
                      data_bath[:,1] < Haha_edges[1])
mask = np.logical_and(mask,data_bath[:,0] > Haha_edges[2])
mask = np.logical_and(mask,data_bath[:,0] < Haha_edges[3])

Haha_coordinates = data_bath[mask,:]

fig, ax = plt.subplots()
scatter = ax.scatter(Haha_coordinates[:,1],Haha_coordinates[:,0],s = 10,c = Haha_coordinates[:,2],cmap = parula_map)

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(scatter,cax = cax)
cbar.set_label(r'$H \; \mathrm{(m)}$')

ax.set_xlim([Haha_edges[0], Haha_edges[1]])
ax.set_ylim([Haha_edges[2], Haha_edges[3]])

#%% Collect GPS coordinates from Stephane's GPS 

date = '0211'
path2data = path2records = f'U:/Data/{date}/GPS/'
gpx_filelist = glob.glob(path2data + '*.gpx')

gpx_file = gpx_filelist[1]
# get gpx
gpx = get_gpx(gpx_file)

# get waypoints position
wpts = {}
wpts['name'] = [waypoint.name for waypoint in gpx.waypoints]
Long,Lat = waypoints_loc(gpx)
wpts['lat'] = Lat
wpts['long'] = Long

# keep only waypoints of sources
name_pattern = 'S'
source_wpts = {'name':[],'lat':[],'long':[]}
for i in range(len(wpts['name'])):
    if 'S' in wpts['name'][i]:
        source_wpts['name'].append(wpts['name'][i])
        source_wpts['lat'].append(wpts['lat'][i])
        source_wpts['long'].append(wpts['long'][i])

#%% Superpose sources with bathymetry 

fig, ax = plt.subplots()
scatter = ax.scatter(Haha_coordinates[:,1],Haha_coordinates[:,0],s = 10,c = Haha_coordinates[:,2],cmap = parula_map)
ax.plot(source_wpts['long'],source_wpts['lat'],'r^')

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(scatter,cax = cax)
cbar.set_label(r'$H \; \mathrm{(m)}$')

ax.set_xlim([Haha_edges[0], Haha_edges[1]])
ax.set_ylim([Haha_edges[2], Haha_edges[3]])

#%% Get bathymetry close to a GPS point 

# get bathymetry for first source
lat,long,H = closest_bathymetry((source_wpts['lat'][0],source_wpts['long'][0]), data_bath)
print(H)

#%% Get bathymetry from interpolation

path2interpolate_H = 'U:/General/Bathymetrie/linear_interpolator_bathymetry.pkl'
H = get_bathymetry_GPS((source_wpts['lat'][0],source_wpts['long'][0]), path2interpolate_H)    
    

#%% Collect tide data

path2tides = f'U:/Data/{date}/Marees/'
filelist = glob.glob(f'{path2tides}*array.pkl')

file2load = filelist[0]
with open(file2load,'rb') as pf:
    tide = pickle.load(pf)

# interpolate tide height using datetime 
# create datetime object

test = datetime(2025,2,11,18,42,43,tzinfo = pytz.timezone('UTC'))
tide_height = weather.tide_from_datetime(test)

#%% Interpolate bathymetry 

interp_H = scipy.interpolate.LinearNDInterpolator(Haha_coordinates[:,:2], Haha_coordinates[:,2])
lat = np.linspace(Haha_coordinates[:,0].min(),Haha_coordinates[:,0].max(),100)
long = np.linspace(Haha_coordinates[:,1].min(),Haha_coordinates[:,1].max(),100)

Long,Lat = np.meshgrid(long,lat)
H = interp_H(Lat,Long)

fig, ax = plt.subplots()
ax.pcolormesh(Long,Lat,H,cmap = parula_map)





