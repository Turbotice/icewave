# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:54:32 2024

@author: sebas
"""

#%% ------------------------- MODULE IMPORTATION --------------------------

import numpy as np 
import glob
import os
import matplotlib.pyplot as plt
from datetime import datetime, time 
import gpxpy
import pickle
import csv 
import cartopy.crs as ccrs

#%% Extract longitude and latitude of waypoints : 
def waypoints_loc(gpx):
    """ Extract longitude and latitude of waypoints stored in a GPX file :
        - Inputs : gpx, a gpxpy.gpx.GPX object
        - Outputs : Long, array of longitude and Lat, array of latitude """
        
    Long,Lat = [],[]
    for waypoint in gpx.waypoints: # loop over all waypoints 
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)
            
        
    return Long, Lat    
    

#%% Show waypoints / segments and routes stores in a .gpx file 

path2data = 'D:/PC Seb/These PMMH/Arctic_refuge_2024/Data/Test/GPS/'
gpx_filelist = glob.glob(path2data + '*.gpx')

# gpx_file = gpx_filelist[0]
for gpx_file in gpx_filelist:
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

#%%
fig, ax = plt.subplots(figsize = (8,8), dpi = 400)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
long, lat = [],[]
for gpx_file in gpx_filelist:
    with open(gpx_file,'r') as file:
    
        gpx = gpxpy.parse(file) 
        current_long,current_lat = waypoints_loc(gpx)
        long += current_long
        lat += current_lat
        ax.scatter(current_long,current_lat,s = 0.5,transform = ccrs.Geodetic(),)
        plt.show()

# ax = plt.axes(projection=ccrs.Mollweide())














#%% Tasks : 
# - Plot waypoints in longitude/latitude coordinates 
# - Extract labels from datetime for instance or from name 
# - Project waypoints on a plane 



        # BBox = box_data(Long,Lat,scale=1)
        # ext = extent(BBox)
        # t = tmp_connect()
        # fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
        # ax,figs = display_map(ext,t,title=date,ax=ax,width=600)

        # X,Y = [],[]
        # for waypoint in gpx.waypoints:
        #     if True:#int(waypoint.name)>155 and int(waypoint.name)<250:
        # #if int(waypoint.name) in wpts[key]:
        #         x,y = project(waypoint.longitude,waypoint.latitude)
        #         X.append(x)
        #         Y.append(y)
        #         plt.text(x,y-2*10**(-7),int(waypoint.name))
        # ax.plot(X,Y,'bo')

        # if save==True:
        #     savefolder = os.path.dirname(filename)+'/'
        #     graphes.save_figs(figs,savedir=savefolder,suffix='_gpx')
        
        # return ax,figs


