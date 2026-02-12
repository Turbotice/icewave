# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 09:40:14 2025

@author: sebas

Gathers few functions useful to extract Garmin gps data, without stephane package (sorry stephane)
"""

import gpxpy
import numpy as np

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

#-----------------------------------------------------------------------------------------------------

def waypoints_loc(gpx):
    """ Extract longitude and latitude of waypoints stored in a GPX file :
        - Inputs : gpx, a gpxpy.gpx.GPX object
        - Outputs : Long, array of longitude and Lat, array of latitude """
        
    Long,Lat = [],[]
    for waypoint in gpx.waypoints: # loop over all waypoints 
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)
            
        
    return Long, Lat 

#----------------------------------------------------------------------------------------------------


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

