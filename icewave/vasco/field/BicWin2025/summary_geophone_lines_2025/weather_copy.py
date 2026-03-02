# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 11:40:48 2025

@author: sebas

Gathers all useful functions concerning field weather (tide, bathymetry)
"""

import numpy as np
import pickle

import csv
from datetime import datetime
import pytz
import glob 
import scipy
import sys

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave')
import icewave.tools.datafolders as df


def get_bathy_interpolator(disk = 'Backup25',year = '2025'):
    """ Load bathymetry interpolator and return it"""
    
    base = df.find_path(disk,year)
    path2bathy = f'{base}Bathymetrie/linear_interpolator_bathymetry.pkl'
    
    with open(path2bathy,'rb') as pf:
        interp_H = pickle.load(pf)
        
    return interp_H
    
#-------------------------------------------------------------------------------------------

def get_bathymetry_GPS(GPS_coords,interp_H):
    """ Get bathymetry at a given GPS coordinate from linear interpolation 
    Input : - GPS_coords, tuple or array, [lat,long]
            - interp_H, bathymetry interpolator
    Output : - H, bathymetry at the given GPS position, gives water level at lowest tide"""
        
    H = interp_H(GPS_coords[0],GPS_coords[1])
    return H

#-------------------------------------------------------------------------------------------

def tide_from_datetime(UTC_datetime,disk = 'Backup25',year = '2025'):
    """ Compute tide height from a datetime object based on UTC schedule """
    monthday = UTC_datetime.strftime('%m%d')
    
    directory = df.find_path(disk = disk,year = year )
    path2tides = f'{directory}/{monthday}/Marees/'
    filelist = glob.glob(f'{path2tides}*tides_array.pkl')
    
    file2load = filelist[0]
    with open(file2load,'rb') as pf:
        tide = pickle.load(pf)
    
    tide_interp = scipy.interpolate.interp1d(tide['timestamp_UTC'], tide['tide_height'])
    tide_height = tide_interp(UTC_datetime.timestamp())
    
    return tide_height

#--------------------------------------------------------------------------------------------------------

def get_water_height(GPS_coords,UTC_datetime,interp_bathy = None,disk = 'Backup25', year = '2025'):
    """ Compute water height from bathymetry data and tide height.
    Input : - GPS_coords, tuple or array, (lat,long)
            - UTC_datetime, datetime object, UTC based 
    Output : - water_height, water depth in meter"""
    
    if interp_bathy == None:    
        interp_bathy = get_bathy_interpolator(disk,year)
    
    H = get_bathymetry_GPS(GPS_coords, interp_bathy)
    tide_height = tide_from_datetime(UTC_datetime,disk,year)
    
    water_height = H + tide_height
    
    return water_height
    

#--------------------------------------------------------------------------------------

def get_tidefromcsvfile(file_tide,format_datetime = '%Y-%m-%dT%H:%M:%S.%f%z'):
    """ Create a tide dictionnary from a csv file of tides
    Inputs : - file_tide, string, path to a csv file containing two headers : 'date' and 'tide_height'
             - format_datetime, string, format with which date is saved in the csv file
    Output : - tide, dictionnary, with keys : + date : string containing local date
                                              + tide_height : height of tide 
                                              + datetime_UTC : UTC datetime object """
                                              
    tide = {'date':[],'tide_height':[]}
    with open(file_tide,'r') as file:
        csv_file = csv.DictReader(file)
        for line in csv_file:
            tide['date'].append(line['date'])
            tide['tide_height'].append(float(line['tide_height']))

    # Convert a date in a datetime object
    UTC_timearea = pytz.timezone('UTC')

    tide['datetime_UTC'] = []
    for d in tide['date']:
        localdatetime = datetime.strptime(d, format_datetime)
        UTC_time = localdatetime.astimezone(UTC_timearea)
        tide['datetime_UTC'].append(UTC_time)
        
    tide['timestamp_UTC'] = [UTC_time.timestamp() for UTC_time in tide['datetime_UTC']]
    return tide 



#--------------------------------------------------------------------------------------------------------

def generate_tide_array(disk = 'Backup25',year = '2025'):
    """ Generate tide array for all available tide data. Tide data must be saved as a csv file, using suffixe :
        'tides.csv' 
        Inputs : - disk, string, optional, details disk on which data are saved 
                 - year, string, optional, details year """
                 
    directory = df.find_path(disk,year)
    print(directory)
    
    filelist = glob.glob(directory + '**/Marees/*tides.csv',recursive = True)

    for file2tide in filelist:
        tide = get_tidefromcsvfile(file2tide)
    
        file2save = file2tide.replace('.csv','_array.pkl')
        print(file2save)
        with open(file2save,'wb') as pf: 
            pickle.dump(tide,pf)
            
            
    
if __name__ == '__main__':
    generate_tide_array(disk = 'Backup25')

