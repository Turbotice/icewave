# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:07:16 2024

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
import pickle
import cv2 as cv
from datetime import datetime, time , timedelta
import pytz

import icewave.field.drone as drone 
import icewave.field as field 
import icewave.tools.rw_data as rw_data
import icewave.field.drone_projection as dp 

#%% 

date = '0223'
drone_ID = 'Bernache'
flight_ID = '17-waves_012'

tz_dic = {}
tz_dic['Bernache'] = pytz.timezone('America/Montreal')
tz_dic['mesange'] = pytz.timezone('Europe/Paris')

feet2meter = 0.3048
focale = 2700
Lx = 3840
Ly = 2160
t0 = '13:04:06'

#%%

def find_SRT_idx(List,t0,dt_format,tz = pytz.timezone('UTC')):

    relevant_SRT = None
    
    for i in range(len(List)):
        start_t = List[i]['time'][0]
        end_t = List[i]['time'][-1]
        
        # convert time to datetime object 
        start_t = full_date + '-' + start_t
        end_t = full_date + '-' + end_t
        start_t = datetime.strptime(start_t, dt_format)
        end_t = datetime.strptime(end_t, dt_format)
        start_t = start_t.replace(tzinfo = tz)
        end_t = end_t.replace(tzinfo = tz)
        
        # compare t0 to start and end time of each SRT file 
        if t0 >= start_t and t0 < end_t:
            relevant_SRT = i
        
    if relevant_SRT == None :
        print('Selected time t0 does not correspond to any SRT file')
        
    return relevant_SRT


# def get_vertices_GPS(date,drone_ID,flight_ID,t0,Lx = 3840, Ly = 2160,focale = 2700):


#%% Get drone flightrecords 
path2flightrecords = 'K:/Share_hublot/Data/' + date + '/Drones/' + drone_ID + '/flightrecords/'

file2load = path2flightrecords + 'Flightrecord_dict.pkl'
with open(file2load,'rb') as pfile :
    flight_record = pickle.load(pfile)
    
    
#%% Load records

path2records = 'K:/Share_hublot/Data/'+ date + '/Summary/'
file2load = path2records + 'records_' + date + '.pkl'

with open(file2load, 'rb') as pfile :
    records = pickle.load(pfile)
        
#%% Define a t0 and extract relevant parameters 

full_date = flight_record['CUSTOM.date [local]'][0]
t0 = full_date + '-' + t0

t0 = datetime.strptime(t0,'%m/%d/%Y-%H:%M:%S')
t0 = tz_dic[drone_ID].localize(t0)
t0 = t0.astimezone(pytz.timezone('UTC'))

dt_format = '%m/%d/%Y-%I:%M:%S.%f %p'
tz = pytz.timezone('UTC') # timezone used in flight records 
idx_closest_time = dp.closest_time(flight_record['CUSTOM.updateTime [local]'], t0, dt_format, full_date,tz)


h = flight_record['OSD.altitude [ft]'][idx_closest_time]*feet2meter # better to get altitude from .SRT files
alpha_0 = abs(flight_record['GIMBAL.pitch'][idx_closest_time])*np.pi/180
yaw = flight_record['OSD.yaw [360]'][idx_closest_time]*np.pi/180

#%% Find which .SRT file we need to choose
dt_format = '%m/%d/%Y-%H:%M:%S'

relevant_SRT = find_SRT_idx(records['drones'][drone_ID][flight_ID],t0,dt_format)

        
#%% get drone GPS coordinates 

dt_format = '%m/%d/%Y-%H:%M:%S'
tz = pytz.timezone('UTC') # timezone used in flight records 
idx = dp.closest_time(records['drones'][drone_ID][flight_ID][relevant_SRT]['time'], t0, dt_format, full_date,tz)


Lat_drone = records['drones'][drone_ID][flight_ID][relevant_SRT]['latitude'][idx]
Long_drone = records['drones'][drone_ID][flight_ID][relevant_SRT]['longitude'][idx]


#%% 

vertices_real = dp.get_FOV_vertices(Lx, Ly, h, alpha_0, focale)

# get GPS coordinates of camera center
distance_todrone = h/np.tan(alpha_0)
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(Lat_drone, Long_drone, yaw, distance_todrone)

(rho,phi) = dp.cart2pol(vertices_real[:,0,:],vertices_real[:,1,:])
vertices_azimuth = phi*180/np.pi 

vertices_lat,vertices_long = dp.LatLong_coords_from_referencepoint(Lat0, Long0, vertices_azimuth, rho)





