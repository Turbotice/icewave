# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 09:03:05 2024

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

#%% Function definition 

def projection_real_space(x,y,x0,y0,h,alpha_0,focale):
    """Definition of x and y in real framework, camera sensor center is taken as a reference 
       Inputs : 
        - x: array of x-coordinates in pixels
        - y: array of y-coordinates in pixels
        - x0 : x-coordinate of camera sensor center
        - y0 : y-coordinate of camera sensor center
        - h : drone altitude in meter (above sea level)
        - alpha_0 : inclination angle of the camera, angle to the horizontal (rad) 
        - focale : camera focal length (pixels) """

    yreal = (y - y0)*h/np.sin(alpha_0)/(focale*np.sin(alpha_0) + (y - y0)*np.cos(alpha_0))
    xreal = (x - x0)*h/(focale*np.sin(alpha_0) + (y - y0)*np.cos(alpha_0))
    
    yreal = -yreal
    
    return xreal,yreal

def get_FOV_vertices(Lx,Ly,h,alpha_0,focale):
    """ Return coordinates of the 4 vertices of a given image. Coordinates are given in the local coordinate system of the drone, 
    choosing point corresponding to the camera center as the origin 
    Inputs : 
        - Lx : horizontal dimension of the studied image (pixels)
        - Ly : vertical dimension of the studied image (pixels)
        - h : drone altitude (meters)
        - alpha_0 : camera pitch angle (to the horizontal) (rad)
        - focale : camera focal length (pixels)
        
    Outputs : vertices_real : 4 x 2 numpy array, 
        #dim 1 corresponds to vertices index : bottom-left / top-left / bottom-right / top-right 
        #dim 2 local coordinates (meters) : x,y 
        #dim 3 corresponds to time index, if a sequence of drone height is given 
        
        /!\ Extremities of the FOV are chosen !! And not the bottom left corner of each pixel """

    x0 = Lx/2
    y0 = Ly/2

    # get pixels associated to vertices 
    vertices_pix = np.array([[0 , Ly], [0 , 0], [Lx, Ly], [Lx , 0 ]], dtype = np.float64)
    vertices_pix = np.tile(vertices_pix,(np.size(h),1,1))
    vertices_pix = np.transpose(vertices_pix,(1,2,0))
    vertices_real = np.zeros((4,2,np.size(h)))
    vertices_real[:,0,:],vertices_real[:,1,:] = projection_real_space(vertices_pix[:,0,:],vertices_pix[:,1,:],x0,y0,h,alpha_0,focale)
    
    return vertices_real

def closest_time(strtime_list,t0,dt_format,full_date):
    time_string = full_date + '-' + strtime_list[0]
    time_object = datetime.strptime(time_string,dt_format)
    
    min_delta = abs(t0 - time_object)
    for idx , time_string in enumerate(strtime_list):
        time_string = full_date + '-' + time_string
        time_object = datetime.strptime(time_string, dt_format) 

        delta = abs(t0 - time_object)
        if delta < min_delta:
            min_delta = delta
            idx_min = idx
            
    return idx_min 

def LatLong_coords_from_referencepoint(Lat0,Long0,azimuth,d,R_earth = 6371e3):
    """" Compute GPS coordinates using distance and orientation from a point of known GPS coordinates
    Inputs : 
        - Lat0 : scalar or numpy array, reference point latitude (deg)
        - Long0 : scalar or numpy array, reference point longitude (deg)
        - azimuth : scalar or numpy array, azimuth (deg), defined clockwise from North direction 
        - d : scalar or numpy array, distance to the reference point (meters)"""

    # conversion of angles in rad
    Lat0 = Lat0*np.pi/180
    Long0 = Long0*np.pi/180
    azimuth = azimuth*np.pi/180
    
    Lat = Lat0 + np.cos(azimuth)*d/R_earth
    Long = Long0 + np.sin(azimuth)*d/R_earth/np.cos(Lat0)
    
    # conversion of Latitude and Longitude in deg
    Lat = Lat*180/np.pi
    Long = Long*180/np.pi

    return Lat,Long



def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return (rho,phi)

def pol2cart(rho,phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x,y)
    


#%%

date = '0211'
drone_ID = 'bernache'
flight_ID = '18-stereo_001'


# path to an image 
path2data = 'K:/Share_hublot/PIV_images/'  + date +  '/Drones/' + drone_ID + '/' + flight_ID + '/' + flight_ID +'/'

filelist = os.listdir(path2data)
#%%
idx_frame = 11500
file2load = path2data + filelist[idx_frame]

img = cv.imread(file2load,cv.IMREAD_GRAYSCALE)

fig, ax = plt.subplots()
ax.imshow(img)

#%%
Lx = np.size(img,1)
Ly = np.size(img,0)

# create list of x and y sensor coordinates 
x = np.arange(0,Lx + 1)
y = np.arange(0,Ly + 1)

# X_pix = np.tile(x,(Ly,1))
# Y_pix = np.tile(y,(Lx,1)).T

[X,Y] = np.meshgrid(x,y)

x0 = (x[0] + x[-1])/2
y0 = (y[0] + y[-1])/2

focale = 2700 #camera focale 
h = 89.9 # drone altitude in meter

alpha_0 = 60.1*np.pi/180 # camera pitch angle in rad

#%% Georeferencing 

Xreal,Yreal = projection_real_space(X,Y,x0,y0,h,alpha_0,focale)


#%% test of pcolormesh 

# flip all arrays 
img_flipped = np.flip(img,0)
Xreal_flipped = np.flip(Xreal,0)
Yreal_flipped = np.flip(Yreal,0)


fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal_flipped,Yreal_flipped,img_flipped,shading = 'auto')
fig.colorbar(c,ax = ax)
ax.set_title('Image using flipped arrays ')


#%% Get four corners of the image 
# h = [0.5 , 1 , 2]

# vertices_pix = np.array([[0 , Ly], [0 , 0], [Lx, Ly], [Lx , 0 ]], dtype = np.float64)
# vertices_pix = np.tile(vertices_pix,(np.size(h),1,1))
# vertices_real = np.zeros((np.size(h),4,2))

h = 89.9
vertices_real = get_FOV_vertices(Lx, Ly, h, alpha_0, focale)




#%% Get drone flightrecords 
path2flightrecords = 'K:/Share_hublot/Data/' + date + '/Drones/' + drone_ID + '/flightrecords/'
csv_file = path2flightrecords + 'DJIFlightRecord_2024-02-11_[15-26-43].csv'

record = drone.parse_csv_flightrecord(csv_file,drone = drone_ID)

file2load = path2flightrecords + 'Flightrecord_dict.pkl'
with open(file2load,'rb') as pfile :
    s = pickle.load(pfile)
    
#%% Get Drone orientation at the corresponding image 
full_date = s['CUSTOM.date [local]'][0]

# datetime_local = datetime.strptime(current_date, '%Y-%m-%d %H:%M:%S.%f')

time_string = full_date + '-' + s['CUSTOM.updateTime [local]'][0]
# convert string into datetime object 
time_object = datetime.strptime(time_string, '%d/%m/%Y-%I:%M:%S.%f %p') 

t0 = full_date + '-' + '15:35:11.690'
t0 = datetime.strptime(t0,'%d/%m/%Y-%H:%M:%S.%f')

# find closest time to start of the analysis 

timedelta_min = abs(t0 - time_object) # initialize timedelta
for idx,time_string  in enumerate(s['CUSTOM.updateTime [local]']):
    time_string = full_date + '-' + time_string
    time_object = datetime.strptime(time_string, '%d/%m/%Y-%I:%M:%S.%f %p') 

    delta = abs(t0 - time_object)
    if delta < timedelta_min :
        timedelta_min = delta
        idx_min = idx

dt_format = '%d/%m/%Y-%I:%M:%S.%f %p'
idx_min = closest_time(s['CUSTOM.updateTime [local]'], t0, dt_format, full_date)
#%% Load records for 0211

path2records = 'K:/Share_hublot/Data/'+ date + '/Summary/'
file2load = path2records + 'records_' + date + '.pkl'

with open(file2load, 'rb') as pfile :
    records = pickle.load(pfile)
    


#%% Get Drone GPS position for the selected t0

SRTfile_idx = 2
full_date = records['drones']['Bernache'][flight_ID][SRTfile_idx]['date'][0]
dt_format = '%Y-%m-%d-%H:%M:%S'

t0_UTC = full_date + '-' + '20:35:11'
t0_UTC = datetime.strptime(t0_UTC, dt_format)
idx_GPS = closest_time(records['drones']['Bernache'][flight_ID][SRTfile_idx]['time'], t0_UTC, dt_format, full_date)

Lat_drone = records['drones']['Bernache'][flight_ID][SRTfile_idx]['latitude'][idx_GPS]
Long_drone = records['drones']['Bernache'][flight_ID][SRTfile_idx]['longitude'][idx_GPS]

#%% Compute (Lat,Long) of the point associated to the camera center 

yaw = s['OSD.yaw [360]'][idx_min] # yaw angle in deg
distance_todrone = h/np.tan(alpha_0)
Lat,Long = LatLong_coords_from_referencepoint(Lat_drone, Long_drone, yaw, distance_todrone)

#%% Extract GPS coordinates of each vertices 

(rho,phi) = cart2pol(vertices_real[:,0],vertices_real[:,1])
vertices_azimuth = phi*180/np.pi 

vertices_lat,vertices_long = LatLong_coords_from_referencepoint(Lat, Long, vertices_azimuth, rho)


#%% Extract latitude / longitude of cellphones and buoys 
fig, ax = plt.subplots()

ax.plot(Long,Lat,'k*',label = drone_ID)
ax.plot(Long_drone,Lat_drone,'kd', label = 'Camera center')

phones_coord = np.zeros((len(records['phones']),2))
for i,cellphone in enumerate (records['phones']):
    for exp_ID in records['phones'][cellphone].keys():
        
        phones_coord[i,0] = records['phones'][cellphone][exp_ID]['latitude']
        phones_coord[i,1] = records['phones'][cellphone][exp_ID]['longitude']
        
        ax.plot(float(records['phones'][cellphone][exp_ID]['longitude']),float(records['phones'][cellphone][exp_ID]['latitude']),'bo',
                label = 'phones')


for i,buoy in enumerate(records['buoys']):
    for exp_ID in records['buoys'][buoy].keys():

        ax.plot(records['buoys'][buoy][exp_ID]['longitude'],records['buoys'][buoy][exp_ID]['latitude'],'ro',label = 'buoys')


# ax.plot(records['gps']['garmin_sp']['Canot_02']['longitude'],records['gps']['garmin_sp']['Canot_02']['latitude'],'ks')
ax.plot(records['gps']['garmin_sp']['Canot_01']['longitude'],records['gps']['garmin_sp']['Canot_01']['latitude'],'ks',label = 'Canot_01')

ax.plot(vertices_long,vertices_lat,'go',label = 'camera vertices')
ax.set_xlim([-70.093,-70.09])
ax.set_ylim([48.250,48.256])

#%%

figname = path2records + 'GPS_mutli_instrument_map_' + date 
plt.savefig(figname + '.pdf', dpi = 600, bbox_inches = 'tight')
plt.savefig(figname + '.png', dpi = 600, bbox_inches = 'tight')





#------------------------------------------------------------------------------------------
#%% Get drone GPS position 

path2geoposition = 'H:/Rimouski_2024/Data/2024/0211/Drones/bernache/stereo_001/'
file2load = path2geoposition + '20240211153234_0273_drone_geoposition.pkl'
with open(file2load, 'rb') as pfile :
    d = pickle.load(pfile)


#%% Test drone package functions

date = '0211'

srtfiles = drone.get_srtfiles(date)
drone_records = drone.get_records(date)

    
    