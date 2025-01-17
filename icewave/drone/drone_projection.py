# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:01:28 2024

@author: sebas
"""

import pylab as plt
import numpy as np 
from datetime import datetime
import pytz
import time
from scipy.interpolate import RegularGridInterpolator
 
import icewave.tools.rw_data as rw_data
import icewave.tools.datafolders as df



global feet2meter
feet2meter = 0.3048

def get_projected_image(drone,key,num,frame):
    #only work currently for first image of the movie
    
    #Load flight record
    base = df.find_path()
    filename =  base+f'{date}/drones/{drone}/flightrecords/Flightrecord_dict.pkl'
    flight = rw_data.load_pkl(filename)
    print(flight.keys())

    #Load global record
    filename =  base+f'{date}/Summary/records_{date}.pkl'
    records = rw_data.load_pkl(filename)
    record = records['drones'][drone][key][num]

    #Cut flight record
    flight_p = drone.cut_flightrecord(record,flight)

    #Load image
    im = get_exemple_image(record)

    #Find coordinates (Lat,Lon) of all pixels 
    Lats,Lons = project_image(record,flight_p,im,0)

    return Lats,Lons,im

def get_exemple_image(record):
    base = df.find_path()
    filename = base+record['path']+'_exemple.tiff'
    im = plt.imread(filename)
    return im

def project_image(record,flight,im,frame,focale=2700):
    # warning ! only work for frame=0 now.
    # sample frequency is different from mp4, srt and flightrecord movies
    frameF = frame
    frameR = frame
    #focale is set to its default value
    
    h = flight['OSD.altitude [ft]'][frameF]*feet2meter # better to get altitude from .SRT file ?
    alpha_0 = abs(flight['GIMBAL.pitch'][frameF])*np.pi/180
    yaw = flight['OSD.yaw [360]'][frameF]#*np.pi/180

    print(h,alpha_0,yaw)
    Lat_drone = record['latitude'][frameR]
    Long_drone = record['longitude'][frameR]

    # get GPS coordinates of camera center
    distance_todrone = h/np.tan(alpha_0)
    Lat0,Long0 = LatLong_coords_from_referencepoint(Lat_drone, Long_drone, yaw, distance_todrone)

    Ly,Lx,nc = im.shape
    x = np.linspace(0,Lx-1,Lx)
    y = np.linspace(0,Ly-1,Ly)
    [x,y]=np.meshgrid(x,y)
    x0 = Lx/2
    y0 = Ly/2

    X,Y = projection_real_space(x,y,x0,y0,h,alpha_0,focale)
    [ds,azimuths] = cart2pol(X,Y)
    azimuths = 90-azimuths*180/np.pi

    Lats,Lons = LatLong_coords_from_referencepoint(Lat0, Long0, azimuths,ds)
    return Lats,Lons

def projection_real_space(x,y,x0,y0,h,alpha_0,focale):
    """Definition of x and y in real framework, camera sensor center is taken as a reference 
       Inputs : 
        - x: array of x-coordinates in pixels
        - y: array of y-coordinates in pixels
        - x0 : x-coordinate of camera sensor center
        - y0 : y-coordinate of camera sensor center
        - h : drone altitude in meter (above sea level)
        - alpha_0 : inclination angle of the camera, angle to the horizontal (rad) 
        - focale : camera focal length (pixels)
        Outputs :
        - xreal : array of x-coordinates of the positions to the center of the image (in meter)
        - yreal : array of y-coordinates of the positions to the center of the image (in meter)
        """

    yreal = (y - y0)*h/np.sin(alpha_0)/(focale*np.sin(alpha_0) + (y - y0)*np.cos(alpha_0))
    xreal = (x - x0)*h/(focale*np.sin(alpha_0) + (y - y0)*np.cos(alpha_0))
    
    yreal = -yreal
    
    return xreal,yreal



def projection_pixel_space(xreal,yreal,x_0,y_0,h,alpha_0,f):

    # Definition of x and y in pixel framework

    # Inputs : 
    # - xreal: array of x-coordinates in metric system
    # - yreal: array of y-coordinates in metric system
    # - x_0 : x-coordinate of camera sensor center (pixel system)
    # - y_0 : y-coordinate of camera sensor center (pixel system)
    # - h : drone altitude in meter (above sea level)
    # - alpha_0 : inclination angle of the camera, angle to the horizontal 
    # - f : camera focal length

    xreal = xreal
    yreal = -yreal
    
    y = yreal*f*np.sin(alpha_0)/(h/np.sin(alpha_0) - yreal*np.cos(alpha_0)) + y_0
    x = xreal/h*(f*np.sin(alpha_0) + (y - y_0)*np.cos(alpha_0))+x_0

    return x, y

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
    vertices_pix = np.array([[0 , Ly], [0 , 0],[Lx , 0 ], [Lx, Ly]], dtype = np.float64)
    vertices_pix = np.tile(vertices_pix,(np.size(h),1,1))
    vertices_pix = np.transpose(vertices_pix,(1,2,0))
    vertices_real = np.zeros((4,2,np.size(h)))
    vertices_real[:,0,:],vertices_real[:,1,:] = projection_real_space(vertices_pix[:,0,:],vertices_pix[:,1,:],x0,y0,h,alpha_0,focale)
    
    return vertices_real

def closest_time(strtime_list,t0,dt_format,full_date,tz = pytz.timezone('UTC')):
    time_string = full_date + '-' + strtime_list[0]
    time_object = datetime.strptime(time_string,dt_format)
    time_object = time_object.replace(tzinfo = tz)
    
    min_delta = abs(t0 - time_object)
    for idx , time_string in enumerate(strtime_list):
        time_string = full_date + '-' + time_string
        time_object = datetime.strptime(time_string, dt_format) 
        time_object = time_object.replace(tzinfo = tz)

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

def backward_projection(fun,points,y0,h,alpha_0,focale,fps,Dt):
    """ Computation of vertical velocity field from the pixel displacement along vertical of the camera. 
    
    INPUTS : - fun : function of a tuple of pixel coordinates (x,y,t), or a list of tuples. It can typically be an interpolator of 
                pixel vertical displacement 
             - points : list of tuples. Each tuple corresponds to pixel coordinates (x,y,t) /!\ must keep this order of dimensions
             - y0 : float, y-coordinate of the middle of the camera sensor
             - alpha_0 : float, angle (in rad) of the camera axis to the horizontal 
             - focale : float, camera focal length (pixels)
             - fps : float, frame rate used (frame/s)
             - Dt : float, time step between two image compared to computer pixel displacement field 
             
    OUTPUT : - Fp : function of a tuple of pixel coordinates (x,y,t), or a list of tuples, which computes the vertical velocity
                associated to a given tuple
    """
    
    Fp = h*focale*fun(points)/((points[:,1] - y0)*np.cos(alpha_0) + 
                              focale*np.sin(alpha_0))/((points[:,1] - y0 + fun(points))*np.sin(alpha_0) - focale*np.cos(alpha_0))*fps/Dt
    
    return Fp 


def vertical_velocity_from_pixel_displacement(Vy,x_pix,y_pix,t,y0,h,alpha_0,focale,fps,Dt):
    """ Compute vertical velocity field from velocity field of pixel displacement on the camera sensor. 
    
    INPUTS : - Vy : numpy array, velocity field of pixel displacement along the vertical of the camera sensor (y-axis)
             - x_pix : numpy 1D array, x-coordinate of boxes center used to compute Vy
             - y_pix : numpy 1D array, y-coordinate of boxes center used to compute Vy
             - t : numpy 1D array, frame index array 
             - y0 : float, y-coordinate of the middle of the camera sensor
             - alpha_0 : float, angle (in rad) of the camera axis to the horizontal 
             - focale : float, camera focal length (pixels)
             - fps : float, frame rate used (frame/s)
             - Dt : float, time step between two image compared to computer pixel displacement field 
             
    OUTPUT : - Vz : numpy array, vertical velocity field (scaled in m/s), same shape as Vy
    
    """
    
    start_time = time.perf_counter()
    
    # Create the meshgrid of new points
    new_x, new_y, new_t = np.meshgrid(x_pix, y_pix, t, indexing='ij')
    
    if new_x.shape != Vy.shape :
        raise TypeError(f"Vy has shape {Vy.shape} while x, y and t have respective shapes : {x_pix.shape}, {y_pix.shape} and {t.shape}")
    
    # compute interpolation of pixel displacement vertical velocity field 
    Fy = RegularGridInterpolator((x_pix,y_pix,t), Vy)
    
    # Combine the new coordinates into a single list of points to interpolate
    new_points = np.array([new_x.ravel(), new_y.ravel(), new_t.ravel()]).T
    
    # Perform the interpolation
    Vz = backward_projection(Fy,new_points,y0,h,alpha_0,focale,fps,Dt)
    
    # Reshape the interpolated field back to the new grid shape
    Vz = Vz.reshape(len(x_pix), len(y_pix), len(t))
    
    end_time = time.perf_counter()
    # show computation time 
    elapsed_time = end_time - start_time
    print('Elapsed time : ',elapsed_time, ' s')
    
    return Vz
