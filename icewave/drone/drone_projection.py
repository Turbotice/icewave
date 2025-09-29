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

#----------------------------------------------------------------------------------------

def get_exemple_image(record):
    base = df.find_path()
    filename = base+record['path']+'_exemple.tiff'
    im = plt.imread(filename)
    return im

#-----------------------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------------------------------

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

#---------------------------------------------------------------------------------------------

def projection_pixel_space(xreal,yreal,x_0,y_0,h,alpha_0,f):

    """ Definition of x and y in pixel framework

    # Inputs : 
    # - xreal: array of x-coordinates in metric system
    # - yreal: array of y-coordinates in metric system
    # - x_0 : x-coordinate of camera sensor center (pixel system)
    # - y_0 : y-coordinate of camera sensor center (pixel system)
    # - h : drone altitude in meter (above sea level)
    # - alpha_0 : inclination angle of the camera, angle to the horizontal 
    # - f : camera focal length """

    xreal = xreal
    yreal = -yreal
    
    y = yreal*f*np.sin(alpha_0)/(h/np.sin(alpha_0) - yreal*np.cos(alpha_0)) + y_0
    x = xreal/h*(f*np.sin(alpha_0) + (y - y_0)*np.cos(alpha_0))+x_0

    return x, y

#----------------------------------------------------------------------------------------------------------------------------

def georectify_image(img,H,alpha_0,focale):
    """Georectify image taken from a camera with a focal f, 
    at a height H from the filmed plane, with an angle alpha_0 to the horizontal 
    Inputs : - img, numpy.ndarray : image [nx,ny,nc]
             - H, float : height of the camera to the filmed horizontal plane (in meter)
             - alpha_0, float : angle to the horizontal (in rad)
             - focale, float : camera focal length 
    Outputs : - Xreal,Yreal, numpy.ndarray meshgrids : coordinates of each pixel edges to be used by the function plt.pcolormesh
    """
    
    [ny,nx,nc] = np.shape(img) 

    x_edges = np.arange(0,nx + 1)
    y_edges = np.arange(0,ny + 1)

    x0 = (nx + 1) / 2
    y0 = (ny + 1) / 2

    Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

    # compute real coordinates for each pixels of the image 
    Xreal,Yreal = projection_real_space(Xedges,Yedges,x0,y0,H,alpha_0,focale)
    
    return Xreal,Yreal
    
#------------------------------------------------------------------------------------------------------------------------------

# def get_angles(xpix,ypix,x0,y0,focal):
#     """ Compute angles with which a pixel (xpix,ypix) is seen. Upper left corner of the image is chosen as the origin 
#     of the camera sensor coordinate system
#     Inputs : - xpix, float or numpy array, x-coordinate on camera sensor 
#              - ypix, float or numpy array, y-coordinate on camera sensor 
#              - x0, float, x-coordinate of camera sensor center 
#              - y0, float, y-coordinate of camera sensor center 
#              - focal, float, camera focal length (in pixels) 
    
#     Outputs : - beta_x, float or numpy array, angle, with respect to optical axis going throught center of 
#     camera sensor, along x-axis
#              - beta_y, float or numpy array, angle, with respect to optical axis going throught center of 
#     camera sensor, along y-axis
#     - beta_r, float or numpy array, angle, with respect to optical axis going through center of 
#     camera sensor without any projection""" 
            
#     beta_x = np.arcsin((xpix - x0)/focal)
#     beta_y = np.arcsin((ypix - y0)/focal)
#     beta_r = np.arcsin(np.sqrt((xpix - x0)**2 + (ypix - y0)**2)/focal)

#     return beta_x,beta_y,beta_r

#----------------------------------------------------------------------------------------------------------------------------

def get_theta_phi(xpix,ypix,x0,y0,focal):
    """ Compute angles with which a pixel (xpix,ypix) is seen. Upper left corner of the image is chosen as the origin 
    of the camera sensor coordinate system
    Inputs : - xpix, float or numpy array, x-coordinate on camera sensor 
             - ypix, float or numpy array, y-coordinate on camera sensor 
             - x0, float, x-coordinate of camera sensor center 
             - y0, float, y-coordinate of camera sensor center 
             - focal, float, camera focal length (in pixels) 
    
    Outputs : - theta, float or numpy array, angle, with respect to optical axis going throught center of 
    camera sensor without any projection 
              - phi, float or numpy array, angle with respect to the line passing through the camera sensor center 
              and the point of coordinate (0,1) (cf np.arctan2)""" 
              
    theta = np.arctan(np.sqrt((xpix - x0)**2 + (ypix - y0)**2)/focal)
    phi = np.arctan2((ypix-y0),(xpix-x0))

    return theta,phi




#--------------------------------------------------------------------------------------------------------------------------

def get_FOV_vertices(Lx,Ly,h,alpha_0,focale):
    """ Return coordinates of the 4 vertices of a given image. Coordinates are given in the local coordinate system 
    of the drone, 
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

#-----------------------------------------------------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------------------------------------------

def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return (rho,phi)

#------------------------------------------------------------------------------------------------------------------------------

def pol2cart(rho,phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x,y)

#------------------------------------------------------------------------------------------------------------------------

def XY2GPS(X,Y,Lat0,Long0,azimuth):
    """ Convert cartesian coordinates X,Y to GPS coordinates
    Inputs : - X,Y : numpy array, containing cartesian (X,Y) coordinates with respect to a reference point O 
             - Lat0,Long0 : GPS coordinates of reference point (in degrees)
             - azimuth : orientation of the cartesian system with respect to the geographic north (in degrees, between 
                                                                                                   0 and 360)
             This azimuth angle is measured between geographic north direction and axis Y of cartesian system
    Outputs : - Lat,Long : numpy array, GPS coordinates of the selected points"""
    
    # convert cartesian to polar coordinates
    rho,theta = cart2pol(X,Y)
    theta = theta*180/np.pi
    # compute azimuth of all selected points
    local_azimuth = azimuth + 90 - theta
    Lat,Long = LatLong_coords_from_referencepoint(Lat0,Long0,
                                                    local_azimuth,rho)
    return Lat,Long

#-----------------------------------------------------------------------------------------------------------------------------

def GPS2XY(lat,long,Lat0,Long0,azimuth,R_earth = 6371e3):
    """Compute cartesian coordinates (X,Y) from GPS coordinates
    Inputs : - lat,long : numpy array, GPS coordinates (in degrees)
             - Lat0,Long0 : GPS coordinates of the point chosen as a reference for the cartesian coordinates system (in degrees)
             - azimuth : angle with respect to geopgraphic north with which the cartesian system will be oriented (in degrees,
                                                                                                            between 0 and 360)
             This azimuth angle is measured between geographic north direction and axis Y of cartesian system
    Outputs : - X,Y : numpy array, cartesian coordinates (in meter) """
    
    rho = R_earth*np.sqrt(((lat - Lat0)*np.pi/180)**2 + np.cos(Lat0*np.pi/180)**2 * ((long - Long0)*np.pi/180)**2)

    psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(long - Long0)*np.pi/180,(lat - Lat0)*np.pi/180)*180/np.pi

    theta = azimuth - psi + 90
    # compute euclidean coordinates
    X,Y = pol2cart(rho,theta*np.pi/180)
    
    return X,Y

def get_height_from_record(record,indice=0,typ='elev'):
    if typ=='elev':
        s = record['params'][indice].split('rel_alt: ')[1].split(' ')[0]
    return float(s)

def georeference(key,record,image,flight):
    """ Project an image into a georeferenced coordinate system, by using the flight meta data stored
        in the record file (position, height, and pitch), and the flight file (drone orientation).
        Currently, only work for pitch = 90°
    
    INPUTS : - key : identification key common to record, flight and image 
             - record : dictionnary with parameters extracted from the data folder
             - image : image to be projected
             - flight : flight record of the drone
             
    OUTPUT : - Lat
             - Long
             - im
    """
    if type(record[key])==list:
        rec=record[key][0]
    else:
        rec = record[key]
        
    im = image[key]
    H = get_height_from_record(rec)
    print(f'Hauteur : {H}m')
    alpha_0 = np.pi/2
    focale = 2700 #sama focale for everybody
    im = image[key]
    Xr,Yr = georectify_image(im,H,alpha_0,focale)

    Lat0=rec['latitude'][0]
    Long0=rec['longitude'][0]

    azimuth = flight[key]['OSD.yaw [360]']
    print(f'Azimuth : {azimuth}°')
    pitch = float(flight[key]['GIMBAL.pitch'])
    print(f'Pitch : {pitch}°')

    if pitch<-80:
        Lat,Long = XY2GPS(Xr,Yr,Lat0,Long0,azimuth)
    else:
        print(pitch)
        return None,None,im
    
    return Lat,Long,im
#-------------------------------------------------------------------------------------------------------------------------------

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
                              focale*np.sin(alpha_0))/((points[:,1] - y0 + fun(points))*np.sin(alpha_0) 
                                                       - focale*np.cos(alpha_0))*fps/Dt
    
    return Fp 

#------------------------------------------------------------------------------------------------------------------------------

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


#--------------------------------------------------------------------------------------------------------------------------------

def change_XY_reference_system(X,Y,param_projected,param_ref,R,translation,scaling):
    """ Return metric coordinates (X_proj,Y_proj) in a new reference system associated
    to an other drone (ref_drone). Metric coordinates are converted to GPS coordinates, then rotation, translation
    and scaling operations are perfomed (these operations are based on tracking of buoys). 
    
    Inputs : - X,Y : numpy array, metric coordinates
             - param_projected : dictionnary, contain parameters of the drone to be projected
             - param_ref : dictionnary, containe parameters of the drone chosen as a reference
             These parameters must contain following keys : 
                     + h : drone height
                     + alpha_0 : pitch angle (radian)
                     + latitude : drone latitude (deg)
                     + longitude : drone longitude (deg)
                     + azimuth : drone azimuth (deg), between 0 and 360
                 
             - R : numpy array, rotation matrix shape (2,2)
             - translation : numpy array, translation vector (Lat,Long), shape (2,)
             - scaling : scalar, scaling factor 
    Outputs : - X_proj,Y_proj : numpy array, metric coodrinates in the new reference system """
    
    # get GPS coordinates 
    # coordinates of central point in projected_drone coordinates system
    dist2drone = param_projected['h']/np.tan(param_projected['alpha_0'])
    LatD = param_projected['latitude']
    LongD = param_projected['longitude']
    azimuth_drone = param_projected['azimuth']
    Lat0,Long0 = LatLong_coords_from_referencepoint(LatD,LongD,
                                                    azimuth_drone,dist2drone)

    Lat,Long = XY2GPS(X,Y,Lat0,Long0,azimuth_drone)
    latlong = np.stack([Lat,Long],axis = -1)

    # reshape matrix
    reshape_latlong = np.reshape(latlong,(-1,2))

    # apply procruste operations to Lat,Long
    rotated = R @ reshape_latlong.T
    latlong_transfo = scaling * rotated.T + translation[None,...]

    # reshape matrix 
    latlong_transfo = np.reshape(latlong_transfo,(X.shape[0],X.shape[1],2))

    # get GPS coordinates of reference drone
    # coordinates of central point in ref_drone coordinates system
    dist2drone = param_ref['h']/np.tan(param_ref['alpha_0'])
    Lat_ref = param_ref['latitude']
    Long_ref = param_ref['longitude']
    azimuth_ref = param_ref['azimuth']
    Lat0,Long0 = LatLong_coords_from_referencepoint(Lat_ref,Long_ref,
                                                    azimuth_ref,dist2drone)

    X_proj,Y_proj = GPS2XY(latlong_transfo[:,:,0],latlong_transfo[:,:,1],Lat0,Long0,azimuth_ref)
    
    return X_proj,Y_proj    

#--------------------------------------------------------------------------------------------------------------------

def get_uz_from_Vy(Vy,theta,phi):
    """ Compute vertical velocity profile from apparent velocity field computed from Digital Image Correlation. 
    This projection only works if waves are parallel to the coordinate system (x,y).
    Inputs : - Vy, array like, dimensions [nx,ny,nt], apparent velocity field along vertical axis of camera sensor
             - theta, array like, dimensions [ny,nx], angle, with respect to optical axis going through center of 
    camera sensor without any projection 
              - phi, array like, dimensions [ny,nx], angle with respect to the line passing through the camera sensor center 
              and the point of coordinate (x,y) = (1,0) (cf np.arctan2) 
    Output : - uz, array like, dimensions [nx,nt], vertical velocity profile """

    uz = np.zeros((Vy.shape[0],Vy.shape[2])) # create array for uz
    err_uz = np.zeros((Vy.shape[0],Vy.shape[2]))
    for frame in range(uz.shape[1]):
        field = Vy[:,:,frame]
        for idx_x in range(uz.shape[0]):
            f_thetaphi = np.tan(theta[:,idx_x])*np.sin(phi[:,idx_x])
            p,V = np.polyfit(f_thetaphi,field[idx_x,:],1,cov = True)
            err_coeffs = np.sqrt(np.diag(V))
            uz[idx_x,frame] = p[0]
            err_uz[idx_x,frame] = err_coeffs[0]

    return uz, err_uz

def get_ux_from_VxVy(Vx,Vy,theta,phi):
    """ Compute horizontal velocity profile from apparent velocity field computed from Digital Image Correlation. 
    This projection only works if waves are parallel to the coordinate system (x,y).
    Inputs : - Vx, array like, dimensions [nx,ny,nt], apparent velocity field along horizontal axis of camera sensor
             - Vy, array like, dimensions [nx,ny,nt], apparent velocity field along vertical axis of camera sensor
             - theta, array like, dimensions [ny,nx], angle, with respect to optical axis going through center of 
    camera sensor without any projection 
              - phi, array like, dimensions [ny,nx], angle with respect to the line passing through the camera sensor center 
              and the point of coordinate (x,y) = (1,0) (cf np.arctan2) 
    Output : - ux, array like, dimensions [nx,ny,nt], horizontal velocity field. Velocity is supposed to be parallel to sea ice plane """
    
    uz, err_uz = get_uz_from_Vy(Vy,theta,phi)
    uz_xyt = np.transpose(uz[:,np.newaxis,:],(2,0,1))*np.tan(theta.T)*np.cos(phi.T)
    uz_xyt = np.transpose(uz_xyt,(1,2,0))
    ux = Vx - uz_xyt
    
    return ux

def get_uz_ux_from_structure(Vx,Vy,S):
    """ Compute vertical velocity profile as well as horizontal velocity field (uz,ux) from apparent velocity field obtained from 
    Digital Image Correlation (Vx,Vy).
    Inputs : - S, dictionnary, structure returned by matlab code main_data_structuration.m. 
             - Vx, array like, dimensions [nx,ny,nt], apparent velocity field along horizontal axis of camera sensor
             - Vy, array like, dimensions [nx,ny,nt], apparent velocity field along vertical axis of camera sensor
    Outputs : - uz, array like, dimensions [nx,nt], real scaled vertical velocity profile. Vertical velocity is supposed to be 
    perpendicular to sea ice plane 
              - ux, array like, dimensions [nx,ny,nt], horizontal velocity field. Velocity is supposed to be parallel to sea ice plane """
    
    xpix,ypix = np.meshgrid(S['PIXEL']['x_pix'],S['PIXEL']['y_pix'],indexing = 'xy')
    theta,phi = get_theta_phi(xpix,ypix,S['PIXEL']['x0'], S['PIXEL']['y0'],S['DRONE']['focale'])
    
    uz,err_uz = get_uz_from_Vy(Vy,theta,phi)
    ux = get_ux_from_VxVy(Vx,Vy,theta,phi)
    
    return uz,ux,err_uz

