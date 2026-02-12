# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:18:51 2025

@author: Banquise
"""

import icewave.drone.drone_projection as dp
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import h5py
import icewave.tools.matlab2python as m2p


def open_mat(path, file):
    f = h5py.File(path + file,'r') 
    matdata = open(path + file, 'r')
    data = m2p.mat_to_dict(f, f )
    return data


def px_to_real_line(data, xpix_0, ypix_0, d ,theta, size = 1000):
    
    #prend point d'interet x/y pix0 (en pixel) et fait les coordonnées des point à une distance d en m avec un angle theta en radians. renvoit aussi dist_p0 la distance au point d'origine en m
    
    xm_m = data['m']['X'][ypix_0, xpix_0]
    ym_m = data['m']['Y'][ypix_0, xpix_0]

    x_m = np.linspace(xm_m - d * np.sin(theta) ,xm_m + d * np.sin(theta), size)
    y_m = np.linspace(ym_m - d * np.cos(theta) ,ym_m + d * np.cos(theta), size)
    x0 = data['m']['PIXEL']['x0']
    y0 = data['m']['PIXEL']['y0']
    h_drone = data['m']['DRONE']['h_drone']
    alpha_0 = data['m']['DRONE']['alpha_0']
    focale = data['m']['DRONE']['focale']
    
    x_pix, y_pix = dp.projection_pixel_space(x_m, y_m, x0, y0, h_drone, alpha_0, focale)
    
    # print('xpix = ' + str(np.quantile(x_pix, 0.5) / 16))
    # print('ypix = ' + str(np.quantile(y_pix, 0.5) / 16))
    
    Lx = 3840
    Ly = 2160
    mask = (x_pix > 0) & (x_pix < Lx) & (y_pix > 0) & (y_pix < Ly)
    x_pix = x_pix[mask]
    y_pix = y_pix[mask]
    x_m = x_m[mask]
    y_m = y_m[mask]
    
    dist_p0 =   np.sqrt( (x_m - xm_m)**2 + (y_m- ym_m) **2)* np.sign( (y_m- ym_m) +0.0000000001) * np.sign( (x_m - xm_m) + 0.0000000001)
    #* np.sign( (y_m- ym_m) + 0.0000000001)    
    return x_pix, y_pix, dist_p0


def real_to_px(data, x_m, y_m) :
    x0 = data['m']['PIXEL']['x0']
    y0 = data['m']['PIXEL']['y0']
    h_drone = data['m']['DRONE']['h_drone']
    alpha_0 = data['m']['DRONE']['alpha_0']
    focale = data['m']['DRONE']['focale']
    
    x_pix, y_pix = dp.projection_pixel_space(x_m, y_m, x0, y0, h_drone, alpha_0, focale)
    
    # print('xpix = ' + str(np.quantile(x_pix, 0.5) / 16))
    # print('ypix = ' + str(np.quantile(y_pix, 0.5) / 16))
    
    Lx = 3840
    Ly = 2160
    mask = (x_pix > 0) & (x_pix < Lx) & (y_pix > 0) & (y_pix < Ly)
    x_pix = np.asarray(np.round(x_pix[mask]/16), dtype = int)
    y_pix = np.asarray(np.round(y_pix[mask]/16), dtype = int)

    return x_pix, y_pix



def extract_rectangle(data, x_0, y_0, t_0, eta, d_x, d_y, size = 100):
    theta = 90* np.pi / 180
    
    x_pix_x, y_pix_x, dist_p0_x = px_to_real_line(data, x_0, y_0, d_x, theta, size)
    theta = 0
    x_pix_y, y_pix_y, dist_p0_y = px_to_real_line(data, x_0, y_0, d_y, theta, size)
    
    x_pix_x_pix = np.asarray(x_pix_x/16, dtype = int)
    y_pix_y_pix = np.asarray(y_pix_y/16, dtype = int)
            
    rectangle = np.zeros((size,size))
    for i in range(len(x_pix_x_pix)) :
        for j in range(len(y_pix_y_pix)) :
            rectangle[i,j] = eta[x_pix_x_pix[i], y_pix_y_pix[j], t_0]
            
    return dist_p0_x, dist_p0_y, rectangle


def extract_line_Vz(data, Vz, x_pix, y_pix, t) :
    
    #depuis une ligne de pixel, return Vz sur cette ligne en interpollant
    
    field = Vz[:,:,t]
    Fz = RegularGridInterpolator((data['m']['PIXEL']['x_pix'],data['m']['PIXEL']['y_pix']),field)
    Vz_line = Fz( (x_pix,y_pix) )
    
    return Vz_line


def mesure_kappa(forme, x, a, imax) :
    #forme : les données à regarder
    #x : abscisse des pts
    #a : largeur du fit (en pixel)
    #imax : l'indice autour duquel calculer la courbure
    yfit = forme[imax-a:imax+a]
    xfit = x[imax-a:imax+a]
    
    popt_max = np.polyfit(xfit,yfit, 2, full = True)
    
    yth = np.polyval(popt_max[0], xfit)
    
    return np.abs(popt_max[0][0]*2)


def drone_pix2real(x,y,t,param) :# [xreal,yreal,treal] = #,X5)

    # Converts pixel coordinates (x,y,t) to real framework
    # (xreal,yreal,treal). Takes the following arguments : 
    #  - x,y,t : coordinates in pixels, and image numbers of piv frame
    #  (check reference)
    #  - param : parameters structure that correspond to a given drone with
    #  fields
    #   alpha_0 : angle between horizontal and the camera optical axis
    #   (rad)
    #   h : drone height in meter
    #   f : focal length in pixel
    #   fps : frame per second in frame per second
    #   ilag : delay in number of frames
    #   x_0 : center of x-axis camera sensor (pix)
    #   y_0 : center of y-axis camera sensor (pix)

    # For the transformations:
    #   fh : homotethie factor
    #   mirrorX : boolean, mirror operation along x
    #   mirrorY : boolean, mirror operation along y
    #   theta : rotation angle, performed always using the optical center.
    

    x_0 = param['x_0'] #pix
    y_0 = param['y_0']
    h = param['h'] # meter
    alpha_0 = param['alpha_0'] # rad
    f = param['f'] # pix , focal length

    # translation en temps (Fulmar référence)
    ilag = param['ilag']
    treal = (t + ilag)/param['fps'] # real time in sec
    
    # Definition of X and Y in real framework
    
    yreal = (y - y_0)*h/np.sin(alpha_0)/(f*np.sin(alpha_0) + (y - y_0)*np.cos(alpha_0))
    xreal = (x - x_0)*h/(f*np.sin(alpha_0) + (y - y_0)*np.cos(alpha_0))
    
    xreal = xreal
    yreal = -yreal # y-axis upward 

    #homothetie
    fh = param ['fh']
    xreal = fh*xreal
    yreal = fh*yreal
    
    #mirror
    if param['mirrorX'] :
        xreal = -xreal

    
    if param['mirrorY'] :
        yreal = -yreal

    
    #rotation
    theta = param['theta']
    
    xreal = xreal*np.cos(theta) + yreal*np.sin(theta)
    yreal = yreal*np.cos(theta) - xreal*np.sin(theta)
    
    #translation
    #load buoy 5 as reference
    #X_ref = X5(t,:);
    #xreal = xreal - X_ref(:,1);
    #yreal = yreal - X_ref(:,2);
    
    return xreal,yreal,treal

def drone_real2pix(xreal,yreal,treal,param) :

    # Converts pixel coordinates (x,y,t) to real framework
    # (xreal,yreal,treal). Takes the following arguments : 
    #  - x,y,t : coordinates in pixels, and image numbers of piv frame
    #  (check reference)
    #  - param : parameters structure that correspond to a given drone with
    #  fields
    #   alpha_0 : angle between horizontal and the camera optical axis
    #   (rad)
    #   h : drone height in meter
    #   f : focal length in pixel
    #   fps : frame per second in frame per second
    #   ilag : delay in number of frames
    #   x_0 : center of x-axis camera sensor (pix)
    #   y_0 : center of y-axis camera sensor (pix)

    # For the transformations:
    #   fh : homotethie factor
    #   mirrorX : boolean, mirror operation along x
    #   mirrorY : boolean, mirror operation along y
    #   theta : rotation angle, performed always using the optical center.



    x_0 = param['x_0'] #pix
    y_0 = param['y_0']
    h = param['h'] # meter
    alpha_0 = param['alpha_0'] # rad
    f = param['f'] # pix , focal length

    # translation en temps (Fulmar référence)
    ilag = param['ilag']
    t = treal*param['fps'] - ilag # real time in sec


        #rotation
    theta = -param['theta']
    
    xreal = xreal*np.cos(theta) + yreal*np.sin(theta)
    yreal = yreal*np.cos(theta) - xreal*np.sin(theta)

        #mirror
    if param['mirrorX'] :
        xreal = -xreal

    
    if param['mirrorY']:
        yreal = -yreal

    
    #homothetie
    fh = param['fh']
    xreal = xreal/fh
    yreal = yreal/fh

    xreal = xreal
    yreal = -yreal # y-axis upward 

    y = yreal*f*np.sin(alpha_0)/(h/np.sin(alpha_0) - yreal*np.cos(alpha_0)) + y_0
    x = xreal/h*(f*np.sin(alpha_0) + (y - y_0)*np.cos(alpha_0))+x_0


    return x,y,t



