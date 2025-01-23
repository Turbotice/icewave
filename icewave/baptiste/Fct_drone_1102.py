# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:18:51 2025

@author: Banquise
"""

import icewave.drone.drone_projection as dp
from scipy.interpolate import RegularGridInterpolator
import numpy as np


def px_to_real_line(data, xpix_0, ypix_0, d ,theta):
    
    #prend point d'interet x/y pix0 (en pixel) et fait les coordonnées des point à une distance d en m avec un angle theta en radians. renvoit aussi dist_p0 la distance au point d'origine en m
    
    xm_m = data['m']['X'][ypix_0, xpix_0]
    ym_m = data['m']['Y'][ypix_0, xpix_0]

    x_m = np.linspace(xm_m - d * np.sin(theta) ,xm_m + d * np.sin(theta), 1000)
    y_m = np.linspace(ym_m - d * np.cos(theta) ,ym_m + d * np.cos(theta), 1000)
    x0 = data['m']['PIXEL']['x0']
    y0 = data['m']['PIXEL']['y0']
    h_drone = data['m']['DRONE']['h_drone']
    alpha_0 = data['m']['DRONE']['alpha_0']
    focale = data['m']['DRONE']['focale']
    
    x_pix, y_pix = dp.projection_pixel_space(x_m, y_m, x0, y0, h_drone, alpha_0, focale)
    
    print('xpix = ' + str(np.quantile(x_pix, 0.5) / 16))
    print('ypix = ' + str(np.quantile(y_pix, 0.5) / 16))
    
    Lx = 3840
    Ly = 2160
    mask = (x_pix > 0) & (x_pix < Lx) & (y_pix > 0) & (y_pix < Ly)
    x_pix = x_pix[mask]
    y_pix = y_pix[mask]
    x_m = x_m[mask]
    y_m = y_m[mask]
    
    dist_p0 =  np.sign( (x_m - xm_m) + 0.0000000001)  * np.sqrt( (x_m - xm_m)**2 + (y_m- ym_m) **2)
    #* np.sign( (y_m- ym_m) + 0.0000000001)
    
    
    return x_pix, y_pix, dist_p0


def extract_line_Vz(data, Vz, x_pix, y_pix, t) :
    
    #depuis une ligne de pixel, return Vz sur cette ligne en interpollant
    
    field = Vz[:,:,t]
    Fz = RegularGridInterpolator((data['m']['PIXEL']['x_pix'],data['m']['PIXEL']['y_pix']),field)
    Vz_line = Fz( (x_pix,y_pix) )
    
    return Vz_line