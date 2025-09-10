

import math
import numpy as np

def sph2cart(r,theta,phi):
    #takes list rthetaphi (single coord)
    x = r * sin( theta ) * cos( phi )
    y = r * sin( theta ) * sin( phi )
    z = r * cos( theta )
    return [x,y,z]

def cart2sph(x,y,z):
    #takes list xyz (single coord)
    r       =  np.sqrt(x**2 + y**2 + z**2)
    theta   =  np.arccos(z/r)*180/ np.pi #to degrees
    phi     =  np.arctan2(y,x)*180/ np.pi
    return [r,theta,phi]

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
