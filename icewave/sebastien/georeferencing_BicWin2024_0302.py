# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:02:31 2024

@author: sebas
"""


import os
import numpy as np 
import matplotlib.pyplot as plt 
import pickle
import cv2 as cv
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio

import icewave.field.drone as drone 
import icewave.field as field 
import icewave.tools.rw_data as rw_data
import icewave.drone.drone_projection as dp 
import icewave.tools.datafolders as df

#%%

# import an image from a single day 
date = '0302'
year = '2024'
drone_ID = 'bernache'
flight_ID = '6-Calib_oblique_002'

base = df.find_path('Hublot24')
base_images = base.replace('Data','PIV_images')
path2tiff = f'{base_images}/{date}/Drones/{drone_ID}/{flight_ID}/'

film_list = os.listdir(path2tiff)
film_ID = film_list[0]
filelist_tiff = glob.glob(f'{path2tiff}/{film_ID}/*.tiff')

#%%

frame = filelist_tiff[0]
# load the image 

img = cv.imread(frame)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
# img = np.transpose(img,axes = [1,0,2])

#%% show image 
fig, ax = plt.subplots()
ax.imshow(img)

#%% 

focale = 2700;
h = 11.800; # drone height in meters 
alpha_0 = 17.7*np.pi/180; # drone pitch angle to horizontal 

[ny,nx,nt] = np.shape(img) 

x_edges = np.arange(1,nx + 2)
y_edges = np.arange(1,ny+ 2)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

# define limits of the image (pixel indices)
x_borders = [1,nx] 
y_borders = [500,ny]

mask_xedges = np.where((x_edges >= x_borders[0]) & (x_edges <= x_borders[1] + 1))[0]
mask_yedges = np.where((y_edges >= y_borders[0]) & (y_edges <= y_borders[1] + 1))[0]

x = x_edges[mask_xedges]
y = y_edges[mask_yedges]


sub_img = img[y_borders[0] - 1 : y_borders[1] , x_borders[0] - 1 : x_borders[1] ,:]
sub_img = np.transpose(sub_img,[1,0,2])

Xedges,Yedges = np.meshgrid(x,y,indexing = 'ij')

# compute real coordinates 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,h,alpha_0,focale)


#%%

fig, ax = plt.subplots()
c = ax.pcolormesh(Xedges,np.flip(Yedges,1),sub_img[:,:,0],shading = 'auto')
fig.colorbar(c,ax = ax)

#%% Display georeferenced image 

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,sub_img[:,:,0],shading = 'auto', cmap = 'viridis')
fig.colorbar(c,ax = ax)
ax.set_xlabel(r'$x \quad \mathrm{(m)}$')
ax.set_ylabel(r'$y \quad \mathrm{(m)}$')

#%% 








