#%%
import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import time 

import h5py
import os

sys.path.append('..')
#from vasco.cold_room.post_traitement.piv_grenoble.piv_oblique_reconstruction import *
from piv_oblique_reconstruction import *

#%% inputs params camera_40300722
input_params1 = {}
input_params1['date'] = '0513'
input_params1['acq_num'] = 4
input_params1['f_exc'] = 150
input_params1['freq_acq'] = 157.001
input_params1['cam_SN'] = 40300722
input_params1['delta_x'] = 8
input_params1['delta_y'] = 500
input_params1['imwidth'] = 1900
input_params1['imheight'] = 400
input_params1['sensor_width'] = 1920
input_params1['sensor_height'] = 1200
input_params1['pixel_size'] = (11.33/1920)*1e-3
input_params1['focale_mm'] = 16
input_params1['h_camera'] = (35-27.5) * 1e-2
input_params1['alpha_0_deg'] = 6.9
input_params1['alpha_0'] = np.radians(input_params1['alpha_0_deg'])
input_params1['unit_input_velocity'] = 'pxpersec'
input_params1['box_position_information'] = True
input_params1['a'] = 0
input_params1['W'] = 64
input_params1['Dt'] = 50
input_params1['i0'] = 0
input_params1['N'] = 0


input_params1['x0'] = input_params1['sensor_width']/2 - input_params1['delta_x']
input_params1['y0'] = input_params1['sensor_height']/2 - input_params1['delta_y']
input_params1['focale'] = input_params1['focale_mm']*1e-3/input_params1['pixel_size'] # distance focale convertie en nombre de pixels du capteur
input_params1['alpha_0'] = input_params1['alpha_0_deg']*np.pi/180

mat_dict1,t1,Xreal1,Yreal1,Vz1,pc1,Xreal_centers1,Yreal_centers1 = reconstruct_Vz_field_from_pivdata(input_params = input_params1)

#%% inputs params cam 40307970
input_params2 = {}
input_params2['date'] = '0513'
input_params2['acq_num'] = 4
input_params2['f_exc'] = 150
input_params2['freq_acq'] = 157.001
input_params2['cam_SN'] = 40307970
input_params2['delta_x'] = 8
input_params2['delta_y'] = 500
input_params2['imwidth'] = 1900
input_params2['imheight'] = 400
input_params2['sensor_width'] = 1920
input_params2['sensor_height'] = 1200

input_params2['pixel_size'] = (11.33/input_params2['sensor_width'])*1e-3
input_params2['focale_mm'] = 16
input_params2['h_camera'] = (35-27.5) * 1e-2
input_params2['alpha_0_deg'] = 6.9
input_params2['alpha_0'] = np.radians(input_params2['alpha_0_deg'])
input_params2['unit_input_velocity'] = 'pxpersec'
input_params2['box_position_information'] = False
input_params2['a'] = 0
input_params2['W'] = 64
input_params2['Dt'] = 50
input_params2['i0'] = 0
input_params2['N'] = 0

input_params2['x0'] = input_params2['sensor_width']/2 - input_params2['delta_x']
input_params2['y0'] = input_params2['sensor_height']/2 - input_params2['delta_y']
input_params2['focale'] = input_params2['focale_mm']*1e-3/input_params2['pixel_size'] # distance focale convertie en nombre de pixels du capteur
input_params2['alpha_0'] = input_params2['alpha_0_deg']*np.pi/180

mat_dict2,t2,Xreal2,Yreal2,Vz2,pc2,Xreal_centers2,Yreal_centers2 = reconstruct_Vz_field_from_pivdata(input_params = input_params2)
# %%

dist_between_cameras = 0.18
plt.figure()
plt.pcolormesh(Xreal1,Yreal1,Vz1[:,:,0],vmin=np.min(Vz1[:,:,0])/5,vmax=np.max(Vz1[:,:,0])/5)
plt.colorbar()
plt.figure()
plt.pcolormesh(Xreal2,Yreal2,Vz2[:,:,0],vmin=np.min(Vz2[:,:,0])/5,vmax=np.max(Vz2[:,:,0])/5)
plt.colorbar()
plt.show()












# %%
# test de l'échelle en ouvrant une image d'échelle prise par une caméra

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')

import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp

from PIL import Image
#impath = 'R:/Gre25/Data/0512/cameras/ref_matin/Image__2025-05-12__14-06-34.tiff'
impath = 'R:/Gre25/Data/0514/ref/Image__2025-05-14__15-47-53.tiff'
impath = '/media/turbots/GreDisk/Gre25/Data/0514/ref/Image__2025-05-14__15-47-53.tiff'

im = Image.open(impath)
#im.show()

img_array = np.array(im)


#%%
x = np.arange(len(img_array[0,:])) + 0.5 # pixel position of each pixel
y = np.arange(len(img_array[:,0])) + 0.5 # pixel position of each pixel 

X,Y = np.meshgrid(x,y)

"""
x_edges = np.resize(x,x.size + 1)
y_edges = np.resize(y,y.size + 1)
x_edges[-1] = x_edges[-2]
y_edges[-1] = y_edges[-2]
x_edges = x_edges
y_edges = y_edges
"""

#Xedges,Yedges = np.meshgrid(x_edges,y_edges, indexing = 'ij')

Xreal_centers,Yreal_centers = dp.projection_real_space(X, Y, 1920/2, 1200/2, (35-27.5) * 1e-2, 6.9, 2711)


# forme : dp.projection_real_space(Xedges, Yedges, x0, y0, h_camera
          #                              ,alpha_0,focale)

#%%
#%matplotlib qt
plt.figure()
plt.pcolormesh(Xreal_centers,Yreal_centers,img_array)
plt.show()
# %%
