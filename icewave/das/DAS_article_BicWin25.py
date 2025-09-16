# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:14:33 2025

@author: sebas
"""


import os
import numpy as np 
import matplotlib.pyplot as plt 
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import imageio as iio
import cv2 as cv
import h5py
from scipy.io import loadmat

import icewave.tools.matlab2python as mat2py
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather

#%% Situation picture 

main_path = 'U/Data/'
date = '0205'
drone_ID = 'mesange'
exp_ID = '18-doc_010'

path2img = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/'
filelist = glob.glob(f'{path2img}*.jpg')
file2read = filelist[0]

# fig_folder = f'U:/Data/{date}/Drones/{drone_ID}/Figures/{exp_ID}/'
# if not os.path.isdir(fig_folder):
#     os.mkdir(fig_folder)
    
img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

# img_name = filename.split('.')[0]

#%%

fig, ax = plt.subplots(figsize = (12,9))
ax.imshow(img)


