
#%% import libraries

import os
import sys
import time
from time import sleep
import  threading
from multiprocessing import Queue
import platform
from PyQt5.QtGui import *

from pypylon import pylon
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import datetime
#if platform.system() =='Windows':
import pyvisa


#%% find cameras and cameras objects

# Create a CameraFactory instance
tl_factory = pylon.TlFactory.GetInstance()

# Choose the devices to use
device_count = tl_factory.EnumerateDevices()
NUM_MAX_CAMERAS = len(list(device_count))

####################################################################

list_camera_SN = ['22458101','23029848'] # choose the camera SN here

####################################################################

#cam_array_temp = pylon.InstantCameraArray(NUM_MAX_CAMERAS)
NUM_CAMERAS = len(list_camera_SN)
cam_array = pylon.InstantCameraArray(NUM_CAMERAS)

cam_info = tl_factory.GetInstance().EnumerateDevices()
selected_cameras_indices = []
for i in range(len(list(cam_info))):
    camera = cam_info[i]
    if camera.GetSerialNumber() in list_camera_SN:
        selected_cameras_indices.append(i)

#print(selected_cameras_indices)

j=0
for i in selected_cameras_indices:
    cam_array[j].Attach(tl_factory.CreateDevice(device_count[i]))
    cam_array[j].Open()
    print('camera with index '+str(j) + ' has seral number : '+ cam_array[j].DeviceSerialNumber())
    cam_array[j].Close()
    j+=1

if len(device_count)==0:
    print('0 camera found')

#%%
##################################################### INPUTS & SETTINGS ######################################################
freq_acq = 15

list_exposure_time = [3000,3000]  # Replace with your desired exposure time in microseconds
list_width =  [2040,696] # Replace with your desired width (px) (doit etre multiple de 16)
list_height = [500,2040]   # Replace with your desired height (px) (peut etre de nimporte quelle taille)
list_offset_x = [8,676]  # Replace with your desired X offset
list_offset_y = [1000,8]  # Replace with your desired Y offset

# choose number of images to grab and the name of the daily folder
numberOfImagesToGrab = 9000 # choisir nombre d'images
max_num_buffer = 9000 # choisir buffer (pour manip fracture 9000 c'est pas mal normalement)


dailyfolder = 'adoeicnoezijcepoaicjezoiceaoijecin' # choisir le daily folder (ou le créer ?)
#############################################################################################################################






list_str_camera_serial_number = []
list_resulting_frame_rates = []

dtn = datetime.datetime.now()
filename = dailyfolder + 'cameras_settings_'+ str(dtn.year)+str(dtn.month).zfill(2)+str(dtn.day).zfill(2)+str(dtn.hour).zfill(2)+str(dtn.minute).zfill(2)+str(dtn.second).zfill(2) +'.txt'

for idx_cam in range(len(list(cam_array))):
    camera = cam_array[idx_cam]
    print(camera.Width.GetValue())
    camera.AcquisitionFrameRate.SetValue(freq_acq)  # Replace with your desired frame rate
    camera.AcquisitionFrameRateEnable.SetValue(True)            
    #print('frame rate : ',camera.ResultingFrameRate.GetValue())
    camera.ExposureTime.SetValue(list_exposure_time[idx_cam])  # Replace with your desired exposure time in microseconds

    camera.Width.SetValue(list_width[idx_cam])  # Replace with your desired width
    camera.Height.SetValue(list_height[idx_cam])  # Replace with your desired height
    camera.OffsetX.SetValue(list_offset_x[idx_cam])  # Replace with your desired X offset
    camera.OffsetY.SetValue(list_offset_y[idx_cam])  # Replace with your desired Y offset

    camera.MaxNumBuffer.SetValue(max_num_buffer)
    str_camera_serial_number = camera.DeviceSerialNumber()
    
    list_str_camera_serial_number.append(str_camera_serial_number)
    #list_freq_acq_effective.append(freq_acq_effective)
    """
    savefolder = acquisition_folder+'camera_'+str_camera_serial_number+'/' + str(f_exc)+'Hz_'+ str(np.round(freq_acq,4))+'Hz/'
    ## create directory to save images
    if save_images:
        if os.path.exists(acquisition_folder)==False:
            os.mkdir(acquisition_folder)
        if os.path.exists(acquisition_folder+'camera_'+str_camera_serial_number+'/')==False:
            os.mkdir(acquisition_folder+'camera_'+str_camera_serial_number+'/')
        if os.path.exists(savefolder)==False:
            os.mkdir(savefolder)
        if os.path.exists(savefolder + 'images/')==False:
            os.mkdir(savefolder + 'images/')
    """
    ## write camera settings in a txt file    

    file = open(filename,'a+')
    file.write('Camera Serial Number : '+ list_camera_SN[idx_cam]+' \n')
    file.write('Resulting Acquisition Frame Rate :'+str(camera.ResultingFrameRate.GetValue())+' \n')
    file.write('Exposure Time : ' + str(list_exposure_time[idx_cam])+' microseconds \n')
    file.write('Width : ' + str(camera.Width.GetValue())+' \n')  # Replace with your desired width
    file.write('Height : '+str(camera.Height.GetValue())+' \n')  # Replace with your desired height
    file.write('OffsetX : '+str(camera.OffsetX.GetValue())+' \n')  # Replace with your desired X offset
    file.write('OffsetY : '+str(camera.OffsetY.GetValue())+' \n')  # Replace with your desired Y offset
    file.write('Max Buffer Size : '+str(camera.MaxNumBuffer.GetValue())+' \n')
    file.write(' \n')
    file.close()

cam_array.Close()
