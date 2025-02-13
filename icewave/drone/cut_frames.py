# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 17:30:41 2025

@author: sebas
"""

""" 
This script enables to decimate frames from a folder containing all frames associated to a given video. 
"""


#%% Import functions

import os
import glob 
import shutil 


def drone_decimation(base,date,drone_ID,exp_ID,fps,final_fps):
    folder_images = f'{base}/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/'

    folder_name = exp_ID + '_decimation'
    destination_folder = f'{base}/PIV_images/{date}/Drones/{drone_ID}/{folder_name}/'
    frames_decimation(folder_images,destination_folder,fps,final_fps)
    


def frames_decimation(folder_images,destination_folder,fps,final_fps):
    """"Perform decimation of frames associated to a given movie.
        Inputs : 
            - folder_images : folder where images are saved
            - destination_folder : folder where decimated images will be saved
            - fps : initial frame rate of the movie
            - final_fps : final fps you want after decimation
        """
    filelist = glob.glob(folder_images + '*.tiff')
    
    if not os.path.isdir(destination_folder):
        os.mkdir(destination_folder)
    
    if int(fps/final_fps) - fps/final_fps != 0 :
        raise ValueError(f'Initial fps ({fps}) is not a multiple of the final fps ({final_fps})')
    
    for i in range(0,len(filelist),int(fps/final_fps)):
        current_file = filelist[i]
        
        shutil.copy2(current_file,destination_folder)
    
    print(f'Decimation done, images are saved in folder : {destination_folder}')
    
    
if __name__ == '__main__':
    base = 'E:/'
    date = '0208'
    drone_ID = 'mesange'
    exp_ID = '02-waypoints_001'
    
    fps = 30 # intial movie fps
    final_fps = 0.5 # final fps after decimation
    
    drone_decimation(base,date,drone_ID,exp_ID,fps,final_fps)
    