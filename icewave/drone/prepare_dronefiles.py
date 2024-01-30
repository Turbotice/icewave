# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:47:06 2024

@author: sebas
"""
import glob
import pylab
import imageio
import cv2

import os


#%%
def sort_files_hour(path2data,base_storage,start_h,end_h,idx_experiment):
    
    """ Sort files from drone DJI Mavic 3 Pro
    
    Moves files from path2data to folders in base_storage depending on the time indicated in the files name
    - start_h : hour at which the experiment started fmt : hhmmss
    - end_h : hour at which the experiment ended fmt : hhmmss
    - idx_experiment : index of the experiment, should be a string 
    
    """

    filelist = os.listdir(path2data)
    
    for f in filelist :
        
        if f.startswith('DJI_'):
            # extract the date from the file 
            file_date = f[4:18]
            file_idx = f[19:23]
            
            year = file_date[0:4]
            day = file_date[4:8]
            hour = file_date[8:]
            
            directory = base_storage + year + '/' + day + '/' + idx_experiment + '/'
            if not os.path.exists(directory):
                os.makedirs(directory)
          
            x = float(hour)

            if (start_h < x) & (x < end_h):
                old_path = os.path.join(path2data,f)
                new_path = os.path.join(directory,f)
                os.replace(old_path,new_path)
                
                print(f"Renamed {old_path} to {new_path}")

def sort_files_indices(path2data,base_storage):
    """ Sort files from drone DJI Mavic 3 Pro
    
    Moves files from path2data to folders in base_storage depending on the index indicated in the files name
    
    """
    
    filelist = os.listdir(path2data)
    
    for f in filelist :
        
        if f.startswith('DJI_'):
            # extract the date from the file 
            file_date = f[4:18]
            file_idx = f[19:23]
            
            year = file_date[0:4]
            day = file_date[4:8]
            hour = file_date[8:]
            
            directory = base_storage + year + '/' + day + '/' + file_idx + '/'
            if not os.path.exists(directory):
                os.makedirs(directory)

            old_path = os.path.join(path2data,f)
            new_path = os.path.join(directory,f)
            os.replace(old_path,new_path)
                
            print(f"Renamed {old_path} to {new_path}")

    
    
#%% Sort files by hour

path2data = 'W:/SagWin2024/Data/Drone/Test 02_10_2023/DJI_001/'
filelist = os.listdir(path2data)
path_storage = 'W:/SagWin2024/Data/Drone/'

start_h = 180000
end_h = 182000
idx_experiment = 'DJI_007'

sort_files_hour(path2data,path_storage,start_h,end_h,idx_experiment)

#%% Sort files by index
path2data = 'W:/SagWin2024/Data/Drone/Test 02_10_2023/DJI_001/'
filelist = os.listdir(path2data)
path_storage = 'W:/SagWin2024/Data/Drone/'

sort_files_indices(path2data,path_storage)


