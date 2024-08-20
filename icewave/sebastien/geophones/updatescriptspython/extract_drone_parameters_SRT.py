# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 09:49:39 2024

@author: sebas


This script decripts how to extract parameters of a flight from a .SRT file.
Relevant parameters we are interested in are : 
    - UTC_t0 (beginning of recording)
    - Total number of recorded frames 
    - latitude and longitude of drone with time 
    - initial position of drone
    - relative altitude with time 
    - initial relative altitude 
"""


#%% Import module 

import numpy as np
import os
import csv
import matplotlib.pyplot as plt 
import pickle
from datetime import datetime, date, time 
import pytz
import re

#%% Load SRT file 

path = 'F:/Rimouski_2024/Data/2024/0211/Drones/bernache/stereo_001/'
filename = 'DJI_20240211152847_0272_D.SRT'
file2load = path + '/' + filename


timearea = pytz.timezone('America/Montreal')
s = {} # initialize dictionnary 
year = '2024'
month = '02'
day = '11'

key_UTC = year + '-' + month + '-' + day
print(key_UTC)
#%%
with open(file2load,'r') as f:
    full_txt = f.read()
    # lines = f.readlines()
    # char_local_t0 = lines[3].split('\n')
    # char_local_t0 = char_local_t0[0] + '000'
    # local_t0 = datetime.strptime(char_local_t0,'%Y-%m-%d %H:%M:%S.%f') #local time of video start
    # local_t0 = timearea.localize(local_t0) # set timezone of extracted time to the corresponding timearea 
    # print("Local t0 is : " , local_t0)
    
    # UTC_t0 = local_t0.astimezone(pytz.timezone('UTC'))
    # UTC_t0_str = UTC_t0.strftime("%Y-%m-%d %H:%M:%S.%f")
    # print('UTC t0 is :', UTC_t0_str)
    
    # line_data = lines[4]
    # data = line_data.split('[')
    
    
    
    
    latitudes = re.findall(r'\[latitude:\s*(\d+\.\d+)',full_txt) 
    longitudes = re.findall(r'\[longitude:\s*(-\d+\.\d+)',full_txt)
    altitudes = re.findall(r'\[rel_alt:\s*(\d+\.\d+)',full_txt)
    UTC_dates = re.findall(r'2024-02-11\s*(\d+:\d+:\d+.\d+)',full_txt)
    # I need to convert each date to a datetime format
    
    # Next steps : 
        # save data in a pickle file 
        # write data in a .txt file 
    
    
    # s['UTC_t0'] = 
  
#%% 
  
with open(file2load,'r') as f:  
    txt = f.read(1000)
    latitudes = re.findall(r'\[latitude:\s*(\d+\.\d+)',txt)
    longitudes = re.findall(r'\[longitude:\s*(-\d+\.\d+)',full_txt)
    print(longitudes)
    