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
filename_key = re.findall(r'DJI_\s*(\d+_\d+)',filename)
file2load = path + '/' + filename

s = {} # initialize dictionnary 
local_timearea = pytz.timezone('America/Montreal')
UTC_timearea = pytz.timezone('UTC')
year = '2024'
month = '02'
day = '11'

key_UTC = year + '-' + month + '-' + day
print(key_UTC)

#%% Find relevant data 

with open(file2load,'r') as f:
    full_txt = f.read()
    
    s['latitude'] = re.findall(r'\[latitude:\s*(\d+\.\d+)',full_txt) 
    s['longitude'] = re.findall(r'\[longitude:\s*(-\d+\.\d+)',full_txt)
    s['altitude'] = re.findall(r'\[rel_alt:\s*(\d+\.\d+)',full_txt)
    date_regex = re.compile(key_UTC)
    date_pattern = date_regex.pattern + r'\s*(\d+:\d+:\d+.\d+)' 
    local_dates = re.findall(r'\s*(\d+-\d+-\d+ \d+:\d+:\d+.\d+)',full_txt)
    # UTC_dates = re.findall(date_pattern,full_txt) # find dates according to date pattern 

  
#%% 
# Add dates to dictionnary s 
s['UTC_t'] = local_dates
s['UTC_str'] = local_dates

N = np.size(s['latitude'])
for i in range(N):
    s['latitude'][i] = float(s['latitude'][i])
    s['longitude'][i] = float(s['longitude'][i])
    s['altitude'][i] = float(s['altitude'][i])
    
    # convert date string to datetime object
    current_date = local_dates[i] + '000'
    datetime_local = datetime.strptime(current_date, '%Y-%m-%d %H:%M:%S.%f')
    # convert local datetime object to UTC datetime object 
    datetime_local = local_timearea.localize(datetime_local)
    s['UTC_t'][i] = datetime_local.astimezone(UTC_timearea)
    s['UTC_str'][i] = s['UTC_t'][i].strftime("%Y-%m-%d %H:%M:%S.%f")
print('Data structuration Done.')


#%% Save data in a pickle file
pkl_file = path + filename_key[0] + '_drone_geoposition.pkl'
with open (pkl_file,'wb') as f:
    pickle.dump(s,f)

print('Data saved as a pickle file')

#%% Check if data can be loaded
with open (pkl_file,'rb') as f :
    a = pickle.load(f)

print('Data loaded from file : ' + pkl_file)

#%% Write data in a csv file 

csv_file = path + 'drone_geoposition_' + filename.strip('.SRT') + '.csv'
with open(csv_file, 'w') as csvfile: 
    writer = csv.DictWriter(csvfile, fieldnames = s.keys()) 
    writer.writeheader() 
    writer.writerow(s)


    