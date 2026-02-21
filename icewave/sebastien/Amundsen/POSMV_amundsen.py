# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 11:00:58 2024

@author: sebas
"""

""" This script enables to extract several parameters associated to the Amundsen navigation as : 
    heading, heaving, picthing, rolling motion, GPS position 
"""


import numpy as np
import os 
import glob 
from datetime import datetime

import icewave.tools.rw_data as rw_data

def get_datetime(line):
    raw_t = line[1]
    t = raw_t.split('.')[0]
    datetime_str = f'{year}{date}{t}'
    # datetime_obj = datetime.strptime(datetime_str, datetime_format)
    # print(datetime_obj)

    return datetime_str

def get_latlong(line):
    
    lat_str = line[2].split('.')
    lat_deg = lat_str[0][:2]
    lat_min =  lat_str[0][2:] + '.' + lat_str[1]
    lat = float(lat_deg) + float(lat_min)/60

    long_str = line[4].split('.')
    long_deg = long_str[0][:3]
    long_min = long_str[0][3:] + '.' + long_str[1]
    long = float(long_deg) + float(long_min)/60
    
    return lat,long

def get_height(line):
    
    height = float(line[9])
    
    return height

#%% Read POSMV file 

date = '0921'
year = '2024'

path2posmv = 'F:/Amundsen_RA_2024/Data/Amundsen/POSMV/LEG04/'
posmv_pattern = f'{path2posmv}posmv_{year}{date}*.log'

filelist = glob.glob(posmv_pattern)

i = 8
file2load = filelist[i]

#%%

datetime_format = '%Y%m%d%H%M%S'
new_format = ''

data = {}
data['time'] = []
data['lat'] = []
data['long'] = []
data['height'] = []
data['heading'] = []
data['roll'] = []
data['pitch'] = []
data['heave'] = []
data['roll_prec'] = []
data['pitch_prec'] = []
data['head_prec'] = []
data['speed'] = []

lines = rw_data.read_csv(file2load)
j = 1 
line = lines[j]
header = line[0]


print(header)
if header == '$INGGA': # collect vessel position
    data['time'].append(get_datetime(line))
    lat,long = get_latlong(line)
    data['lat'].append(lat)
    data['long'].append(long)
    data['height'].append(get_height(line))

    
#%% get latitude longitude
elif header == '$PASHR':
    heading = float(line[2])
    roll = float(line[4])
    pitch = float(line[5])
    heave = float(line[6])
    roll_accu = float(line[7])
    pitch_accu = float(line[8])
    head_accu = float(line[9])
    
elif header == 'INVTG':
    speed = float(line[7])
    


# lat = float(lines[j][2])/100
# print(lat)

# long = float(lines[j][4])/100
# print(long)









