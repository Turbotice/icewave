# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 11:28:05 2024

@author: sebas
"""

import numpy as np 
import pandas as pd
import csv
import glob
import os 
import time
import datetime as dt


#%% FUNCTION SECTION

def convert_date(w):
    """ convert a date 'mm/jj/yyyy' to a format : 'yyyymmjj' """
    print(w)
    w_length = len(w)
    year = w[w_length - 4:]
    day = w[w_length - 7 : w_length - 5]

    if w_length == 10 :
        month = w[:2]
    else :
        month = '0' + w[0]
    f = year + month + day
     
    return f

def convert_time(w):
    
    """ convert a time '3:20:40.50 PM' to a format : 'hhmmssddd' """

    a = w[-2:]
    h = w[0]
    m = w[2:4]
    s = w[5:7]
    millis = w[8:10] + '0'

    if a == 'PM':

        hour = int(h) + 12
        h = str(hour)

    
    f = h + m + s + millis
    
    return f 

def convert_UTC(date,t):
    
    """ Merge date 'yyyymmjj' and a time 'hhmmssddd' to UTC time 'yyyymmjj'T'hhmmssddd'Z format, used by Dany Dumont team"""
    
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:8])

    hour = int(t[:2])
    minute = int(t[2:4])
    sec = int(t[4:6])
    millisec = t[6:]
    print(millisec)
    # convert a time year, month, day, hour, min, sec, microsec to a timme in seconds since the epoch
    t = dt.datetime(year,month,day,hour,minute,sec).timestamp()

    # convert a time to UTC time 
    UTC_t = time.gmtime(t)
    
    y_txt = str(UTC_t[0])
    UTC_txt = y_txt
    
    for i in range(1,6):
        
        a_txt = str(UTC_t[i])
        if UTC_t[i] < 10:
            a_txt = '0' + a_txt 
        
        UTC_txt += a_txt     
        
        if i == 2:
            UTC_txt += 'T'
        
    UTC_txt += millisec + 'Z'
    
    return UTC_txt
    

def get_timeline_row(data,idx,ref,objet,comments,facq):
    """ This function takes the index of a given panda DataFrame, and builds a row to be written in a timeline csv file"
    The timeline has the following headers :
    'Index','Objet','T_0','T_F','X','Y','Z','Latitude','Longitude','Elevation','Commentaire','Facq','Theta'
    
    The function takes as arguments : 
    - data, a pandas DataFrame
    - idx, the indices of the DataFrame in which we are interested, array 2x1 [idx_start,idx_end]
    - ref : index in the csv file 
    - objet : object column in the csv file 
    - comments, envetual comments to add, must be a string
    - facq, acquisition frequency, must be a float
    """
    
    # Get the date and time
    date = data['CUSTOM.date [local]'][idx[0]]
    date = convert_date(date) # convert date in stuitable format
    
    UTC_t = np.array(['a','a'])
    for i in range(len(idx)):
        t = data['CUSTOM.updateTime [local]'][idx[i]]
        t = convert_time(t) # Convert time in suitable format
    
        # Convert date and time in UTC format
        UTC_t[i] = convert_UTC(date,t)
    
    # Get other variables
    Z = float(data['OSD.height [ft]'][idx[0]])*0.308
    latitude = data['OSD.latitude'][idx[0]]
    longitude = data['OSD.longitude'][idx[0]]
    elevation = float(data['OSD.altitude [ft]'][idx[0]])*0.308
    theta = data['GIMBAL.pitch'][idx[0]]

    # Build a line to be written in a csv 
    line_csv = [ref,objet,UTC_t[0],UTC_t[1],'','',Z,latitude,longitude,elevation,comments,facq,theta]
    
    return line_csv

#########################################
#%%%%% MAIN %%%%%
#########################################

path = 'G:/Rimouski_2024/Data/2024/0211/Drones/bernache/flightrecords/'

filename = glob.glob(path + '*[15-26-43].csv')

print(filename)
data = pd.read_csv(filename[0],header = 1, low_memory = False)

# detect beginning of the movie
mask = np.where(data['CAMERA.isVideo'] == 1)[0]
idx_start = np.min(mask)

# detect end of the movie 

idx_end = idx_start
cam_bool = data['CAMERA.isVideo'][idx_end]

while cam_bool :
    idx_end += 1
    cam_bool = data['CAMERA.isVideo'][idx_end]
    
indices = np.array([idx_start,idx_end])

# Collect data associated to the movie 
ref = 'D001'
objet = 'Stereo 001'
comments = ''
facq = 30
line_csv = get_timeline_row(data,indices,ref,objet,comments,facq)

# Generate headers for the csv file 
header_0 = ['Instrument','','Temps','','Geometrie','','','Position','','','Texte','Variables','']
header_1 = ['Index','Objet','T_0','T_F','X','Y','Z','Latitude','Longitude','Elevation','Commentaire','Facq','Theta']

# Create the csv file and fill it 
csvname = 'Timeline_Drone_stereo001_bernache.csv'
fullname = path + csvname

with open(fullname, 'w', newline="") as file:
    csvwriter = csv.writer(file) # 1. create a csvwriter object
    csvwriter.writerow(header_0) # 2. write the header
    csvwriter.writerow(header_1) # 3. write second header
    csvwriter.writerow(line_csv) # 4. write data





