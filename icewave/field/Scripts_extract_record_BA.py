# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:25:04 2024

@author: Banquise
"""
import pandas
import pickle 
import cv2 
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters
from skimage import feature
from scipy.signal import convolve2d
from scipy.signal import savgol_filter, gaussian
from scipy.signal import medfilt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy import stats
import scipy.fft as fft
import os
from PIL import Image

import baptiste.display.display_lib as disp
import baptiste.experiments.import_params as ip
import baptiste.signal_processing.fft_tools as ft 
import baptiste.image_processing.image_processing as imp
import baptiste.math.fits as fits
import baptiste.math.RDD as rdd
import baptiste.files.save as sv
import baptiste.files.dictionaries as dic
import baptiste.tools.tools as tools


#%% read file line by line

# 

path = 'K:\Share_hublot\\Data\\0221\\Drones\\bernache\\pos.txt'




file = open( path, "r")
lines = file.readlines()
lines = np.array(lines)

file.close()


lines[int(len(lines)/ 2)]
#%% Read str file (position drone) : lat and long for all videos at half movie

pathh = 'K:\Share_hublot\\Data\\0220\\Drones\\mesange\\'

dossiers = os.listdir(pathh)[2:-4]
long = []
lat = []
for i in dossiers :
    fichiers = os.listdir(pathh + i)
    for j in fichiers :
        if j[-3:] == 'SRT' :
            file = open( pathh + i + "\\" + j ,"r")
            lines = file.readlines()
            lines = np.array(lines)
    
            file.close()
            if lines[int(len(lines) / 6 / 2) + 1][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 1][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 1][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 2][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 2][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 2][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 3][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 3][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 3][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 4][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 4][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 4][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 5][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 5][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 5][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 6][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 6][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 6][-88:-79]]


for i in range (len(lat)) :
    u = lat[i].find(']')
    if u > 0 :
        lat[i] = lat[i][:-1]

for i in range (len(long)) :
    u = long[i].find(']')
    if u > 0 :
        long[i] = long[i][:-1]

print(lat)
print(long)

#%% read one specfic str
    
pathh = 'K:\Share_hublot\\Data\\0220\\Drones\\mesange\\'

dossiers = os.listdir(pathh)

i = dossiers[30]
fichiers = os.listdir(pathh + i)
for j in fichiers :
    if j[-3:] == 'SRT' :
        file = open( pathh + i + "\\" + j, "r")
        lines = file.readlines()
        lines = np.array(lines)

        file.close()

        print(lines[int(len(lines) / 6 / 2) + 1])
        print(lines[int(len(lines) / 6 / 2) + 2])
        print(lines[int(len(lines) / 6 / 2) + 3])
        print(lines[int(len(lines) / 6 / 2) + 4])
        print(lines[int(len(lines) / 6 / 2) + 5])
        print(lines[int(len(lines) / 6 / 2) + 6])
 
    
#%% lat et long moyenne

path = 'K:\Share_hublot\\Data\\0226\\Drones\\mesange\\pos.txt'

pos = np.loadtxt(path)

mean_lat = np.mean(pos[:,0])
        
mean_long = np.mean(pos[:,1])




#%% path drone

path_drone = 'K:\Share_hublot\\Summary\\'

a = pandas.read_csv(path_drone + '0220_path_drone.txt', sep = ' ', header = None)

path_data = 'K:\Share_hublot\\Data\\0220\\Drones\\bernache'

files = os.listdir(path_data)

a = np.asarray(a)

datas = np.array([])

for i in a[:,1] :
    if i in files :
        datas = np.append(datas, path_data + '\\' + i)
        
np.savetxt(path_drone + 'dronesss.txt', datas)


#%% Timetable

date = '0223'

path_drone = 'K:\Share_hublot\\Summary\\'

datas = pandas.read_csv(path_drone + date + '_path_drone.txt', header= None)
datas = np.asarray(datas)[:,0]

time = dict()
tf = np.zeros(len(datas), dtype = object)
t0 = np.zeros(len(datas), dtype = object)

u = 0
blbl = np.array([], dtype = object)
for file in datas :
    print(file)
    
    a = os.listdir(file)
    blbl = np.array([])
    height = []
    lat = []
    long = []
    type_file = np.array([])
    for i in a :
        if i[-3:] == 'SRT' :  #for a video
            strr = open( file + '\\' + i, "r")
            lines = strr.readlines()
            lines = np.array(lines)

            strr.close()
            height += [lines[-2][-45:-38]]
            
            #Find t0 and tf
            if lines[-2][-45:-38] != '' :
                if file[37] == 'm' :                     #ATTTENTION : 'm' pour 0221, 0226, 0220, 'M' pour 0223
                    t0[u] = int(i[-17:-11]) - 10000
                    tf[u] = int(lines[-3][-13:-11]+ lines[-3][-10:-8]+lines[-3][-7:-5]) - 10000
                    blbl = np.append(blbl, t0[u])
                    blbl = np.append(blbl, tf[u])
                    type_file = np.append(type_file, 'video')
                    type_file = np.append(type_file, 'video')
                    
                if file[37] == 'b':
                    t0[u] = int(i[-17:-11]) + 50000
                    tf[u] = int(lines[-3][-13:-11]+ lines[-3][-10:-8]+lines[-3][-7:-5]) + 50000
                    blbl = np.append(blbl, t0[u])
                    blbl = np.append(blbl, tf[u])
                    type_file = np.append(type_file, 'video')
                    type_file = np.append(type_file, 'video')
            else :
                print('marche pas')
                blbl = np.append(blbl, [0,0])
                type_file = np.append(type_file, 'video')
                type_file = np.append(type_file, 'video')
                
            #Find latitude and longitude at the middle of the video
            if lines[int(len(lines) / 6 / 2) + 1][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 1][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 1][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 2][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 2][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 2][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 3][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 3][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 3][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 4][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 4][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 4][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 5][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 5][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 5][-88:-79]]
            if lines[int(len(lines) / 6 / 2) + 6][-65:-56] != '' :
                lat += [lines[int(len(lines) / 6 / 2) + 6][-66:-56]]
                long += [lines[int(len(lines) / 6 / 2) + 6][-88:-79]]
                
            
        # t0 and tf for a picture
        elif i[-3:] == 'JPG' :
            if file[37] == 'm':           #ATTTENTION : 'm' pour 0221, 0226, 0220, 'M' pour 0223
                blbl = np.append(blbl, int(i[-17:-11]) - 10000)
                type_file = np.append(type_file, 'photo')
            if file[37] == 'b':
                blbl = np.append(blbl, int(i[-17:-11]) + 50000)
                type_file = np.append(type_file, 'photo')
        else :
            blbl = np.append(blbl,[np.nan,np.nan])
            type_file = np.append(type_file, 'autre')
            type_file = np.append(type_file, 'autre')
            
    for j in range (len(lat)) :
        uuu = lat[j].find(']')
        if uuu > 0 :
            lat[j] = lat[j][:-1]

    for j in range (len(long)) :
        uuu = long[j].find(']')
        if uuu > 0 :
            long[j] = long[j][:-1]
            
    
    tf[u] = np.nanmax(blbl)
    t0[u] = np.nanmin(blbl)
        
    time[file] = {}
    time[file]['time'] = blbl
    time[file]['h_drone'] = height
    time[file]['latitude'] = lat
    time[file]['longitude'] = long
    time[file]['type'] = type_file


    u += 1
    
timetable = np.vstack((datas, t0))
timetable = np.vstack((timetable, tf))

if date == '0220' :
    timetable[1,25] = 193623
    timetable[2,25] = 194232
    
#%% search flight record
date = '0223'
path_drone = 'K:\Share_hublot\\Summary\\'

path_m = 'K:\Share_hublot\\Data\\' + date + '\\Drones\\mesange\\flightrecords\\'

path_b = 'K:\Share_hublot\\Data\\' + date + '\\Drones\\bernache\\flightrecords\\'

list_flight_m = os.listdir(path_m)
list_flight_b = os.listdir(path_b)

list_flight = []
time_flight = []

for j in list_flight_b :
    if j [-3:] == 'csv' :
        list_flight += [path_b + j]
        time_flight += [int(j[-7:-5]) + int(j[-10:-8]) * 100 + (int(j[-13:-11]) + 5) * 10000]
        
for i in list_flight_m :
    if i [-3:] == 'csv' :
        list_flight += [path_m + i]
        time_flight += [int(i[-7:-5]) + int(i[-10:-8]) * 100 + (int(i[-13:-11]) - 1) * 10000]

time_flight = np.array(time_flight)
list_flight = np.array(list_flight)

#%% Open flight record

datas = pandas.read_csv(path_drone + date + '_path_drone.txt', header= None)
datas = np.asarray(datas)[:,0]

def find_flight_record(list_flight, time_flight, path, t) :
    drone = path[37]
    u = []
    tt = []
    j = 0
    for i in range (len(list_flight)) :
        if path[37] == list_flight[i][33]:
            if time_flight[i] > t :
                u += [i]
                tt += [time_flight[i]]    
    j = u[np.argmin(tt)]
    return list_flight[j-1]


def tall (path, tt) :
    flight = pandas.read_csv(path,low_memory=False)    
    heure = []
    latitude = []
    longitude = []
    angle = []
    
    
    
    for i in range (flight.shape[0]-3) :
        uu = flight[i+1:i+2]['sep='][0:1]
        j = np.array(uu.keys()[0])
        h = (int(j[1][0]) + 12) * 10000 + int(j[1][2:4])*100 + int(j[1][5:7])
        if h == tt :
            heure += [h]
            latitude += [j[4]]
            longitude += [j[5]]
            angle += [j[56]]
        
    angle = np.asarray(angle, dtype = float)
    latitude = np.asarray(latitude, dtype = float)
    longitude = np.asarray(longitude, dtype = float)           
            
    
    
    # angle = np.asarray(angle, dtype = float)
    return angle, latitude, longitude

def tall_v2 (path, tt) :
    flight = pandas.read_csv(path,low_memory=False)    
    
    uu = flight[1:2]['sep='][0:1]
    j = np.array(uu.keys()[0])
    
    tt = int(str(tt)[:2]) * 36000 + int(str(tt)[2:4]) * 600 + int(str(tt)[4:6]) * 10
    t0 = (int(j[1][0]) + 12) * 36000 + int(j[1][2:4])*600 + int(j[1][5:7]) *10
    print('t0 = ', t0)
    
    indice_t = (round(tt - t0, -2) * 60 + tt - t0 - round(tt - t0, -2) ) * 10
    
    uu = flight[indice_t:indice_t+1]['sep='][0:1]
    j = np.array(uu.keys()[0])
    latitude = j[4]
    longitude = j[5]
    angle = j[56]
    
    # for i in range (flight.shape[0]-3) :
    #     uu = flight[i+1:i+2]['sep='][0:1]
    #     j = np.array(uu.keys()[0])
    #     h = (int(j[1][0]) + 12) * 10000 + int(j[1][2:4])*100 + int(j[1][5:7])
    #     if h == tt :
    #         heure += [h]
    #         latitude += [j[4]]
    #         longitude += [j[5]]
    #         angle += [j[56]]
        
    # angle = np.asarray(angle, dtype = float)
    # latitude = np.asarray(latitude, dtype = float)
    # longitude = np.asarray(longitude, dtype = float)           
            
    
    
    # angle = np.asarray(angle, dtype = float)
    return angle, latitude, longitude


flight_record = []
for j in range (1) : # (len(datas)) :
    angle = []
    latitude = []
    longitude = []
    for i in range(len( time[datas[j]]['time'])) :
        tt = time[datas[j]]['time'][i]
        if time[datas[j]]['type'][i] == 'photo' :
            print(tt)
            flight_record += [find_flight_record(list_flight, time_flight, datas[i], tt)]
            b,c,d = tall_v2(flight_record[-1], tt)
            angle += [b]
            latitude += [c]
            longitude += [d]
    print(angle)
    time[datas[j]]['angle'] = angle
    time[datas[j]]['latitude'] = latitude
    time[datas[j]]['longitude'] = longitude
    time[datas[j]]['flight_record'] = flight_record
    




    
    
#%% plot timetable

t = range(int(np.nanmin(t0)), int(np.nanmax(tf)), 1)

disp.figurejolie(width = 15)
if date == '0220' :
    timetable[1,25] = 193623
    timetable[2,25] = 194232

for j in range (timetable.shape[1]) :
    t00 = timetable[1,j]
    tff = timetable[2,j]
    tt = np.array([])
    x_t = np.array([])
    if timetable[0,j][37] =='m' :
        for i in t :
            if i <= tff and i >= t00 :
                tt = np.append(tt,i)
                x_t = np.append(x_t,1)
        disp.joliplot('t (UTC)', 'Instrument', tt,x_t, exp = False, color = 5, linewidth = 10)
        plt.annotate(timetable[0,j][46:48], [np.min(tt), 1])
    if timetable[0,j][37] =='b' :                                       #ATTENTION des fois 'b' des fois 'B', dépend de la date
        for i in t :
            if i <= tff and i >= t00 :
                tt = np.append(tt,i)
                x_t = np.append(x_t,0)
        disp.joliplot('t (UTC)', 'Instrument', tt,x_t, exp = False, color = 2, linewidth = 10)
        plt.annotate(timetable[0,j][47:49], [np.min(tt), 0])
    
plt.xlim(min(t), max(t))  #le temps est en 195030 (hhmmss) ce qui n'a aucun sens, il faudrait plot en UTC bien mais je sais pas faire
plt.yticks([0,1,2,3, 4], ['Bernache', 'Mésange', 'G', 'WB', 'T'])


#%% save timetable
path_save = 'K:\\Share_hublot\\Summary\\Timeline\\'

to_save = np.rot90(timetable)

df = pandas.DataFrame(to_save)

df.to_csv(path_save + date + "_timetable.txt", sep='\t', index=False)

dic.save_dico(time, path = path_save + date + '_times.pkl')
