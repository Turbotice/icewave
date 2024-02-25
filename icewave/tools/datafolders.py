import numpy as np
import pylab as plt
import glob
import os
import shutil
import time

#import garmin
#import fitdecode


import icewave.tools.browse as browse

base = '/Volumes/'

base = '/media/'

folderactivity = '/Volumes/GARMIN/Garmin/Activities/' # on macOS X
folderwpts = '/Volumes/GARMIN/Garmin/GPX/' # on macOS X

#destfolder = '/Users/stephane/Documents/git/icewave/icewave/field/Data/'
#serveurfolder = '/Volumes/labshared2/Banquise/Rimouski 2023/Data/GPS/'


import glob
import os
#Global variables0

import platform
import socket

global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()

def date_folder(date=''):
    #date should be in format 2024_0205
    base = find_path()
    print(base)
    if date=='':
        date = get_current_date()
    year,day = date.split('_')

    #browse.create_folder(base+year)
    browse.create_folder(base+day)
    return base+day+'/'

def get_current_date():
    t = time.localtime()
    date = str(t[0])+'_'+ndigit(t[1])+ndigit(t[2])

    print(t)
    return date    

def date_GPS(date):
    year,day = date.split('_')
    print(year,day)
    month = day[:2]
    day = day[2:]
    return year+'-'+month+'-'+day

def date_phox(date):
    year,day = date.split('_')
    print(year,day)
    month = day[:2]
    day = day[2:]
    return year+'-'+month+'-'+day

def ndigit(num,n=2):
    s = str(num) #no . 
    if len(s)<n:
        for i in range(n-len(s)):
            s='0'+s
    return s

def instrument_folder(key='G',date=''):
    instruments = {'G':'Geophones','T':'Telephones','GPS':'GPS','B':'Buoys','D':'Drone'}
    base = date_folder(date=date)
    if key in instruments.keys():
        browse.create_folder(base+instruments[key])
        return base + instruments[key]+'/'
    else:
        print('Instrument key does not exist yet, please add it in icewave/tools/datafolders')
        return None

def path_GPSdata():
    if 'Linux' in ostype:
        base = '/media/turbots/'
    if 'Darwin' in ostype:
        base = '/Volumes/'
    return base


def find_path(disk='BicWin2024'):
    print('OS type : '+str(ostype))
    print('Computer name : '+str(osname))

    if disk=='BicWin2024':
        base = 'Bicwin2024/Share/Data/'
    else:
        print('please specify the folders path')
    if 'Linux' in ostype:
        #if 'oural' in osname:
        base = '/media/turbots/'+base#home/turbots/Documents/'+base
        #elif 'adour' in osname:
        #    serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base
        #elif 'spi201711-Latitude-5480' in osname:
        #    serveurfolder = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/'+disk+'/'+base
        #elif 'thiou' in osname:
        #    serveurfolder = '/volume3/'+disk+'/'+base #praise UNIX system    		
        #else:
        #    base = '/home/turbots/Documents/'+base

        #    serveurfolder = '/media/turbots/DATA/thiou/'+disk+'/'+base #praise UNIX system
    if 'Windows' in ostype:    #fucking windows OS : beware of the mounting disk you used
        if disk=='labshared2':
            base = 'W:/'+base 
        if disk=='labshared1':
            base = 'Y:/'+base 
        if disk=='storageshared':
            base = 'X:/'+base 
        if disk=='homes':
            base = 'Z:/'+base
        if disk=='BicWin2024':
            base = '//192.168.1.70/Share/Data/'

    if 'Darwin' in ostype:
        if 'laita' in osname:
            if disk =='BicWin2024':
                base = '/Volumes/Share/Data/'
            else:
                base = '~/Documents/'+base #praise UNIX system    
        else:
            base = '/Volumes/Share/Data/'

    browse.create_folder(base)
    return base

if __name__=='__main__':
    print(date_folder())
