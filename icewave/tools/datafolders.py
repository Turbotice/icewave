import numpy as np
import pylab as plt
import glob
import os
import shutil
import time
from pprint import pprint

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

print(ostype)
print(osname)

def find_path(disk='Backup25',year='2025',smb=False,date='0211'):
    if year=='2025':
        if 'macOS' in ostype:
            pass#disk='F-1'
        elif 'Linux' in ostype:
            disk='F'#Shack25'
    #print('OS type : '+str(ostype))
    #print('Computer name : '+str(osname))
    #if year=='2025':
    #    dateday = (int(date[:2])-1)*31+int(date[2:])
    #    if dateday>(31+18):
    #        disk = 'F'
    #    else:
    #        disk = 'Shack25'
    if disk=='BicWin2024':
        base = disk+'/Share/Data/'
    elif disk=='Hublot24' or disk=='Elements':
        base = disk+'/Share_hublot/Data/'
    elif disk=='Shack25' or 'F':
        base = disk+'/Data/'
    else:
        print('please specify the folder path')
    if 'Linux' in ostype:
        if 'oural' or 'saguenay' in osname:
            if disk=='Shack25':
                base = '/media/turbots/BlueDisk/Shack25_local/Data/'
                if smb:
                    base = '/home/turbots/Documents/Bic25/Data/'
            else:
                base = '/media/turbots/'+base#home/turbots/Documents/'+base
        else:
            print('computer unknown, define the path folder')
            base = ''
        #if 'saguenay' in osname:

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
        print('ostype : Windows')
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
        if disk=='Hublot24' or disk == 'Elements':
            base = 'K:/Share_hublot/Data/'
        if disk=='Shack25':
            base = 'F:/Data/'
        if disk == 'Backup25':
            base = 'U:/Data/'
        if disk=='F':
            base = 'F:/Data/'

    if 'macOS' in ostype:
        base = '/Volumes/'+base
        #if 'laita' in osname:
        #    if disk =='BicWin2024':
        #        base = '/Volumes/Share_hublot/Data/'
        #    elif disk=='Hublot24':
        #        base = disk+'/Share_hublot/Data/'
        #    else:
        #        base = '~/Documents/'+base #praise UNIX system    
        #else:
        #    base = '/Users/stephane/Documents/git/icewave/icewave/field/Bicwin2024/Data/2024/'#/Volumes/Share_hublot/Data/'

    #browse.create_folder(base)
    #print(base)
    return base

def date_folder(date=''):
    #date should be in format 2024_0205
    base = find_path()
    #print(base)
    if date=='':
        date = get_current_date()
    year,day = date.split('_')

    #browse.create_folder(base+year)
    browse.create_folder(base+day)
    return base+day+'/'

def find_instruments(date=''):
    folder = date_folder(date=date)

    folderlist = glob.glob(folder+'*')
    instruments = [os.path.basename(f) for f in folderlist]
    #pprint(instruments)
    return instruments

def get_current_date():
    t = time.localtime()
    date = str(t[0])+'_'+ndigit(t[1])+ndigit(t[2])

    #print(t)
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
    print(ostype)
    if 'Linux' in ostype:
        base = '/media/turbots/'
    if 'macOS' in ostype:
        base = '/Volumes/'
    return base

if __name__=='__main__':
    print(date_folder())
