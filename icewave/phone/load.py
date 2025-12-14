import glob
import os
import numpy as np
import csv

import zipfile

global table
global exclude
table = {'Accelerometer':'a','Gyroscope':'g','Location':'l','Magnetometer':'m','Usb':'u'}
table.update({'accelerometer':'a','gyroscope':'g','magnetic_field':'m','usb':'u'})
exclude = ['coords'] # field to not be concatenated
coords = ['x','y','z']

def get_number(folder):
    #find phone number from the name of the folder
    s = os.path.basename(folder)
    print(s)
    if '_' in s:
        return s 
    else:
        return int(s)
    
def get_pĥonelist(folder):
    folderlist = get_folderlist(folder)
    phonelist = [get_number(f) for f in folders]
    return list(set(phonelist))

def get_folderlist(folder):
    folderlist = glob.glob(folder+'/*/')
    folderlist = [f for f in folderlist if len(glob.glob(f+'/*.csv'))>0]
    return folderlist

def loads(folder,header=False):
    folderlist = get_folderlist(folder)
    data = {}
    for folder in folderlist:
        phone = get_number(folder)
        print(f"Load data for phone {phone}")
        if header:
            data[phone] = load_header(folder,max_rows=1)
        else:
            data[phone] = load_gobfile(folder)
    return data

def concatenate(data1,data2):
    #data1 : dictionnary containing data. Typical keys : ta, ax, gy, mz, tgps, gpslat, ...
    #data2 :  dictionnnary containing data. same type as data1, key may be different
    #add data2 at the end of data1, for each key present in both dictionnnary
    for key in data1.keys():
        if key in exclude:
            continue
        if key in data2.keys():
            data1[key] = np.asarray(list(data1[key])+list(data2[key2]))#concatenate both dataset
    for key in data2.keys():
        if key in exclude:
            continue
        if not key in data1.keys():
            data1[key]=data2[key]
    return data1

def load_gobfile(datafile,max_rows=None):
    data = {}
    if "-gps-" in datafile:
        var = 'gps'
        raw = np.loadtxt(datafile,usecols=(0,1,2,3),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        else:
            i0 = len(raw[:,0])
        data['t'+var]=raw[:i0,0]
        data[var+'lat']=raw[:i0,1]
        data[var+'lon']=raw[:i0,2]
        data[var+'elev']=raw[:i0,3] # is it useful ? keep it 
        return data
    elif "-Usb-" in datafile:
        name = datafile.split('_')[-1].split('-')[1]
        print(name)
        var = table[name]
        raw = np.loadtxt(datafile,usecols=(0,1),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        else:
            i0 = len(raw[:,0])
        data['t'+var]=raw[:i0,0]
        data[var+'_raw']=raw[:i0,1]
    else:
        data['coords']=['x','y','z']
        name = datafile.split('.')[-2].split('-')[0]
        #print(name)
        var = table[name]
        raw = np.loadtxt(datafile,usecols=(0,1,2,3),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        else:
            i0 = len(raw[:,0])
        data['t'+var]=raw[:i0,0]
        for i,c in enumerate(data['coords']):
            data[var+c]=raw[:i0,i+1]
    return data
    
def get_time(data):
    if 'tgps' in data.keys():
        key = 'tgps'
    elif 'ta' in data.keys():
        print('no gps time, use accelerometer time')
        key = 'ta'
    else:
        print('no time stamp found, check if the data are empty')
        return np.nan,np.nan,np.nan
    t0 = data['tgps'][0]
    t1 = data['tgps'][-1]
    Dt = t1 - t0
    return t0,t1,Dt

def get_mean_position(data):
    if not 'loc' in data.keys():
        data = sort(data)
    if 'loc' in data.keys():    
        Lat = np.round(np.mean(data['loc']['lat']),decimals=7)
        Lon = np.round(np.mean(data['loc']['lon']),decimals=7)
    else:
        print('No location available')
        return np.nan,np.nan
    return Lat,Lon
