import glob
import os
import numpy as np
import csv

import zipfile

global table
table = {'Accelerometer':'a','Gyroscope':'g','Location':'l','Magnetometer':'m'}
coords = {'x','y','z'}

def extract_all(folder):
    filelist = glob.glob(folder+'*.zip')
    print(f"Extract : {filelist}")

    for filename in filelist:
        with zipfile.ZipFile(filename,"r") as zip_ref:
            foldersave = filename.split('.')[0]
            zip_ref.extractall(foldersave)

def get_pĥonelist(folder):
    folders = glob.glob(folder+"*")
    phonelist = [get_number(f) for f in folders]
    return list(set(phonelist))

def get_folderlist(folder):
    folderlist = glob.glob(folder+'/*/')
    folderlist = [folder[:-1] for folder in folderlist]
    print(folderlist)
    return folderlist

def get_number(folder):
    print(os.path.basename(folder))
    if os.path.basename(folder).split('_')[0]=='test':
        return int(os.path.basename(folder).split('_')[-1])-100
    else:
        return int(os.path.basename(folder).split('_')[1])

def loads(folderlist,header=False):
    data = {}

    for folder in folderlist:
        phone = get_number(folder)
        print(f"Load data for phone {phone}")
        if header:
            data[phone] = load_header(folder)
        else:
            data[phone] = load(folder)
    return data

def load_header(folder):
    metalist = glob.glob(folder+'/meta/*.csv')
    datalist = glob.glob(folder+'/Location.csv')

    data = {}
    for meta in metalist:
        key = os.path.basename(meta).split('.')[0]
        data[key]  = read_meta(meta)

    for datafile in datalist:
        key = os.path.basename(datafile).split('.')[0]
        
        if key in table:
            key = table[key]
       #try:
        print(key)
        print(datafile)
        data[key]={}
        data[key]['d'] = np.loadtxt(datafile, delimiter=',',usecols=(1,2,3),skiprows=1)
            #get times
        data[key]['t'] = np.loadtxt(datafile, delimiter=',',dtype=str,usecols=(0),skiprows=1)
        #except:
            #print('data cannotbe converted to float')
    return data   

def load(folder):
    metalist = glob.glob(folder+'/meta/*.csv')
    datalist = glob.glob(folder+'/*.csv')

    data = {}
    for meta in metalist:
        key = os.path.basename(meta).split('.')[0]
        data[key]  = read_meta(meta)

    for datafile in datalist:
        key = os.path.basename(datafile).split('.')[0]
        if key in table:
            key = table[key]
       #try:
        print(key)
        print(datafile)
        data[key]={}
        data[key]['d'] = np.loadtxt(datafile, delimiter=',',usecols=(1,2,3),skiprows=1)
            #get times
        data[key]['t'] = np.loadtxt(datafile, delimiter=',',dtype=str,usecols=(0),skiprows=1)
        #except:
            #print('data cannotbe converted to float')
    return data

def sort(data):
    newdata = {}
    keys = {'time','device'}
    for key in keys:
        if key in data.keys():
            newdata[key]=data[key]
        else:
            print(f'no {key} key in the data')
    
    keys = {'a','m','g'}
    
    for key in keys:
        if key in data.keys():
            for i,k in enumerate(coords):
                newdata[key+k]=data[key]['d'][:,i]
            newdata['t'+key] = data[key]['t'].astype(float)

    key = 'l'
    if key in data.keys():
        newdata['loc'] = {}
        newdata['loc']['lat']=data[key]['d'][:,0]
        newdata['loc']['lon']=data[key]['d'][:,1]
        newdata['loc']['elev']=data[key]['d'][:,2]
        newdata['loc']['t']= data[key]['t'].astype(float)
    
    data['coords']=coords
    return newdata

def get_time(data):
    if 'time' in data.keys():
        print(data['time'].keys())

        if 'system_START' in data['time'].keys():
            t0 = data['time']['system_START']
        else:
            print('no t0 stamp for this phone')
            t0 = np.nan
        if 'system_PAUSE' in data['time'].keys():
            t1 = data['time']['system_PAUSE']
        else:
            print('no t1 stamp for this phone')
            t1 = np.nan
        if 'experi_PAUSE' in data['time'].keys():
            Dt = np.round(float(data['time']['experi_PAUSE']),decimals=2)
        else:
            Dt = np.nan
        return t0,t1,Dt
    else:
        print('no time stamp for this phone')
        return np.nan,np.nan,np.nan

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

def read_meta(csvfile):
    d={}
    with open(csvfile) as f:
        csv_reader = csv.reader(f, delimiter=',')
        line=0
        for row in csv_reader:
            if line==0:
                headers = row
            else:
                #print(headers,row)
                for i,header in enumerate(headers):
                    d[header[:6]+'_'+row[0]]=row[i]
            line+=1

    return d

