import glob
import os
import numpy as np
import csv

import zipfile

global table
table = {'Accelerometer':'a','Gyroscope':'g','Location':'l','Magnetometer':'m'}
coords = ['x','y','z']

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

    print('toto')
    data = {}
    for meta in metalist:
        key = '_'.join(os.path.basename(meta).split('.')[:-1])
        if key=='device':
            data[key]  = read_meta_device(meta)
        elif key=='time':
            data[key]  = read_meta_time(meta)
        else:
            print(key)
            print('meta csv file unrecognized')

    for datafile in datalist:
        key = '_'.join(os.path.basename(datafile).split('.')[:-1])
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

    newdata['coords']=coords
    return newdata

def split(data):
    #split into several recordings multiple recordings that are stacked together
    #test if multiple START times':
    if 'system time_START_1' in data['time']:
        print('multiple recordings in the same .csv files')
        #find the number of recordings
        found = True
        c = 1
        while found:
            found = 'system time_START_'+str(c) in data['time']
            c+=1
            print(c)
        c=c-1
        print(f'Number of recordings : {c}')

        datas = []
        for i in range(c):
            d = {}
            d['device'] = data['device']
            d['time']={}

            for key in data['time'].keys():
                if i==0:
                    if key[-5:]=='PAUSE' or key[-5:]=='START':
                        d['time'][key]=data['time'][key]
                else:
                    if f'_{i}' in key:
                        k = '_'.join(key.split('_')[:-1])
                        d['time'][k]=data['time'][key]

            t0 = float(d['time']['experiment time_START'])
            tf = float(d['time']['experiment time_PAUSE'])
            print(t0,tf)

            for key in data.keys():
                if not key in ['device','time']:
                    d[key]=data[key]
            datas.append(d)
    else:
        print('no multiple recordings found')
        return [data]
    return datas

def get_time_interval(data):
    t0 = float(data['time']['experiment time_START'])
    tf = float(data['time']['experiment time_PAUSE'])
    return t0,tf
    
def get_time(data):
    if 'time' in data.keys():
        #print(data['time'].keys())

        found = False
        for key in ['system_START','system time_START']:
            #print(key,data['time'].keys())
            if key in data['time'].keys():
                t0 = data['time'][key]
                found = True
        if not found:
            print('no t0 stamp for this phone')
            t0 = np.nan

        for key in ['system_PAUSE','system time_PAUSE']:
            if key in data['time'].keys():
                t1 = data['time'][key]
                found = True
        if not found:
            print('no t1 stamp for this phone')
            t1 = np.nan

        found=False
        for key1,key2 in zip(['experi_START','experiment time_START'],['experi_PAUSE','experiment time_PAUSE']):
            if key1 in data['time'].keys():
                Dt = np.round(float(data['time'][key2])-float(data['time'][key1]),decimals=2)
                found = True
        if not found:
            print('no duration stamp for this phone')
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

def read_meta_device(csvfile):
    d={}
    with open(csvfile) as f:
        csv_reader = csv.reader(f, delimiter=',')
        line=0
        for j,row in enumerate(csv_reader):
            if line==0:
                headers = row
            else:
                #print(headers,row)
                for i,header in enumerate(headers):
                    key = header[:6]+'_'+row[0]
                    if not key in d.keys():
                        d[key]=row[i]
            line+=1
    return d

def read_meta_time(csvfile):
    d={}
    with open(csvfile) as f:
        csv_reader = csv.reader(f, delimiter=',')
        line=0
        for j,row in enumerate(csv_reader):
            if line==0:
                headers = row[1:]
            else:
                for i,header in enumerate(headers):
                    #print(header) 
                    key = header+'_'+row[0]#beware of back compatibility !!! header[:6] replaced by header
                    if not key in d.keys():
                        d[key]=row[i+1]
                    else:
                        num = int((j-1)/2)
                        d[key+'_'+str(num)] = row[i+1]
                        print('key already exist ! Multiple recordings, num = '+str(num+1))
                
            line+=1
    return d



