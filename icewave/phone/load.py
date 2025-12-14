import glob
import os
import numpy as np
import time

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
    #assume that the data for each phone are stored in a separated folder
    folderlist = get_folderlist(folder)
    data = {}
    for folder in folderlist:
        phone = get_number(folder)
        print(f"Load data for phone {phone}")
        data[phone] = load_folder(folder+'/')
    return data

def load_folder(path,header_only=False,date=None):
    if date is not None:
        datafiles = glob.glob(path+f'*{date}*.csv')
    else:
        datafiles = glob.glob(path+'*.csv')
    return load_files(datafiles,header_only=header_only)

def extract_var(filename):
    if "-gps-" in filename:
        var = 'gps'
    elif "-Usb-" in filename:
        var = 'u'
    else:
        name = filename.split('.')[-2].split('-')[0]
        #print(name)
        var = table[name]
    return var

def load_files(datafiles,header_only=False):
    datas={}
    for filename in datafiles:
        print(filename)
        if header_only:
            max_rows=1
        else:
            max_rows=None
        data = load_gobfile(filename,max_rows=max_rows)
        datas = update_datafile(datas,data)
    return datas

def update_datafile(datas,data):
    #data1 : dictionnary containing data. Typical keys : ta, ax, gy, mz, tgps, gpslat, ...
    #data2 :  dictionnnary containing data. same type as data1, key may be different
    #add data2 at the end of data1, for each key present in both dictionnnary
    for key in data:
        if not key in datas:
            datas[key]=data[key]
        else:
            if key in exclude:#don't duplicate coords key
                continue
            print(key)
            datas[key] = np.asarray(list(datas[key])+list(data[key]))#concatenate both dataset
    return datas

def load_gobfile(datafile,max_rows=None):
    data = {}
    data['filename']=[datafile]
    var = extract_var(datafile)

    if var=='gps':
        raw = np.loadtxt(datafile,usecols=(0,1,2,3),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        else:
            i0 = len(raw[:,0])
        if i0>0:
            data['t'+var]=raw[:i0,0]
            data[var+'lat']=raw[:i0,1]
            data[var+'lon']=raw[:i0,2]
            data[var+'elev']=raw[:i0,3] # is it useful ? keep it
        else:
            print('gps data empty')
        return data
    
    elif var=='u':
        raw = np.loadtxt(datafile,usecols=(0,1),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        elif np.sum(raw[:,0]<0)>0:
            print('negative time detected in usb file, cutting')
            i0 = np.where(raw[:,0]<0)[0][0]
        else:
            i0 = len(raw[:,0])
        data['t'+var]=raw[:i0,0]
        data[var+'_raw']=raw[:i0,1]
    else:
        data['coords']=['x','y','z']
        raw = np.loadtxt(datafile,usecols=(0,1,2,3),delimiter=',',skiprows=0,max_rows=max_rows)#
        if np.sum(np.isnan(raw))>0:
            i0 = np.where(np.isnan(raw)[:,1])[0][0]
        else:
            i0 = len(raw[:,0])
        data['t'+var]=raw[:i0,0]
        for i,c in enumerate(data['coords']):
            data[var+c]=raw[:i0,i+1]
    return data

def sync_time(data,tsync=None):
    #provide the path to a Tsync file, if not, use by default the time of creation of the folder ("sync_method"='local')

    if tsync is None:
        #retrieve start time from folder name
        # use first filename as reference (accelerometer ?)
        fileref = data['filename'][0]
        if 'accelerometer' in fileref:
            t0 = data['ta'][0]           
        else:
            print('check data, accelerometer data may be missing')
        s = '20'+fileref.split('_D20')[1]
        strtime = '-'.join(s.split('-')[:3])#extract only the time stamp
        utc_time = time.strptime(strtime, "%Y-%m-%dT%H_%M_%S") #convert to utc time. 
        t1 = time.mktime(utc_time)
        tlag = t1 - t0/1e6
        #print(tlag)
    else:
        #load the time sync folder and do something with it
        tlag = 0
        print('Tsync function to be implemented')
        pass
    key_times,variables = get_times(data)
    for key in key_times:
        data[key]=data[key]/1e6+tlag
    return data

def get_times(data):
    keys = [key for key in data.keys() if key[0]=='t']
    variables = [key[1:] for key in keys]
    return keys,variables

def stat(data,date=False):
    stats = {}
    keys,variables = get_times(data)
    keyref = keys[0]
    print('')
    print('date : '+dispdate(data[keyref][0]))
    for key,var in zip(keys,variables):
        tmin = data[key][0]
        tmax = data[key][-1]
        n = len(data[key])
        print(var,n,disptime(tmin,date=date),disptime(tmax,date=date))

def dispdate(epoch):
    return time.strftime('%Y-%m-%d', time.gmtime(epoch))

def disptime(epoch,date=True):
    if date:
        return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(epoch))
    else:
        return time.strftime('%H:%M:%S', time.gmtime(epoch))

def get_time(data):
    if 'tgps' in data.keys() and len(data['tgps'])>0:
        key = 'tgps'
    elif 'ta' in data.keys():
        print('no gps time, use accelerometer time')
        key = 'ta'
    else:
        print('no time stamp found, check if the data are empty')
        return np.nan,np.nan,np.nan
    t0 = data[key][0]
    t1 = data[key][-1]
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
