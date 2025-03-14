
import glob
import h5py
import numpy as np
from pprint import pprint
import os
import datetime


import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    files = glob.glob(base+date+'/boueeVague/*/mat/*.mat')#/*/*.srt')
    nbase = len(base)
    
    records = {}
    records['buoys'] = {}
    for filename in files:
        record,name = read_matfile(filename)
        path = filename[nbase:]
        if record==None:
            continue
        if not name in records['buoys'].keys():
            records['buoys'][name] = {}
        key = os.path.basename(filename).split('.')[0]
        records['buoys'][name][key]=record
        records['buoys'][name][key]['path']=path
        
    return records

def load_data(record):
    base = df.find_path()
    filename = base + record['path']
    return read_buoy_data(filename)

def read_buoy_data(filename):
    print(filename)
    buoy = h5py.File(filename)
    return buoy

def read_matfile(filename):
    name = filename.split('/')[-3]
    buoy = read_buoy_data(filename)
    
    record={}
        #location
    try:
        record['latitude']= [np.mean(buoy['IMU']['GPS1_POS']['LAT'][0])]#.keys()
        record['longitude']= [np.mean(buoy['IMU']['GPS1_POS']['LONG'][0])]
    except:
        print(buoy['IMU'].keys())
        return None,None
        
    #time 
    hours = buoy['IMU']['UTC_TIME']['HOUR'][:][0]#.keys()
    mins = buoy['IMU']['UTC_TIME']['MIN'][:][0]#.keys()
    secs = buoy['IMU']['UTC_TIME']['SEC'][:][0]#.keys()
    times = [str(int(hour))+':'+str(int(m))+':'+str(sec).replace('.','')[:2] for (hour,m,sec) in zip(hours,mins,secs)]
    record['time']=times
    return record,name 

def get_time(buoy):
    a = buoy['IMU']['UTC_TIME']
    time_stamp = f"{int(a['YEAR'][0,0])}-{int(a['MONTH'][0,0])}-{int(a['DAY'][0,0])} {int(a['HOUR'][0,0])}:{int(a['MIN'][0,0])}:{int(a['SEC'][0,0])}.{int(a['NANOSEC'][0,0])}"
    #time.struct_time(tm_year=a['YEAR'],tm_mon=a['MONTH'],tm_day=a['DAY'],tm_hour=a['HOUR'],tm_min=a['MIN'],tm_sec=a['SEC'])
    t00, t01 = time_stamp.split(".")
    date = datetime.datetime.strptime(t00.split(" UTC")[0], "%Y-%m-%d %H:%M:%S")
    tbuoys = date.timestamp() + int(t01)/1000
    ts = datetime.datetime.fromtimestamp(tbuoys)
    return (tbuoys,ts)
