
import glob
import h5py
import numpy as np
from pprint import pprint

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    files = glob.glob(base+date+'/boueeVague/*/mat/*.mat')#/*/*.srt')

    records = {}
    records['buoys'] = {}
    for filename in files:
        record,name = read_matfile(filename)
        if not name in records['buoys'].keys():
            records['buoys'][name] = [record]
        else:
            records['buoys'][name].append(record)
    return records

def read_matfile(filename):
    print(filename)
    name = filename.split('/')[-3]
    buoy = h5py.File(filename)
    record={}
    #time 
    hours = buoy['IMU']['UTC_TIME']['HOUR'][:][0]#.keys()
    mins = buoy['IMU']['UTC_TIME']['MIN'][:][0]#.keys()
    secs = buoy['IMU']['UTC_TIME']['SEC'][:][0]#.keys()
    times = [str(int(hour))+':'+str(int(m))+':'+str(sec).replace('.','')[:2] for (hour,m,sec) in zip(hours,mins,secs)]
    record['time']=times

    #location
    record['latitude']=np.mean(buoy['IMU']['GPS1_POS']['LAT'][0])#.keys()
    record['longitude']=np.mean(buoy['IMU']['GPS1_POS']['LONG'][0])

    return record,name 