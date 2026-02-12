import numpy as np
import glob
import os
from pprint import pprint
import subprocess

import sqlite3

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw


global base #adjust with OS
base = '/Users/stephane/Downloads/SuntooWatch/'

global database
database = base+'stt-dump/amer_app.db'

global savefolder
savefolder = base+'data/'

if not os.path.isdir(savefolder):
    os.makedirs(savefolder)

def download_data(date=''):
    retrieve()
    rows = read_database()
    data = convert_data(rows)
    save(data,date=date)

def retrieve():
    # get the database in the Suunto app folder using adb
    #only work if only one adb device is connected

    android_folder = '/storage/emulated/0/Android/data/com.stt.android.suunto/files/'
    filename = android_folder+'stt-dump/amer_app.db'
    outputfile = '~/Downloads/SuntooWatch/stt-dump/amer_app.db'
    out = subprocess.run(['adb','pull',filename,outputfile],capture_output=True)
    
def read_database():
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute("SELECT * FROM pois")

    rows = cur.fetchall()
    savefolder = base+'data/'
    return rows

def save(data,date=''):
    if date == '':
        date = df.get_current_date()
        print(date)
    filename = savefolder+f'{date}_POIs_Suunto.h5'
    rw.write_h5(filename,data)

def convert_data(rows):
    keys = ['t0','t1','LON','LAT','h','tag']
    keytag_num = 5
    data = {}
    for row in rows:
        d={}
        for i,key in enumerate(keys):
            d[key]=row[i]
            if i==keytag_num:
                k = row[i]
        print(k)
        data[k]=d
    return data

def ls_date(date=''):
    if date=='':
        date = df.get_current_date()
#    dayGPS = date_GPS(date)
#    fitlist = glob.glob(folderactivity+dayGPS+'*.fit')
#    pprint(fitlist)  
#    wayptslist = glob.glob(folderwpts+'Way*'+dayGPS+'*')
#    pprint(wayptslist)
    return None#fitlist,wayptslist
    
def date_watch(date):
    year,day = date.split('_')
    print(year,day)
    month = day[:2]
    day = day[2:]
    return year+'-'+month+'-'+day


def convert_traj(d,tmin=None,tmax=None):
    """
    Convert trajectory from
    INPUT 
    ----- 
    d : dict with key 'r', with fields position_long and position_lat
    tmin : float
    tmax : float
    OUTPUT
    -----
    X,Y : np 1d arrays
    """
    X,Y = [],[]
    for i,a in enumerate(d['r']):
        x = a.get_value('position_long')
        y = a.get_value('position_lat')
        if x is not None:
            if tmin is not None:
                #print(a.get_value('timestamp').time())
                b1 = a.get_value('timestamp').time()>tmin
                b2 = a.get_value('timestamp').time()<tmax
                if not (b1 and b2):
                    continue
            X.append(x/2**32*360)
            Y.append(y/2**32*360)#2^32/360
    X = np.asarray(X)
    Y = np.asarray(Y)
    return X,Y
    #plt.plot(X,Y,'b.')


def get_time(waypoint):
    import datetime
    t1 = waypoint.time
    t = datetime.time(t1.hour,t1.minute,t1.second)
    return t

def get_traj(filename,tmin=None,tmax=None):
    d = decode(filename)
    X,Y = convert_traj(d,tmin=tmin,tmax=tmax)
    return X,Y,d

if __name__ == '__main__':
    download_data()
