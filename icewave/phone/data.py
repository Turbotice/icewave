import pathlib
import urllib.request
import re
from pprint import pprint
import numpy as np
import time

from multiprocessing import Process
import os

import icewave.tools.rw_data as rw

import phonefleet.connect as connect
global port
port = 8080
global folder
folder = '/media/turbots/BlueDisk/Shack25_local/'

def savefolder(date):
    return folder+f'Data/{date}/Phone/'

def list_files(phone,prefix='',date='',ranges=[]):
    ip = connect.get_adress(phone)
    a = urllib.request.urlopen(f"http://{ip}:{port}/list-files").read()
    s = a.decode('utf-8')
    filelist = s[1:-1].split(', ')

    if len(prefix)>0:
        if date=='':
            date = time.today()
        rx = re.compile(f'-{date}T*')
        filelist = list(filter(rx.search, filelist)) 

    pprint(f"Phone {phone}, number of files : {len(filelist)}")

    if len(ranges)>0:
        file_to_save=[]
        for filename in filelist:
            lis = filename.split('-')
            if 'gps' in filename:
                num = int(lis[5])
            else:
                num = int(lis[5])
            #print(len(lis),num,filename)
            if num>=ranges[0] and num<=ranges[-1] and len(lis)<10:
                #print(filename)
                file_to_save.append(filename)
            #files[phone]=filelist
        #print(file_to_save)
        return file_to_save
    else:
        return filelist

def save_file(phone,filename,date):
    data = get_file(phone,filename)
    name = list(data['raw'].keys())[0]
    filename = savefolder(date)+f'{phone}/{name}.csv'
    rw.write_csv(filename,data['raw'][name])
    print(f"Saved, {filename}")

def save_files_onedate(phone,date,imin=0,year='2025',ranges=[]):
    dategob = f'{year}-{date[:2]}-{date[2:]}'
    file_to_save = list_files(phone,prefix='-',date=dategob,ranges=ranges)
    for i,filename in enumerate(file_to_save):
        if i>=0:
            print(i,filename)
            filename_copy = filename
            #print(data['raw'].keys())
            name = '-'.join(filename.split('.')[:-1])
            filename = savefolder(date)+f'{phone}/{name}.csv'
            print(filename)

            if not os.path.dirname(filename):
                os.makedirs(os.path.dirname(filename))
            path = pathlib.Path(filename)
            if not path.exists():
                a = get_file(phone,filename_copy)
                with open(filename,'wb') as f:
                    f.write(a)
            else:
                print(f"{path} already exists, skipping")
#            np.savetxt(filename, a, delimiter=",")
            #rw.write_csv(filename,dic)

def get_files(phone,filelist):
    data = {}
    data['raw'] = {}
    ip = connect.get_adress(phone)

    for i,filename in enumerate(filelist):
        #print(i,filename)
        a = urllib.request.urlopen(f"http://{ip}:{port}/get-file/"+filename).read()
        name = filename.split('.')[-2]
        data['raw'][name] = a
    return data

def get_file(phone,filename):
    ip = connect.get_adress(phone)
    data = {}
    data['raw'] = {}
    try:
        a = urllib.request.urlopen(f"http://{ip}:{port}/get-file/"+filename).read()
    except urllib.request.http.client.IncompleteRead as e:
        print(f"Phone {phone}, file {filename} incomplete")
        a = e.partial
    name = filename.split('.')[-2]
    #data = decode(a)
    #dic = read_data(data,name)
    return a

def load_data(filename):
    #print(num)
    #phone = int(filename.split('/')[-2])
    num = int(filename.split('-')[1])
    data = rw.read_csv(filename,delimiter=',')
    name = os.path.basename(filename)
    dic = read_data(data,name)
    return dic

def read_data(data,name):
    if 'gps' in name:
        dic = read_gps(data)
        return dic
    elif 'accelerometer' in name:
        var='a'
    elif 'gyroscope' in name:
        var='g'
    elif 'magnetic' in name:
        var='m'
    else:
        print('Data type not recognized')
        print(name)
        return None
    dic = read_raw(data,var)
    return dic
    
def read_gps(data):
    var = 'gps'
    dic = {}
    try:
        float(data[0,0])
        j=0
    except:
        j=1

    dic['t'+var] = np.asarray(data)[j:,0].astype(float)/10**6
    dic[var+'lat'] = (np.asarray(data)[j:,1]).astype(float)
    dic[var+'lon'] = (np.asarray(data)[j:,2]).astype(float)
    dic[var+'elev'] = (np.asarray(data)[j:,3]).astype(float)
    return dic

def decode(data):
    d = data.decode('utf-8')
    lines = d.split('\n')
    d = np.asarray([line.split(', ')[:-1] for line in lines]).astype(float)
    return d

def read_raw(d,var):
    try:
        d = np.asarray(d)
    except:
        print('dimension does not match ! Skipping')
        return {}
    dic = {}
    try:
        float(d[0,0])
        j=0
    except:
        j=1
        
    dic['t'+var] = np.asarray(d[j:,0]).astype(float)/10**6
    dic['coords']=['x','y','z']
    for i,c in enumerate(dic['coords']):
        dic[var+c]= np.asarray(d[j:,i+1]).astype(float)
    return dic

if __name__=='__main__':
    print('load data')
    phonelist = range(21,43)
    ranges = []#[23,47]
    for phone in phonelist:
        print(f'load data for phone {phone}')
        date = '0204'
        t1 = Process(target=lambda phone,date:save_files_onedate(phone,date,imin=0,year='2025',ranges=ranges),args=(phone,date))
        t1.start()
