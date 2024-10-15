
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

import numpy as np
import os

global base
base = df.find_path(disk='Hublot24')

import argparse
def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate multi instruments data")
    parser.add_argument('-date', dest='date', type=str,default='0226',help='select date to process data')
    #parser.add_argument('-step', dest='step', type=int,default=3,help='select Step to be performed')
#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args


def get_records(date):
    srtfiles = get_srtfiles(date)
    records = {}
    records['drones']={}
    for key in srtfiles.keys():
        records['drones'][key]={}
        for i,srtfile in enumerate(srtfiles[key]):
            print(i,srtfile)
            name = srtfile.split('/')[-2]
            record = get_flighrecord(srtfile,drone=key)
            records['drones'][key][name]=record
    return records
    
def get_srtfiles(date):
    srtfiles = {}
    print(base)
    drones = ['mesange','Bernache','Fulmar']
    for key in drones:
        srt = glob.glob(base+date+'/Drones/'+key+'/*/*.SRT')#/*/*.srt')
        pprint(srt)
        if len(srt)>0:
            srtfiles[key] = srt
        else:
            print(f"No data for {key} on {date}")
    return srtfiles

def convert_flightrecords(date):
    csvfiles = get_csvfiles(date)
    records={}
    records['drones']={}
    for drone in csvfiles.keys():
        records['drones'][drone]={}
        for i,csvfile in enumerate(csvfiles[drone]):
            record=parse_csv_flightrecord(csvfile,drone=drone)
            if i==0:
                records['drones'][drone]=record
            else:
                for key in records['drones'][drone].keys():
                     records['drones'][drone][key]=records['drones'][drone][key]+record[key]
        savedict = records['drones'][drone]
        filename = os.path.dirname(csvfile)+'/Flightrecord_dict.pkl'
        rw_data.write_pkl(filename,savedict)
#    return records
            #merge the record files ?
            
def get_csvfiles(date):
    csvfiles = {}
    print(base)
    drones = ['mesange','Bernache','Fulmar']
    for key in drones:
        filelist = glob.glob(base+date+'/Drones/'+key+'/flightrecords/*.csv')#/*/*.srt')
        pprint(filelist)
        if len(filelist)>0:
            csvfiles[key] = filelist
        else:
            print(f"No flightrecord for {key} on {date}")
    return csvfiles

def parse_csv_flightrecord(csvfile,drone='mesange'):
    table = rw_data.read_csv(csvfile)
    data = rw_data.csv2dict(table,headerindex=1)
    record = {}
    keys_date = ['CUSTOM.date [local]', 'CUSTOM.updateTime [local]']
    keys_bool = ['CAMERA.isPhoto', 'CAMERA.isVideo']
    #, 'OSD.flyTime', 'OSD.flyTime [s]', 'OSD.latitude', 'OSD.longitude', 'OSD.height [ft]']
    keys1 = ['OSD.altitude [ft]', 'OSD.mileage [ft]', 'OSD.hSpeed [MPH]', 'OSD.xSpeed [MPH]', 'OSD.ySpeed [MPH]', 'OSD.zSpeed [MPH]']
    keys2 = ['OSD.pitch', 'OSD.roll', 'OSD.yaw', 'OSD.yaw [360]', 'OSD.gpsNum']
    keys3 = ['GIMBAL.pitch', 'GIMBAL.roll', 'GIMBAL.yaw', 'GIMBAL.yaw [360]']            
    keys_float = keys1+keys2+keys3

    record = {}
    for key in keys_float:
#        print(key,data[key][1000])
        record[key]= [float(d) for d in data[key]]
    for key in keys_date:
        record[key]= data[key]
    for key in keys_bool:
        record[key]= [bool(d) for d in data[key]]
    return record

def get_flighrecord(srtfile,step=100,drone='mesange'):
    #convert all times to UTC
    if drone=='mesange':
        h0 = -1
    elif drone=='Bernache':
        h0 = 5
    elif drone=='Fulmar':
        h0 = 5
    else:
        h0=0
        print(drone)
        print('Drone unknown')
        
    data = rw_data.read_csv(srtfile)
    n = int(len(data)/6)
    print('number of records : '+str(n))
    record = {}
    for key in ['record_time','date','frame','latitude','longitude','params']:
        record[key]=[]
    for i in range(0,n-1,step):
        #SRT Files contain 6 rows per timestamp
        event = data[i*6:(i+1)*6]
        if int(event[0][0])==i+1:  
            record['record_time'].append(event[1])
            record['date'].append(event[3][0].split(' ')[0])

            time = event[3][0].split(' ')[1][:8]#[:-4]
            h,m,s = time.split(':')
            hnew = str(int(h)+h0)
            t = f"{hnew}:{m}:{s}"            
            record['time'].append(t)
            params = event[4][0]
            if drone=='Fulmar':
                print(params)
                latitude = float(params.split('latitude : ')[1].split(']')[0])
                longitude = float(params.split('longtitude : ')[1].split(']')[0])
            else:
                latitude = float(params.split('latitude: ')[1].split(']')[0])
                longitude = float(params.split('longitude: ')[1].split(']')[0])

            #print(event[3],latitude,longitude)
            record['latitude'].append(latitude)
            record['longitude'].append(longitude)
            record['params'].append(params)
#pprint(d[6:12])
    return record

def main(args):
    convert_flightrecords(args.date)
    
if __name__ =='__main__':
    args = gen_parser()
    main(args)
