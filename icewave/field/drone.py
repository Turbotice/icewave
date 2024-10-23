
import glob
from pprint import pprint
import numpy as np
import os

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import icewave.drone.drone_timeline as timeline
import icewave.field.Save_extract_record as drone_save

global base
base = df.find_path(disk='Hublot24')
global drones
drones = ['mesange','Bernache','Fulmar']

import argparse
def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate multi instruments data")
    parser.add_argument('-date', dest='date', type=str,default='0226',help='select date to process data')
    parser.add_argument('-step', dest='step', type=int,default=1,help='select step. 1: get_records, 2:convert_flightrecords')

    #parser.add_argument('-step', dest='step', type=int,default=3,help='select Step to be performed')
#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args

def get_record(drone,srtfile)
    nbase = len(base)
    name = srtfile.split('/')[-2]#.split('.')[0]
    record = get_flighrecord(srtfile,drone=drone)
    record['name']=srtfile.split('/')[-1].split('.')[0]
    record['path']=srtfile[nbase:].split('.')[0]
    record['format']='mp4'
    return record,name

def get_records(date,jpg=True):
    srtfiles = get_srtfiles(date)
    
    records = {}
    records['drones']={}
    for key in srtfiles.keys():
        records['drones'][key]={}
        for i,srtfile in enumerate(srtfiles[key]):
            record,name = get_record(key,srtfile)
            print(i,srtfile,name)
            if not name in records['drones'][key]:
                records['drones'][key][name]=[record]
            else:
                records['drones'][key][name].append(record)

    if jpg==True:
        jpgfiles = drone_save.get_jpgfiles(date)
        print(jpgfiles)
        #for now, it does not find any jpg files
        for key in jpgfiles.keys():
            for i,jpgfile in enumerate(jpgfiles[key]):
                name = jpgfile.split('/')[-2]#.split('.')[0]
                print(i,jpgfile,name)
                record = get_jpg_record(jpgfile,drone=key)
                record['format']='jpg'
                if not name in records['drones'][key]:
                    records['drones'][key][name]=[record]
                else:
                    records['drones'][key][name].append(record)
    return records
    
def get_srtfiles(date):
    srtfiles = {}
    print(base)
    for key in drones:
        srt = glob.glob(base+date+'/Drones/'+key+'/*/*.SRT')#/*/*.srt')
        pprint(srt)
        if len(srt)>0:
            srtfiles[key] = srt
        else:
            print(f"No data for {key} on {date}")
    return srtfiles

def get_mp4files(date,save=True):
    import cv2
    import os

    for key in drones:
        mp4files = glob.glob(base+date+'/Drones/'+key+'/*/*.MP4')#/*/*.srt')
        for filename in mp4files:
            if save:
                save_mp4file(drone,filename)
    return mp4files

def  generate_flightrecords(date):
    srtfiles = get_srtfiles(date)
    

def save_mp4file(drone,filename):
    print(key,filename.split('/')[-2])
    cam = cv2.VideoCapture(filename)
    ret,frame = cam.read()
    imagefile = filename.split('.')[0]+'_exemple.tiff'

    print(f"Save image : {imagefile.split('/')[-1]}")
    cv2.imwrite(imagefile, frame) # Save the image

def convert_flightrecords(date):
    csvfiles = get_csvfiles(date)
    records={}
    records['drones']={}
    for drone in csvfiles.keys():
        records['drones'][drone]={}
        for i,csvfile in enumerate(csvfiles[drone]):
            print(csvfile)
            record=parse_csv_flightrecord(csvfile,drone=drone)
            if i==0:
                records['drones'][drone]=record
            else:
                for key in records['drones'][drone].keys():
                     records['drones'][drone][key]=list(records['drones'][drone][key])+list(record[key])
        savedict = records['drones'][drone]
        times = np.asarray([timeline.to_UTC(s,h0=0) for s in savedict['CUSTOM.updateTime [local]']])
        indices = np.argsort(times)
    
        for key in record.keys():
            savedict[key]=np.asarray(savedict[key])[indices]    

        print('save pickle')
        filename = os.path.dirname(csvfile)+'/Flightrecord_dict.pkl'
        rw_data.write_pkl(filename,savedict)
        print('flightrecord saved in pickle format')

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
    keys_date = ['CUSTOM.date [local]', 'CUSTOM.updateTime [local]','OSD.flyTime', 'OSD.flyTime [s]']
    keys_bool = ['CAMERA.isPhoto', 'CAMERA.isVideo']
    keys0 = ['OSD.latitude', 'OSD.longitude', 'OSD.height [ft]']
    keys1 = ['OSD.altitude [ft]', 'OSD.mileage [ft]', 'OSD.hSpeed [MPH]', 'OSD.xSpeed [MPH]', 'OSD.ySpeed [MPH]', 'OSD.zSpeed [MPH]']
    keys2 = ['OSD.pitch', 'OSD.roll', 'OSD.yaw', 'OSD.yaw [360]', 'OSD.gpsNum']
    keys3 = ['GIMBAL.pitch', 'GIMBAL.roll', 'GIMBAL.yaw', 'GIMBAL.yaw [360]']            
    keys_float = keys0+keys1+keys2+keys3

    record = {}
    for key in keys_float:
#        print(key,data[key][1000])
        record[key]= [float(d) for d in data[key]]
    for key in keys_date:
        record[key]= data[key]
    for key in keys_bool:
        record[key]= [bool(d) for d in data[key]]
    return record

def cut_flightrecord(record,flight):# at that stage, all files should already by in UTC time
    import icewave.field.time as fieldtime

    tinit = timeline.to_UTC(record['time'][0],h0=0)
    tend = timeline.to_UTC(record['time'][-1],h0=0)

    times = np.asarray([timeline.to_UTC(s,h0=0) for s in flight['CUSTOM.updateTime [local]']])
    iinit = np.argmin(np.abs(times-tinit))
    iend = np.argmin(np.abs(times-tend))
    #to_UTC(string,h0=-1)
    print(iinit,iend)
    print(fieldtime.display_time([tinit,tend]))#,tend)
    #print(flight['CUSTOM.updateTime [local]'][idx])

    flight_p = {}
    for key in flight.keys():
        flight_p[key]=flight[key][iinit:iend]
    return flight_p
    

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
    for key in ['record_time','time','date','frame','latitude','longitude','params']:
        record[key]=[]
        
    indlist = list(range(0,n-1,step))
    if indlist[-1]<n-2:
        indlist=indlist+[n-2]
    for i in indlist:
        #SRT Files contain 6 rows per timestamp
        event = data[i*6:(i+1)*6]
        if int(event[0][0])==i+1:
            record['frame'].append(i) #starting at frame 0
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
    if args.step==1:
        get_records(args.date)
    if args.step==2:
        convert_flightrecords(args.date)
    if args.step==3:
        get_mp4files(args.date)
#    convert_flightrecords(args.date)
    
if __name__ =='__main__':
    args = gen_parser()
    main(args)
