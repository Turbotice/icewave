
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint
import datetime

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    files = glob.glob(base+date+'/Telephones/*/averages_Summary.csv')
    records = {}
    records['phones'] = {}
    for filename in files:
        record = read_summary(filename)
        for name in record.keys():
            if not name in records['phones'].keys():
                records['phones'][name] = {}
            key = record[name]['name']
            print(key)
            records['phones'][name][key]=record[name]
            #records['phones'][name][key]['path']=filename[nbase:]
    return records

def load_data(record):
    base = df.find_path()
    filename = base + record['path']
    data = rw_data.load_pkl(filename)
    return data
    
def read_summary(filename):
    print(filename)
    data_phone = rw_data.read_csv(filename)
    h0 = -1 #for conversion to UTC
    
    header = data_phone[0]
    keys = header
    phonedict = {}
    records = {}
    for d in data_phone[1:]:
        phone = 'T'+d[0]
        phonedict[phone]={}
        for i,key in enumerate(header):
            try:
                phonedict[phone][key]=d[i]
                incomplete=False
            except:
                incomplete=True
        if incomplete:
            print(f'Header reading failed for {phone}')
            continue
        print(phonedict[phone].keys())
        records[phone]={}
        if 'time_start' in phonedict[phone].keys():
            t0 = phonedict[phone]['time_start'].split(' ')[1][:-4]
            t1 = phonedict[phone]['time_end'].split(' ')[1][:-4]
            times = [t0,t1]
            for j,time in enumerate(times):
                h,m,s = time.split(':')
                hnew = str(int(h)+h0)
                times[j] = f"{hnew}:{m}:{s}"
        else:
            print('No time stamp available')
            times = ['00:00:00','00:00:00']

        records[phone]['time']=times
        lat = float(phonedict[phone]['lat_mean'])
        lon = float(phonedict[phone]['lon_mean'])
        records[phone]['latitude']= [lat,lat]
        records[phone]['longitude']= [lon,lon]
        records[phone]['params']= phonedict[phone]
        records[phone]['name']=phonedict[phone]['name']
        records[phone]['path']=phonedict[phone]['path']
    return records

def get_time_utc(t0):
    tphone,ts = get_time(t0)
    tphone = tphone
    return tphone

def get_time(t0):
    t00, t01 = t0.split(" UTC")[0].split(".")
    date = datetime.datetime.strptime(t00.split(" UTC")[0], "%Y-%m-%d %H:%M:%S")
    tphone = date.timestamp() + int(t01)/1000-3600
    ts = datetime.datetime.fromtimestamp(tphone)
    return tphone,ts
