
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    files = glob.glob(base+date+'/Telephones/*/averages_Summary.csv')

    records = {}
    records['phones'] = {}
    for filename in files:
        record = read_summary(filename)
        for key in record.keys():
            if not key in records['phones'].keys():
                records['phones'][key]=[record[key]]
            else:
                records['phones'][key].append(record[key])
    return records

def read_summary(filename):
    print(filename)
    data_phone = rw_data.read_csv(filename)
    h0 = -1 #for conversion to UTC
    
    header = data_phone[0]
    keys = header
    phonedict = {}
    records = {}
    for d in data_phone[1:]:
        phone = int(d[0])
        phonedict[phone]={}
        for i,key in enumerate(header):
            phonedict[phone][key]=d[i]

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
        records[phone]['latitude']= phonedict[phone]['lat_mean']
        records[phone]['longitude']= phonedict[phone]['lon_mean']
        records[phone]['params']= phonedict[phone]
    return records
    
