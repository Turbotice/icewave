
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    files = glob.glob(base+date+'/Telephones/*/averages_Summary.csv')

    records = {}
    records['phone'] = {}
    for filename in files:
        record = read_summary(filename)
        for key in record.keys():
            if not key in records['phone'].keys():
                records['phone'][key]=[record[key]]
            else:
                records['phone'][key].append(record[key])
    return records

def read_summary(filename):
    print(filename)
    data_phone = rw_data.read_csv(filename)

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
        t0 = phonedict[phone]['time_start'].split(' ')[1][:-4]
        t1 = phonedict[phone]['time_end'].split(' ')[1][:-4]
        records[phone]['time']= [t0,t1]
        records[phone]['latitude']= phonedict[phone]['lat_mean']
        records[phone]['longitude']= phonedict[phone]['lon_mean']
        records[phone]['params']= phonedict[phone]
    return records
    
