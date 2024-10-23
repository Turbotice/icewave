
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

import icewave.geometry.tables as tables

global base
base = df.find_path(disk='Hublot24')
import gpxpy

def get_records(date,year='2024'):
    folder = base+date+'/GPS/'
    filelist = glob.glob(folder+'*.gpx')
    if len(filelist)==0:
        print('No GPS data for this day')
        return {}
    if len(filelist)>1:
        print('Warning : several gpx files found')
    filename=filelist[0]

    with open(filename,'r') as f:
        gpx = gpxpy.parse(f)

    record = get_record_fromgpx(gpx,folder)
    records={}
    records['gps']={}
    records['gps']['garminSP']=record

    return records


def get_record_fromgpx(gpx,folder):
    record = tables.dict_from_gpx(gpx,folder)
    #
    for key in record.keys():
        record[key]['name']=key
        record[key]['time'] = str(record[key]['time']).split(' ')[1].split('+')[0]
    for key in record.keys():
        for k in record[key]:
            record[key][k]=[record[key][k]]
    return record
