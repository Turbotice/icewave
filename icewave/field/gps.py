
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
    filename = glob.glob(folder+'*.gpx')[0]#/*/*.srt')
    with open(filename,'r') as f:
        gpx = gpxpy.parse(f)

    record = tables.dict_from_gpx(gpx,folder)
    records={}
    records['gps']={}
    records['gps']['garmin_sp']=record
    return record
