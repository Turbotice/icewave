import stephane.display.graphes as graphes
import stephane.tools.Smath as smath

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

#import sympy #symoblic python
#import mpmath as math
#import cv2
import glob
import csv
import os
import numpy as np

import icewave.phone.rw_pyphone as rw

import icewave.field.drone as drone
import icewave.field.phone as phone
import icewave.field.geophone as geophone
import icewave.field.buoys as buoys

def get_records(date):
    records = drone.get_records(date)
    records.update(phone.get_records(date))
    records.update(geophone.get_records(date))
    records.update(buoys.get_records(date))
    return records

base = '/media/turbots/Hublot24/Share_hublot/Data/'


def get_time_interval(record):
    times = [record[key]['time'] for key in record.keys()]
    t0 = times[0]
    t1 = times[-1]
    return t0,t1

def get_avg_position(record):
    latitudes = [record[key]['latitude'] for key in record.keys()]
    longitudes = [record[key]['longitude'] for key in record.keys()]
    latitude = np.mean(latitudes)
    longitude = np.mean(longitudes)
    return latitude,longitude

