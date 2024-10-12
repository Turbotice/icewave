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

import icewave.phone.rw_pyphone as rw




base = '/media/turbots/Hublot24/Share_hublot/Data/'


def get_flighrecord(srtfile):
    data = rw_data.read_csv(srtfile)

    n = int(len(data)/6)
    print('number of records : '+str(n))
    record = {}
    for i in range(0,n-1,100):
        event = data[i*6:(i+1)*6]
        if int(event[0][0])==i+1:  
            record[i]={}
            record[i]['time']=event[1]
            record[i]['date']=event[3]
            params =event[4][0]

            latitude = float(params.split('latitude: ')[1].split(']')[0])
            longitude = float(params.split('longitude: ')[1].split(']')[0])

            #print(event[3],latitude,longitude)
            record[i]['latitude']=latitude
            record[i]['longitude']=longitude
            record[i]['params'] = params
#pprint(d[6:12])
    return record



