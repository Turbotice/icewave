import pylab as plt
import numpy as np

#import icewave.field.multi_instruments as multi
#import icewave.geometry.display as display
#import icewave.display.graphes as graphes

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

import icewave.gps.gps as gps

def get_base(disk='Elements',year='2024'):
    return df.find_path(disk=disk,year=year)

def get_record(year='2024'):
    base = get_base(year=year)
    filename =  base+f'{date}/Summary/records_{date}.pkl'
    return rw_data.load_pkl(filename)


def display(records,'all'):

    fig,ax = plt.subplots(figsize=(10,10))
    gps.display_haha(ax)

    for drone in records['drones'].keys():
        rec = records['drones'][drone]
        for name in rec.keys():
            print(name)
            lat = rec[name][0]['latitude'][0]
            lon = rec[name][0]['longitude'][0]

            X,Y = gps.project(lon,lat)
            ax.plot(X,Y,'ko')
            
    markers = {'geophones':'gv','phones':'rs','buoys':'mo'}
    for instru in markers.keys():
        marker = markers[instru]
        for key in records[instru].keys():
            rec = records[instru][key]
            for name in rec.keys():
                print(name)
                print(rec[name].keys())
                lat = rec[name]['latitude']
                lon = rec[name]['longitude']

                X,Y = gps.project(lon,lat)
                ax.plot(X,Y,marker)
            
                      

            
