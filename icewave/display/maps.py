import numpy as np
import pylab as plt

import icewave.display.graphes as graphes
import icewave.geometry.display as disp
import icewave.field.time as fieldtime

import icewave.gps.gps as gps

import argparse


def display_map(measures,remote=True,w=10,h=10,ax=None):
    b = 10**(-7)
    if ax==None:
        fig,ax = plt.subplots(figsize=(w,h))

    Lons,Lats=[],[]
    
    for key in measures.keys():
        m = measures[key]
        instrument = key.split('_')[0]
        sensor = key.split('_')[1]
        #print(instrument,sensor)
        if sensor=='mesange':
            latshiftlabel = -b
        else:
            latshiftlabel = b
            
        label = disp.colors[instrument]
        if instrument=='drones':
            if key.split('_')[-1]=='0':
                text = '_'.join(key.split('_')[1:-1])
            else:
                text = ''
        elif instrument=='gps':
            text=sensor
            latshiftlabel = -b
        elif instrument=='buoys':
            
            if True:#m['recording'].split('_')[-1]=='1700':
                text=sensor
                latshiftlabel = b*1.5
            else:
                text=''
        else:
            text=sensor
        if instrument=='gps':
            k = m['name'][0].split('_')[0]
            label = disp.labels[k]
            text= m['name'][0][0]+m['name'][0][-2:]
#            label = disp.colors[
        excludes = ['buoys_B5_buoy5_sbg_20240223_1700_0']
        name = key.split('_')[1]
        if name=='B5':# or instrument=='phones':# or instrument=='gps':
            continue
        if len(m['longitude'])>0:
            lon = m['longitude'][0]
            lat = m['latitude'][0]
            Lons.append(lon)
            Lats.append(lat)
            
            X,Y = gps.project(lon,lat)

            print(Y)
            ax.plot(X,Y,label)
            ax.text(X,Y+latshiftlabel,text,color='k')
        else:
            print(f'No GPS coordinates for {key}')

    Lonmin = np.min(Lons)
    Lonmax = np.max(Lons)
    Latmin = np.min(Lats)
    Latmax = np.max(Lats)
#    [Lonmin,Lonmax] = ax.get_xlim()
#    [Latmin,Latmax] = ax.get_ylim()

    print(Lonmin,Lonmax,Latmin,Latmax)
    deg,minute,sec = get_range(Lonmax,Lonmin)
    lon_display = get_range_display(deg,minute,sec)    
    lon_ticks = coord2angle(deg,minute,sec,sign=-1)
    
    deg,minute,sec = get_range(Latmin,Latmax)
    lat_display = get_range_display(deg,minute,sec)    
    lat_ticks = coord2angle(deg,minute,sec,sign=1)
#    lat_display = [display_latitude(angle2coord(lat)) for lat in lat_ticks]

    X,_ = gps.project(lon_ticks,np.zeros(len(lon_ticks)))
    _,Y = gps.project(np.zeros(len(lat_ticks)),lat_ticks)
    print(X,Y)
    ax.set_xticks(X,lon_display)
    ax.set_yticks(Y,lat_display)
    #ax.set_xticks(lon_ticks,lon_display)
    #ax.set_yticks(lat_ticks,lat_display)

    #plt.axis('equal')
    figs = graphes.legende('Longitude','Latitude','')
    return figs

def get_range_display(deg,minute,sec):
    #print(deg,minute,sec)
    print(deg,minute,sec)
    print(type(minute))
    if type(sec)==list or type(sec)==np.ndarray:
        lon_display = [display_longitude((deg,minute,s)) for s in sec]#(d,m,s in lon_ticks]
    elif type(minute)==list or type(minute)==np.ndarray:
        lon_display = [display_longitude((deg,m,sec)) for m in minute]#(d,m,s in lon_ticks]
    elif type(deg)==list or type(deg)==np.ndarray:
        lon_display = [display_longitude((d,minute,sec)) for d in deg]#(d,m,s in lon_ticks]
    else:
        print('No spatial coordinate range available')
        lon_display=None
    print(lon_display)
    return lon_display
    
def get_measures(records,tmin,tmax):
    measures = {}

    ti = fieldtime.convert_time(tmin)
    tf = fieldtime.convert_time(tmax)
    trange = [ti,tf]
    for instrument in records.keys():
        for sensor in records[instrument].keys():
            for k in records[instrument][sensor].keys():
                record = records[instrument][sensor][k]
                if not type(record)==list:
                    record = [record]
                for i,rec in enumerate(record):
                    #print(instrument,sensor,k,i)
                    ti = rec['time'][0]
                    tf = rec['time'][-1]
                    #print(ti,tf)
                    ti = fieldtime.convert_time(ti)
                    tf = fieldtime.convert_time(tf)

                    if tf>trange[0] and ti<trange[1]:
                        name = '_'.join([instrument,sensor,k,str(i)])
                        measure={}
                        for key in rec.keys():
                            measure[key]=rec[key]
                        measure['instrument']=instrument
                        measure['sensor']=sensor
                        measure['recording']=k
                        measure['num']=i
                        measures[name]=measure
    return measures
    
def get_range(anglemin,anglemax):
    degm,minutem,secm = angle2coord(anglemin)
    degM,minuteM,secM = angle2coord(anglemax)

    #print(angle2coord(anglemin))
    #print(angle2coord(anglemax))
    if degm==degM:
        if minutem==minuteM:
            n = secM-secm+1
            step = int(np.ceil(n/8))
            secs = np.arange(secm,secM+1,step)
            #print(secs)
            return degm,minutem,secs
        else:
            n = minuteM-minutem+1
            minutes = np.arange(minutem,minuteM+1)
            return degm,minutes,0
    else:
        n = degM-degm+1
        degs = np.arange(degm,degM+1)
        return degs,0,0
    
def convert_longitude(lon):
    if lon<0:
        lon=-lon
        deg = int(np.floor(lon))
        s=str(deg)+r'$^\circ$W'
    else:
        deg = int(np.floor(lon))
        s=str(deg)+r'$^\circ$E'
    minute = int(np.floor((lon-deg)*60))
    sec = int(np.floor(((lon-deg)*60-minute)*60))
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def convert_latitude(lat):
    if lat<0:
        lat=-lat
        deg = int(np.floor(lat))
        s=str(deg)+r'$^\circ$S'
    else:
        deg = int(np.floor(lat))
        s=str(deg)+r'$^\circ$N'
    minute = int(np.floor((lat-deg)*60))
    sec = int(np.floor(((lat-deg)*60-minute)*60))
#    print((lon-deg)*60-minute)
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def str2angle(s):
    if s[-1]=='S' or s[-1]=='W':
        sign=-1
    else:
        sign=1
    deg,minute,sec = np.asarray(s[:-1].split('-')).astype(float)
    return coord2angle(deg,minute,sec,sign=sign)
    
def coord2angle(deg,minute,sec,sign=1):
    return sign*(np.abs(deg)*3600+minute*60+sec)/3600
    
def angle2coord(angle):
    if angle<0:
        angle=-angle
        sign=-1
    else:
        sign=1
    deg = int(np.floor(angle))
    minute = int(np.floor((angle-deg)*60))
    sec = int(np.floor(((angle-deg)*60-minute)*60))

    if sign==-1:
        deg = -deg
    return deg,minute,sec

def display_longitude(tup):
    deg,minute,sec = tup
    if deg<0:
        s=str(-deg)+r'$^\circ$W'
    else:
        s=str(deg)+r'$^\circ$E'
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def display_latitude(tup):
    deg,minute,sec = tup
    if deg<0:
        s=str(-deg)+r'$^\circ$S'
    else:
        s=str(deg)+r'$^\circ$N'
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def display_longitudes(Longitudes):
    return [convert_longitude(Lon) for Lon in Longitudes] 
def display_latitudes(Latitudes):
    return [convert_latitude(Lat) for Lat in Latitudes]
