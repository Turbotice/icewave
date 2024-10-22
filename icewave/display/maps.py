import numpy as np
import pylab as plt

import icewave.display.graphes as graphes

def display_map(measures,remote=True,w=10,h=10):
    plt,ax = plt.subplots(figsize=(w,h))
    for measure in measures:
        ax.plot(measure['longitude'],measure['latitude'],measure['label'])

    [Lonmin,Lonmax] = ax.get_xlim()
    [Latmin,Latmax] = ax.get_ylim()

    deg,minute,sec = get_range(Lonmin,Lonmax)
    lon_ticks = toangle(deg,minute,sec,sign=-1)
    deg,minute,sec = get_range(Latmin,Latmax)
    lat_ticks = toangle(deg,minute,sec,sign=1)
                
    ax.set_xticks(lon_ticks,display_longitude(lon_ticks))
    ax.set_yticks(lat_ticks,display_latitude(lat_ticks))

    figs = graphes.legende('Longitude','Latitude','')


        
def get_range(anglemin,anglemax):
    degm,minutem,secm = angle2coord(anglemin)
    degM,minuteM,secM = angle2coord(anglemax)

    if degm==degM:
        if minutem==minuteM:
            n = secM-secm+1
            secs = np.arange(secm,secM+1,n)
            return degm,minutem,secs
        else:
            n = minuteM-minutem+1
            minutes = np.arange(minutem,minuteM+1,n)
            return degm,minutes,0
    else:
        n = degM-degm+1
        degs = np.arange(degm,degM+1,n)
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

def coord2angle(deg,minute,sec,sign=1):
    return sign*(deg*3600+minute*60+sec)/3600

def angle2coord(angle):
    if angle<0:
        angle=-angle
    deg = int(np.floor(angle))
    minute = int(np.floor((angle-deg)*60))
    sec = int(np.floor(((angle-deg)*60-minute)*60))
    return deg,minute,sec

def display_longitude(deg,minute,sec):
    if deg<0:
        s=str(-deg)+r'$^\circ$W'
    else:
        s=str(-deg)+r'$^\circ$E'
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def display_latitude(deg,minute,sec):
    if deg<0:
        s=str(-deg)+r'$^\circ$S'
    else:
        s=str(-deg)+r'$^\circ$N'
    s = s+' '+str(minute)+'\''+' '+str(sec)+'\"'
    return s

def display_longitudes(Longitudes):
    return [convert_longitude(Lon) for Lon in Longitudes] 
def display_latitudes(Latitudes):
    return [convert_latitude(Lat) for Lat in Latitudes] 


