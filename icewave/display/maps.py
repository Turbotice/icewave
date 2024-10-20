import numpy as np

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


