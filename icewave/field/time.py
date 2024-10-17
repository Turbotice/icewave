

import datetime

def gps_time(waypoint):
    time = str(waypoint.time).replace('-','').replace(' ','T').split('+')[0].replace(':','')+'Z'
    return time

def gps_get_times(gpx):
    Times= []
    for waypoint in gpx.waypoints:
        Times.append(gps_time(waypoint))
    return Times

def today_time(times):
    t0 = times[0]
    tstr = str(datetime.datetime.fromtimestamp(t0))
    tref = tstr.split(' ')[0]+' 00:00:00'
    tref = datetime.datetime.strptime(tref, '%Y-%m-%d %H:%M:%S')

    tdtimes = []
    for t in times:
        t =  datetime.datetime.fromtimestamp(t)#datetime_object = datetime.strptime(datetime_str, '%m/%d/%y %H:%M:%S')
        tdtimes.append((t - tref).total_seconds())
    return tdtimes

def display_time(times):
    times = today_time(times)
    
    timestamps=[]
    for t in times:
        hour = int(t/3600)
        minute = int((t-hour*3600)/60)
        sec = int(t-hour*3600-minute*60)
        if sec<10:
            sec = '0'+str(sec)
        timestamps.append(f'{hour}:{minute}:{sec}')
    return timestamps
