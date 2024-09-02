



def gps_time(waypoint):
    time = str(waypoint.time).replace('-','').replace(' ','T').split('+')[0].replace(':','')+'Z'
    return time

def gps_get_times(gpx):
    Times= []
    for waypoint in gpx.waypoints:
        Times.append(gps_time(waypoint))
    return Times

