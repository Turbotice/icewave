import numpy as np
import pylab as plt
import glob
import os

#import stephane.display.graphes as graphes
import stephane.display.graphes as graphes

import icewave.gps.garmin as garmin
import fitdecode
import gpxpy


# Import the required library
from geopy.geocoders import Nominatim
import tilemapbase


def tmp_connect():
    """
    Connect to tilemapbase. Beware ! tilemapbase unstable
    INPUT 
    ----- 
    None
    OUTPUT
    -----
    t : tilemapbase object
    """
    tilemapbase.start_logging()
    tilemapbase.init(create=True)
    # Use open street map
    t = tilemapbase.tiles.build_OSM()
    return t

def project(Long,Lat,R=6378000,meter=False):
    """
    Naive local projection of Long & Lat coordinates on a plane 
    INPUT 
    ----- 
    Long : list or np array
    Lat : list or np array
    OUTPUT
    -----
    (x,y) : np arrays
    """
    Long = np.asarray(Long)
    Lat = np.asarray(Lat)
    
    xtile = (Long + 180.0) / 360.0
    lat_rad = np.radians(Lat)
    
    ytile = (1.0 - np.log(np.tan(lat_rad) + (1 / np.cos(lat_rad))) / np.pi) / 2.0
    if meter:
        xtile = xtile*R
        ytile = ytile*R
#        xtile = xtile-xtile[0]
#        ytile = ytile-ytile[0]

    return (xtile, ytile)

def extent(BBox,ratio=1.0):
    """
    #to be revised, weird projection system
    """

    ext = tilemapbase.Extent.from_lonlat(BBox[0],BBox[1],BBox[2],BBox[3])
    ext = ext.to_aspect(ratio)

    print(ext)
    return ext

def display_traj(filename,Long,Lat,save=True):
    title = title_std(filename,Long,Lat)
    ax,t,figs = map_traj(Long,Lat,scale=1,save=False,title=title)
    if save==True:
        savefolder = os.path.dirname(filename)+'/'
        graphes.save_figs(figs,savedir=savefolder)
    return ax,figs

def display_map(ext,t,title='',ax=None,width=600):
    import tilemapbase
    plotter = tilemapbase.Plotter(ext, t, width=width)
    plotter.plot(ax, t)

    figs = {}
    figs.update(graphes.legende('Longitude','Latitude',title,cplot=True,ax=ax))
    return ax,figs

def get_wpts(filename,display=True):
    gpx_file = open(filename, 'r')
    gpx = gpxpy.parse(gpx_file)
    if display:
        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    print('Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation))

        for waypoint in gpx.waypoints:
            print('waypoint {0}, {1} -> ({2},{3})'.format(waypoint.name, waypoint.time.ctime(), waypoint.latitude, waypoint.longitude))
    return gpx

def display_gpx(filename,date,gpx,save=True):
    Long,Lat = [],[]
    for waypoint in gpx.waypoints:
        if True:#int(waypoint.name)>155 and int(waypoint.name)<250:
    #print(int(waypoint.name))
            Long.append(waypoint.longitude)
            Lat.append(waypoint.latitude)
            
    BBox = box_data(Long,Lat,scale=1)
    ext = extent(BBox)
    t = tmp_connect()
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    ax,figs = display_map(ext,t,title=date,ax=ax,width=600)

    X,Y = [],[]
    for waypoint in gpx.waypoints:
        if True:#int(waypoint.name)>155 and int(waypoint.name)<250:
    #if int(waypoint.name) in wpts[key]:
            x,y = project(waypoint.longitude,waypoint.latitude)
            X.append(x)
            Y.append(y)
            plt.text(x,y-2*10**(-7),int(waypoint.name))
    ax.plot(X,Y,'bo')

    if save==True:
        savefolder = os.path.dirname(filename)+'/'
        graphes.save_figs(figs,savedir=savefolder,suffix='_gpx')
    
    return ax,figs

def display_dictwpts(filename,date,wpts,save=True,scale=0.75):
    Long,Lat = [],[]
    for key in wpts.keys():
        X,Y = [],[]
        for waypoint in wpts[key]['wpts']:
            Long.append(waypoint.longitude)
            Lat.append(waypoint.latitude)
            
    BBox = box_data(Long,Lat,scale=scale)
    ext = extent(BBox)
    t = tmp_connect()
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    ax,figs = display_map(ext,t,title=date,ax=ax,width=600)

    for key in wpts.keys():
        X,Y = [],[]
        for waypoint in wpts[key]['wpts']:
            x,y = project(waypoint.longitude,waypoint.latitude)
            X.append(x)
            Y.append(y)
            plt.text(x,y-2*10**(-7),int(waypoint.name))
        ax.plot(X,Y,wpts[key]['label'])

    if save==True:
        savefolder = os.path.dirname(filename)+'/'
        graphes.save_figs(figs,savedir=savefolder,suffix='_Summary')
    return ax,figs
#savefolder = os.path.dirname(filename)+'/'
#graphes.save_figs(figs,savedir=savefolder,prefix='wpts_'+date+'_',suffix='summary',frmt='pdf')


def map_traj(Long,Lat,save=False,scale=0.8,title='',ax=None):
    BBox = box_data(Long,Lat,scale=scale)
    print(BBox)
    X,Y = project(Long,Lat)

    print('toto')
    ext = extent(BBox)

    t = tmp_connect()
    
    if ax==None:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    ax,figs = display_map(ext,t,ax=ax,width=600)
    ax.plot(X,Y,'k-')

    figs.update(graphes.legende('Longitude','Latitude',title,cplot=True))

    if save:
        graphes.save_figs(figs,savedir=savefolder,frmt='png')

    return ax,t,ext,figs


def box_data(Longs,Lats,scale=1.2,width=0.02,square=True):
    """
    return a bounding box scaled around Longs et Lats datapoint
    """
    X,Y = Longs,Lats

    Wlong = (np.max(X)-np.min(X))*scale
    Wlat = (np.max(Y)-np.min(Y))*scale
    Xc = (np.max(X)+np.min(X))/2
    Yc = (np.max(Y)+np.min(Y))/2

    if square:
        W = np.max([Wlong,Wlat])
        Wlong = W
        Wlat = W
    BBox = [Xc-Wlong,Xc+Wlong,Yc-Wlat,Yc+Wlat]
    #bound the BBox
    BBox[0] = np.max([BBox[0],-179.9])
    BBox[1] = np.min([BBox[1],179.9])
    BBox[2] = np.max([BBox[2],-89.9])
    BBox[3] = np.min([BBox[3],89.9])    
    return BBox

def box_data_tight(Longs,Lats):
    BBox = np.zeros(4)
    BBox[0] = np.min(Longs)
    BBox[1] = np.max(Longs)
    BBox[2] = np.min(Lats)
    BBox[3] = np.max(Lats)  
    return BBox

def boxes(name):
    """"
    Return a dictionnary of typical boxes coordinates (for scaled display maps)
    Feel free to add standard boxes here
    """
    d={}
    d['haha'] = [-68.8311,-68.7995,48.3389,48.3551]
    d['bic'] = [-68.8681,-68.7995,48.3181,48.3551]
    d['rimouski'] = [-68.5462,-68.5065,48.4374,48.4636]
    d['capelans'] = [ -68.857421,-68.849667,48.327819,48.332625]#[-68.8625,-68.8330,48.3185,48.3363]

    return d[name]

def check_box(Long,Lat,BBox):
    """
    check box dimensions. Return a boolean value
    """
#    BBox = boxes(name)
    bLong = np.logical_and(Long>BBox[0],Long<BBox[1])
    bLat = np.logical_and(Lat>BBox[2],Lat<BBox[3])
    b = np.logical_and(bLong,bLat)
    #print('Box dimension :'+str(b))
    return b

def get_city(Long,Lat):
    """
    Find city name using Nominatim API
    """
    Longc,Latc = np.mean(Long),np.mean(Lat)
    geolocator = Nominatim(user_agent="MyApp")
    location = geolocator.reverse((Latc,Longc))
    return location

def display_city(Long,Lat):
    location = get_city(Long,Lat)
    print(location)
    if location is None:
        return ''
    parse = location.address.split(',')
    s = parse[2][1:]+','+parse[-1] #fucking arbitrary

    print(s)
    return s

def title_std(filename,Long,Lat):
    date = filename.split('/')[-1].split(' ')[0]
    hour = filename.split('/')[-1].split(' ')[1].split('.fit')[0]
    print(date,hour)

    s = display_city(Long,Lat)
    title = date+', '+hour.replace('.','_')+', '+s
    return title

def display_standard_map(ax,name,title=''):
    if name=='capelans':
        print('Use Capelans map')
        ax,figs = display_mercier(ax,title=title)
    elif name=='haha':
        ax,figs = display_haha(ax,title=title)
    else:
        ax,figs = display_haha(ax,title=title)
    return ax,figs

def display_haha(ax,title='',BBox=None,ratio=1):
    if BBox==None:
        BBox = boxes('haha')
        ext = extent(BBox,ratio=2)
    else:
        ext = extent(BBox,ratio=ratio)
    t = tmp_connect()
    ax,figs = display_map(ext,t,title=title,ax=ax,width=600)
    return ax,figs

def display_mercier(ax,title=''):
    BBox = boxes('capelans')
    ext = extent(BBox,ratio=1)
    t = tmp_connect()
    ax,figs = display_map(ext,t,title=title,ax=ax,width=600)
    return ax,figs
