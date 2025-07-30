import numpy as np
import pylab as plt
import scipy.interpolate as interp


import icewave.display.graphes as graphes
import icewave.geometry.display as disp
import icewave.field.time as fieldtime

import icewave.tools.rw_data as rw_data

filename = '/Volumes/Backup25/General/Bathymetrie/NONNA_CHS/Bathymetry/NONNA10_4830N06890W.txt'

bathy = rw_data.read_csv(filename,delimiter='\t')


def get_bathy():
    #only for Haha bay now
    filename = df.find_path(disk='Backup25')+'/Bathymetrie/NONNA_CHS/Bathymetry/NONNA10_4830N06890W.txt'
    bathy = rw_data.read_csv(filename,delimiter='\t')

    bathydic = rw_data.csv2dict(bathy)

    bathy = format_bathy(bathydic)
    return bathy

def format_bathy(bathydic):
    Lats=[]
    Longs=[]
    i0=0
    i1=len(bathydic['Depth (m)'])
    step=1

    bathy={}
    Lat = np.asarray(bathydic['Lat (DMS)'][i0:i1:step])
    Long = np.asarray(bathydic['Long (DMS)'][i0:i1:step])
    Depths = np.asarray(bathydic['Depth (m)'][i0:i1:step]).astype(float)

    for (lat,long,depth) in zip(Lat,Long,Depths):
        lat = maps.str2angle(lat)
        long = maps.str2angle(long)
        Lats.append(lat)
        Longs.append(long)
        
    bathy['Lat']=np.asarray(Lats)
    bathy['Long']=np.asarray(Longs)
    bathy['Depth']=np.asarray(Depths)

    return bathy

def interpolator(bathy):
    ###
    #compute an interpolation function of the depth in both sites (Capelan & Haha)
    #return a function of two variabples fi(lon,lat)
    ###
    HahaBox = gps.boxes('haha')
    b1 = gps.check_box(bathy['Long'],bathy['Lat'],HahaBox)
    MercierBox = gps.boxes('capelans')
    b2 = gps.check_box(bathy['Long'],bathy['Lat'],MercierBox)

    b = np.logical_or(b1,b2)
    points = [(lon,lat) for (lon,lat) in zip(bathy['Long'][b],bathy['Lat'][b])]
    depth = bathy['Depth'][b]

    fi = interp.LinearNDInterpolator(points, depth, rescale=False)

    return fi

def display_bathy(bathy,site):
    if site=='mercier':
        BBox = gps.boxes('capelans')
    elif site=='haha':
        BBox = gps.boxes('haha')

    b = gps.check_box(bathy['Long'],bathy['Lat'],BBox)
        
    fig,ax = plt.subplots(figsize=(10,6))
    cm = plt.colormaps['magma']

    sc = ax.scatter(bathy['Long'][b],bathy['Lat'][b],c=bathy['Depth'][b],cmap=cm,marker='o',vmin=-5,vmax=3)
    cbar = plt.colorbar(sc,pad=0.01)#,ticks=[0, 15, 30, 45, 60])
    cbar.set_label('Depth (m)', labelpad=-0.5)

    figs = graphes.legende('Longitude','Latitude',f'Bathymétrie, {site}',cplot=True)
    return ax,figs

def load_tide_data(date):
    marees = {}

    base = df.find_path(disk='Backup25')
#n = len(marees[date]['date'])
    filename = glob.glob(base+date+'/Marees/*.csv')[0]
    print(filename)
    dat=rw_data.read_csv(filename)
    dic=rw_data.csv2dict(dat)
    
    n = len(dic['date'])
    dic['time']=[]
    for i in range(n):
        t0 = dic['date'][i].split('T')[1].split('-')[0]        
        dic['time'].append(tfield.convert_time(t0))
    marees[date]=dic
    marees[date]['tide_height'] = np.asarray(marees[date]['tide_height']).astype(float)
    marees[date]['time'] = np.asarray(marees[date]['time'])
    return marees


def get_local_height(date,t0):
    # return the tide height in Bic Park, the time is given in UTC time
    marees = load_tide_data(date)
    ind = np.argmin(np.abs(marees[date]['time']+6*3600-t0))
    h = marees[datephone]['tide_height'][ind]
    return h
    
