import numpy as np
import pylab as plt
import glob
import scipy.interpolate as interp


import icewave.display.graphes as graphes
import icewave.geometry.display as disp
import icewave.field.time as tfield

import icewave.tools.rw_data as rw_data
import icewave.tools.datafolders as df
import icewave.gps.gps as gps


import icewave.display.maps as maps


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
    #Note that the new area of interest can be defined in icewave.gps.gps
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
    if site=='capelans':
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

def load_tide_data(date,year):
    marees = {}

    if year=='2025':
        base = df.find_path(disk='Backup25')
    elif year=='2024':
        base = df.find_path(disk='Elements')

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

def display_tide(tide_data,date,year='2025'):
    fig,ax = plt.subplots(figsize=(10,2))
    #manual conversion into UTC time from local time
    ax.plot(tide_data[date]['time']+5*3600,tide_data[date]['tide_height'],'ko')
        
    times = np.linspace(5,29,13)*3600
    #at.set_xlim(18*3600,21*3600)
    #t=[]
    tlist = tfield.display_time(times,precision='min')
    ax.set_xticks(times,tlist)
    figs = graphes.legende('UTC time','Height (m)',f'{date} / {year}')
    return ax,figs

def get_height(date,t0,tide_data=None):
    if tide_data is None:
    # return the tide height in Bic Park,
    # t0 is in UTC time
        tide_data = load_tide_data(date)

    ind = np.argmin(np.abs(tide_data[date]['time']+5*3600-t0))
    h = tide_data[date]['tide_height'][ind]
    return h
    
def get_local_height(lon,lat,date,t0,tide_data=None,bathy_data=None,fi=None):
    htide = get_height(date,t0,tide_data=tide_data)
    if fi is None:
        fi = interpolator(bathy_data)
    hlocal = fi(lon,lat)    
    return htide+hlocal
