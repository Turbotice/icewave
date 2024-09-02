import numpy as np
import pylab as plt
import glob
import os
import shutil
from pprint import pprint


import icewave.gps.gps as gps
import icewave.tools.datafolders as df

def represent_waypoints(gpx,imin,imax,iplus=[],table=None,ax=None,date=''):
    if table is None:
        table = read_table()

    indices = select(gpx,imin,imax)
    if len(iplus)>0:
        indices = list(indices)+iplus
    Long,Lat = [],[]
    waypoints = np.asarray(gpx.waypoints)[indices]
    for waypoint in waypoints:
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)

    BBox = gps.box_data(Long,Lat,scale=0.7)
    ext = gps.extent(BBox)
    t = gps.tmp_connect()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
    ax,figs = gps.display_map(ext,t,title=date,ax=ax,width=600)

    X,Y = [],[]
    for waypoint in waypoints:
        number = int(waypoint.name[-3:])
        x,y = gps.project(waypoint.longitude,waypoint.latitude)
        X.append(x)
        Y.append(y)
        display(x,y,ax=ax,name=waypoint.name,table=table)
    #savefolder = os.path.dirname(filename)+'/'
    return ax,figs

def represent_table(table,gpx,imin,imax,ax=None):
    indices = select(gpx,imin,imax)
    Long,Lat = [],[]
    waypoints = np.asarray(gpx.waypoints)[indices]
    for waypoint in waypoints:
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)

    BBox = gps.box_data(Long,Lat,scale=0.9)
    ext = gps.extent(BBox)
    t = gps.tmp_connect()
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
    ax,figs = gps.display_map(ext,t,title=date,ax=ax,width=600)

    X,Y = [],[]
    numbers = np.asarray([int(waypoint.name[-3:]) for waypoint in gpx.waypoints])
    for (number,elem) in table:
        if number in numbers:
            ind = np.where(np.asarray(numbers)==number)[0][0]
            waypoint = gpx.waypoints[ind]
            x,y = gps.project(waypoint.longitude,waypoint.latitude)
            X.append(x)
            Y.append(y)
            name = elem
            if '_' in elem:
                tag,num = elem.split('_')
                label = norme[tag]
            else:
                label = norme[elem]
            ax.plot(x,y,label)
            
            #name = name.replace('_','_{')
            #name = r'$'+name+'}$'
            plt.text(x,y-10**(-7),name)
    #print(filename)
    savefolder = os.path.dirname(filename)+'/'
    return ax,figs
    
def display(x,y,ax=None,name='',table=None):
    if table==None:
        label = 'bo'
    else:
        number = int(name[-3:])
        numbers = [tab[0] for tab in table]
        if number in numbers:
            ind = np.where(np.asarray(numbers)==number)[0][0]
            key,elem = table[ind]
            if '_' in elem:
                name = elem
                tag,num = elem.split('_')
                label = norme[tag]
            else:
                label = norme[elem]
            #print(number,table[number],label)
        elif number>max(numbers):
            label='bo'
        else:
            label='bo'
        ax.plot(x,y,label)
        
    #name = name.replace('_','_{')
    #name = r'$'+name+'}$'
    plt.text(x,y-10**(-7),name)
    
def get_day(waypoint):
    return f"{waypoint.time.timetuple().tm_mon:02d}{waypoint.time.timetuple().tm_mday:02d}"

def select(gpx,imin,imax):
    indices = []
    for i,waypoint in enumerate(gpx.waypoints):
        number = int(waypoint.name[-3:])
        #print(number)
        #print(waypoint)
        if number>=imin and number<=imax:#True:#int(waypoint.name)>155 and int(waypoint.name)<250:
            indices.append(i)
    return indices

def read_table(folder):
    filename = glob.glob(folder+'map_table.txt')[0]

    #print(filename)
    with open(filename,'r') as f:
        out = f.read()
    
    lines = out.split('\n')
    table = [line.split('\t') for line in lines]
    dtable = [[int(tab[0]),tab[1]] for tab in table]
    return dtable

def read_norme(path):
    #print(glob.glob(path+'Nomenclature_GPS.txt'))
    filename = glob.glob(path+'Nomenclature_GPS.txt')[0]
    print(filename)
    with open(filename,'r') as f:
        out = f.read()
    
    lines = out.split('\n')
    print(lines)
    tab = [line.split('\t') for line in lines]

    pprint(tab)
    table = np.asarray(tab)

    dtable = {tab[0]:tab[1] for tab in table[1:]}
    #pprint(dtable)
    return dtable

global norme
path = '/Users/stephane/Documents/git/Data_local/Nomenclature/'
norme = read_norme(path)

