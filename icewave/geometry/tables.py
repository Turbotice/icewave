import numpy as np
import pylab as plt
import glob
import os
import shutil
from pprint import pprint

import icewave.gps.gps as gps
import icewave.tools.datafolders as df
import icewave.geometry.tables as tables

global norme_folder
norme_folder = '/Users/stephane/Documents/git/Data_local/Nomenclature/'


def read_norme(path):
    #print(glob.glob(path+'Nomenclature_GPS.txt'))
    filename = glob.glob(path+'Nomenclature_GPS.txt')[0]
    #print(filename)
    with open(filename,'r') as f:
        out = f.read()
    
    lines = out.split('\n')
    #print(lines)
    tab = [line.split('\t') for line in lines]

    #pprint(tab)
    #pprint(tab)
    table = np.asarray(tab)
    keys = table[0,1:]
    dtable = {tab[0]:tab[1] for tab in table[1:]}
    details ={}
    for i,t in enumerate(dtable.keys()):
        details[t]={}
        for j,key in enumerate(keys):       
            details[t][key] = table[i+1,j+1]    
    return dtable,details

def dict_from_gpx(gpx,folder,reg=''):
    data={}
    table = read_table(folder,reg=reg)

    print(table)
    tabledict = table_2dict(table)

    pprint(tabledict.keys())
    imin = np.min(list(tabledict.keys()))
    imax = np.max(list(tabledict.keys()))
    print(imin,imax)
    
    indices = select(gpx,imin,imax)
    waypoints = np.asarray(gpx.waypoints)[indices]
    
    for key in tabledict.keys():
        found=False
        for waypoint in waypoints:
            #name = waypoint.name
            number = get_number(waypoint)
            print(number,key)
            if int(number)==int(key):
                print(number)
                print('found!')
                found=True
                tag = tabledict[key]['tag']
                data[tag]={}
                data[tag]['longitude']=waypoint.longitude
                data[tag]['latitude']=waypoint.latitude
                data[tag]['elevation']=waypoint.elevation
                data[tag]['time']=waypoint.time
                if 'value' in tabledict[number].keys():
                    data[tag]['value']=tabledict[number]['value']
        if not found:
            print(f'Waypoint {key} not found in the table')
    return data

def represent_waypoints(gpx,imin,imax,iplus=[],imoins=[],table=None,ax=None,date='',scale=0.7):
    if table is None:
        table = read_table()

    indices = select(gpx,imin,imax)
    if len(iplus)>0:
        indices = list(indices)+iplus
    if len(imoins)>0:
        pass
    Long,Lat = [],[]
    waypoints = np.asarray(gpx.waypoints)[indices]
    for waypoint in waypoints:
        Long.append(waypoint.longitude)
        Lat.append(waypoint.latitude)

    BBox = gps.box_data(Long,Lat,scale=scale)
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
    norme,details = read_norme(norme_folder)

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
                print(elem)
                label = norme[elem]
            ax.plot(x,y,label)
            
            #name = name.replace('_','_{')
            #name = r'$'+name+'}$'
            plt.text(x,y-10**(-7),name)
    #print(filename)
    savefolder = os.path.dirname(filename)+'/'
    return ax,figs
    
def display(x,y,ax=None,name='',table=None):
    norme,details = read_norme(norme_folder)
    
    if table==None:
        label = 'bo'
    else:
        number = int(name[-3:])
        numbers = [tab[0] for tab in table]
        if number in numbers:
            ind = np.where(np.asarray(numbers)==number)[0][0]
            print(table[ind])
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

def get_number(waypoint):
    try:
        return int(waypoint.name[-3:])
    except:
        print(f'cannot read {waypoint.name}')
        return 0

def select(gpx,imin,imax):
    indices = []
    for i,waypoint in enumerate(gpx.waypoints):
        number = get_number(waypoint)
        #print(number)
        #print(waypoint)
        if number>=imin and number<=imax:#True:#int(waypoint.name)>155 and int(waypoint.name)<250:
            indices.append(i)
    return indices

def read_table(folder,reg=''):
    filelist = glob.glob(folder+'Map_Table*'+reg+'.txt')
    print(filelist)
    filename = filelist[0]
    if len(filelist)>1:
        print(f'Warning : several map tables found in {folder}')
    #filename = glob.glob(folder+'Map_Table.txt')[0]

    #print(filename)
    with open(filename,'r') as f:
        out = f.read()
    
    lines = out.split('\n')
    table = [line.split('\t') for line in lines]
    dtable = [[int(tab[0])]+tab[1:] for tab in table]
    return dtable

def table_2dict(table):
    data = {}
    for tab in table:
        key = tab[0]
        data[tab[0]]={}
        data[tab[0]]['tag']=tab[1]
        if len(tab)==3:
            print(f"Value {tab[2]} associated to waypoint {tab[0]} for tag {tab[1]}")
            data[tab[0]]['value']=tab[2]
    return data

def rename_waypoints(gpx,table,root='BIC25'):
    names = [wpt.name.split('_')[0] for wpt in gpx.waypoints]
    for item in table:
        num = format(item[0], '04d')
        target = root+num
        if target in names:
            ind = names.index(target)
            gpx.waypoints[ind].name = gpx.waypoints[ind].name + '_'+item[1]
            #print(gpx.waypoints[ind].name)
    return gpx

def change_symbols(gpx,table,root='BIC25'):
    names = [wpt.name.split('_')[0] for wpt in gpx.waypoints]
    norme,details = read_norme(norme_folder)

    for item in table:
        num = format(item[0], '04d')
        target = root+num
        if target in names:
            ind = names.index(target)
            key = item[1].split('_')[0]
            gpx.waypoints[ind].symbol = details[key]['Garmin_flag']
    return gpx

def save_gpx(filegpx,gpx):
    filename = filegpx.split('.')[0]+'_edited.gpx'
    with open(filename, "w") as f:
        f.write(gpx.to_xml())
    
#global norme
#path = '/Users/stephane/Documents/git/Data_local/Nomenclature/'
#norme = read_norme(path)

