import pylab as plt
import numpy as np

#import icewave.field.multi_instruments as multi
#import icewave.geometry.display as display
#import icewave.display.graphes as graphes

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

import icewave.display.graphes as graphes
import icewave.gps.gps as gps

import timeline


def get_base(disk='Elements',year='2024'):
    return df.find_path(disk=disk,year=year)

def get_record(date,year='2024',disk='Elements'):
    base = get_base(year=year,disk=disk)
    filename =  base+f'{date}/Summary/records_{date}.pkl'
    return rw_data.load_pkl(filename)

def plot(lon,lat,BBox=None,text=None,ax=None,marker='ko',project=True):
    if BBox is not None:
        b = gps.check_box(np.asarray(lon),np.asarray(lat),BBox)
    else:
        #print(np.asarray(lon).shape)
        b = np.ones(np.asarray(lon).shape,dtype=bool)
        
    if project:
        X,Y = gps.project(lon,lat)
    else:
        X,Y = np.asarray(lon),np.asarray(lat)
    ax.plot(X[b],Y[b],marker)
    if text is not None and np.sum(b)>0:
        ax.text(X[b],Y[b],text)

def display_map(records,ax=None,BBox=None,project=True):
    if ax==None:
        fig,ax = plt.subplots(figsize=(10,10))
        gps.display_haha(ax,BBox=BBox)
        #print('display')
        
    for drone in records['drones'].keys():
        rec = records['drones'][drone]
        for name in rec.keys():
            #print(name)
            lat = rec[name][0]['latitude'][0]
            lon = rec[name][0]['longitude'][0]
            plot(lon,lat,BBox=BBox,text=name,ax=ax,marker='ko',project=project)
            
    markers = {'geophones':'gv','phones':'rs','buoys':'mo'}
    for instru in markers.keys():
        marker = markers[instru]
        for key in records[instru].keys():
            rec = records[instru][key]
            for name in rec.keys():
                #print(name)
                #print(rec[name].keys())
                lat = rec[name]['latitude']
                lon = rec[name]['longitude']

                plot(lon,lat,BBox=BBox,text=None,ax=ax,marker=marker,project=project)            
    return ax            

def portfolio(date,images,drone):
    drone = 'mesange'
    nombre = len(list(images[drone].keys()))
    if nombre<6:
        nc=1
    elif nombre<9:
        nc=2
    else:
        nc=3
    nr=int(np.ceil(nombre/nc))

    fig,axs = plt.subplots(figsize=(nc*6.5,nr*4),nrows=nr,ncols=nc)
    axs = np.reshape(axs,(nr*nc,1))
    for key,ax in zip(images[drone].keys(),axs):
        ax = ax[0]
        im = images[drone][key]
        ax.imshow(im)
        figs = graphes.legende('','',key,ax=ax)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.01)
    #prefix = f'{drone}_portfolio_{date}'
    return figs,ax

def synthesis(records,date='0226',project=False):
    #manual BBox, to be improved
    BBox = [-68.8238994794379, -68.80561255147802, 48.33896404326063, 48.35725097122051]
    fig = plt.figure(figsize=(15,8))
    nx = 3
    ny = 3

    axmap = plt.subplot2grid((nx, ny), (0, 0),colspan=3,rowspan=2)

    #gps.display_haha(axmap,BBox=BBox,ratio=1.5)
    display_map(records,BBox=BBox,ax=axmap,project=project)

    #selection = ['08_waves_003','14-waves_007','22-waves_011']#,'20-waves_010','24-waves_013']
    axtime = plt.subplot2grid((nx, ny), (nx-1, 0),colspan=ny)
#axim.append(plt.subplot2grid((nx, ny), (1, 0),colspan=2))
#axim.append(plt.subplot2grid((nx, ny), (1, 2),colspan=2))
    figs = timeline.display_records(records,date,ax=axtime)

    axmap.axis(BBox)
    axmap.axis('equal')
#graphes.save_figs(figs,savedir=savefolder_local,prefix=f'{date}_situationmap_overview_2',overwrite=True,frmt='png')
    return fig,axmap,axtime

def exemple_portfolio():
    date = '0226'

            
