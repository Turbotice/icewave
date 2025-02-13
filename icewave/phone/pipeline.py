import numpy as np
import pylab as plt
import glob
import os
from pprint import pprint
import scipy.interpolate as interp
import scipy.signal as sig

import icewave.tools.datafolders as df
import icewave.gps.gps as gps
import icewave.display.graphes as graphes
import phonefleet.data as dataphone

from matplotlib.colors import LinearSegmentedColormap as LSC


import argparse

def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate multi instruments data")
    parser.add_argument('-date', dest='date', type=str,default='0211',help='select date to process data')
    parser.add_argument('-step', dest='step', type=int,default=1,help='select Step to be performed')
    parser.add_argument('-num', dest='num', type=int,default=10,help='Select recording number for display phone map')

#    print(parser)   

    args = parser.parse_args()
    print(args)
    return args

global base

def get_folder(date):
    base = df.find_path(disk='Shack25')#f'/media/turbots/BlueDisk/Shack25_local/'
    return base +f'{date}/Phone/'

def get_savefolder(date):
    folder = get_folder(date)+'Results/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def get_phonelist(date):
    folder = get_folder(date)
    folders = glob.glob(folder+'*')

    phonelist=[]
    for f in folders:
        try:
            num = int(os.path.basename(f))
        except:
            continue
        phonelist.append(num)
    print(f'Found {len(phonelist)} phone folders : ')
    print(phonelist)
    return phonelist

def parse_filename(filename):
    name = filename.split('/')[-1]
    c0 = name.split('-')[-1]
    c1 = name.split('-')[-2]
    num = int(name.split('-')[-3])
    typ = name.split('-')[-4]

    if '.' in typ:
        typ = typ.split('.')[-1]

    control = not(c0=='1.csv')
    return num,typ,control

def get_filelist(date,keytest='accelerometer',display=False):
    phonelist = get_phonelist(date)
    folder = get_folder(date)

    folder 
    files = {}
    for phone in phonelist:
        files[phone]={}
        f = folder+str(phone)
        filelist = glob.glob(f+'/*')

        for filename in filelist:
            num,typ,control = parse_filename(filename)
            if not control:
                continue
            if not typ in files[phone].keys():
                files[phone][typ]={}
            files[phone][typ][num]=filename
        if display:
            if keytest in files[phone].keys():
                print(phone,len(files[phone][keytest].keys()))
            else:
                print(phone,f'{keytest} not found')
    return files

def load_lvl_0(files,phone,num,keys=None):
    r={}
    r['num']= num
    r['phone']=phone
    if keys==None:
        keys = files[phone].keys()
    for k in keys:
        if num in files[phone][k]:
            filename = files[phone][k][num]
            dic = dataphone.load_data(filename,phone)
            for key in dic.keys():
                r[key] = dic[key]
        else:
            print(f'No data for phone {phone}, num {num}')
            return None
    return r

def is_moving(r,maxmotion=10,meanmotion=3):
    t = r['tgps']
    X,Y = gps.project(r['gpslon'],r['gpslat'],meter=True)
    d = np.sqrt(X**2+Y**2)    
    #plt.plot(t,d,'k.')
    
    dmax = np.round(np.max(d),decimals=1)
    dmean = np.round(np.mean(d),decimals=1)

    if dmax<maxmotion and dmean<meanmotion:
        move=False
    else:
        move=True
    return move

def position(r):
    X,Y = gps.project(r['gpslon'],r['gpslat'])
    Xmean = np.mean(X)
    Ymean = np.mean(Y)
    return Xmean,Ymean

def display_position(r,ax=None,eps = 10**(-6)):
    X,Y = position(r)
    phone = r['phone']
    ax.plot(X,Y,'rs',markersize=10)
    plt.plot(X,Y,'rs',markersize=10)
    
    sign = int((np.mod(phone,2)-0.5)*2)
    ax.text(X-eps/4,Y-eps*sign,str(phone),color='r',fontsize=20)

def situation_map(files,num):
    phonelist = files.keys()
    fig,ax = plt.subplots(figsize=(20,10))
    ax,figs = gps.display_haha(ax)

    for phone in phonelist:
        r = load_lvl_0(files,phone,num,keys=['gps'])
        if r is not None:
            move = is_moving(r)
            #print(phone,move)
            if move==False:
                display_position(r,ax=ax)
    return fig,ax,figs

def moving_table(files,imin=1,imax=None):
    phonelist = files.keys()
    nphone = len(phonelist)
    imax=0
    
    length=[]
    key = 'accelerometer'
    for phone in phonelist:
        if key in files[phone].keys():
            n = len(files[phone][key].keys())
        else:
            print(f'no accelerometer data for phone {phone}')
            n=0
        length.append(n)
    imax = int(np.nanmedian(length))-3
    print(f'Median number of files : {imax}')
    
    indices = range(imin,imax)
    nt = len(indices)
    moving = np.zeros((nphone,nt))
    for i,phone in enumerate(phonelist):
        for j,num in enumerate(indices):
            r = load_lvl_0(files,phone,num,keys=['gps'])
            if r is not None:
                move = is_moving(r)
                #print(phone,move)       
                moving[i,j]=int(move)
            else:
                moving[i,j]=np.nan
    return moving,indices

def moving_map(moving,indices,phonelist):
    Indices,Phones = np.meshgrid(indices,phonelist)
    fig,ax = plt.subplots(figsize=(15,10))

    colors = [(0, 1, 0), (1, 0, 0)]# R -> G -> B
    cmap_name = 'onoff'
    n_bin = 2
    cmap = LSC.from_list(cmap_name, colors, N=n_bin)

    sc = ax.pcolormesh(Indices,Phones,moving,cmap=cmap)
    #plt.colorbar()
    cbar = fig.colorbar(sc, ax=ax,ticks=[0,1])
    cbar.ax.set_yticklabels(['Immobile','Moving'])

    figs = graphes.legende('number of recording','Phone number','')
    ax.set_xticks(indices)
    ax.set_yticks(phonelist)

    return fig,ax,figs

def N0_to_N1(date,imax=None,num=None,save=True):
    files = get_filelist(date)
    phonelist = list(files.keys())
    savefolder = get_savefolder(date)

    moving,indices = moving_table(files,imax=imax)
    fig,ax,figs = moving_map(moving,indices,phonelist)
    if save:
        graphes.save_figs(figs,savedir=savefolder,prefix='Moving_map_',suffix='_'+date,overwrite=True)

    if num is None:
        num = int(np.mean(indices))
    fig,ax,figs = situation_map(files,num)
    if save:
        graphes.save_figs(figs,savedir=savefolder,prefix='Situation_map',suffix='_'+date,overwrite=True)


if __name__=='__main__':
    args = gen_parser()
    N0_to_N1(args.date,num=args.num,save=True)
