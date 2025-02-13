import numpy as np
import pylab as plt
import glob
import os
from pprint import pprint
import scipy.interpolate as interp
import scipy.signal as sig

import icewave.tools.datafolders as df
import icewave.gps.gps as gps
import phonefleet.data as dataphone

global base
base = f'/media/turbots/BlueDisk/Shack25_local/'

def get_folder(date):
    return base +f'Data/{date}/Phone/'

def get_phonelist(date):
    folder = get_folder(date)
    folders = glob.glob(folder+'*')

    phonelist=[]
    for f in folders:
        try:
            num = int(f.split('/')[-1])
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
            dic = dataphone.load_data(filename)
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

def situation_map(files,phonelist,num):
    fig,ax = plt.subplots(figsize=(20,10))
    gps.display_haha(ax)

    for phone in phonelist:
        r = load_lvl_0(files,phone,num,keys=['gps'])
        if r is not None:
            move = is_moving(r)
            #print(phone,move)
            if move==False:
                display_position(r,ax=ax)

def moving_map(phonelist,files):
    nphone = len(phonelist)
    nt = 40
    moving = np.zeros((nphone,nt))
    for i,phone in enumerate(phonelist):
        for j,num in enumerate(range(1,38)):
            r = pl.load_lvl_0(files,phone,num,keys=['gps'])
            if r is not None:
                move = pl.is_moving(r)
                #print(phone,move)       
                moving[i,j]=int(move)
            else:
                moving[i,j]=np.nan
    return moving