
import scipy.interpolate as interp

import icewave.phone.pipeline as pl

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df

import icewave.field.multi_instruments as multi
import icewave.field.time as timest

import numpy as np


def smooth(data):
    for key in data.keys():
        pass

def load_data(date,phonelist,nums,tmin,tmax,orientation,dt=0.02):
    tmin = multi.convert_time(tmin)
    tmax = multi.convert_time(tmax)

    data = {}
    vars = ['a','g','m']
    names = ['z','y','x'] #z : vertical direction, y: transverse direction to the line, x: longitudinal direction of the line
    keylist = []
    for phone in phonelist:
        print(f'Load data N1 for phone {phone}')
        data[phone]={}
        ts = []
        for var in vars:
            for c,name in zip(orientation,names):
                data[phone][var+'i_'+name]=[]
        
        for num in nums:
            hf = pl.load_lvl_1(date,phone,num)
            t = list(hf['t_sync'])
            t = timest.today_time(t)

            for var in vars:
                for c,name in zip(orientation,names):
                    y_high,y_wave,y_trend,err = filtering(hf['a'+c])
                    f = interp.interp1d(t,y_wave)
                    tinit = t[0]
                    tend = t[-1]
                    ti = np.arange(tinit,tend,dt)
                    yi = f(ti)
                    data[phone][var+'i_'+name] = data[phone][var+'i_'+name] + list(yi)

            ts = ts+t
        ts = np.asarray(ts)
        indices = np.logical_and(ts>=tmin,ts<tmax)
        data[phone]['ts'] =  ts[indices]

        for var in vars:
            for c,name in zip(orientation,names):
                data[phone][var+'i_'+name] = np.asarray(data[phone][var+'i_'+name])[indices]
        save_W2_phone(data,phone,timest.display_time(tmin),timest.display(tmax))
    return data

def save_W2_phone(data,phone,tmin,tmax,index):
    tmintag = tmin.replace(':','_')
    tmaxtag = tmax.replace(':','_')

    folder = dataphone.phone_folder(date,phone)
    filesave = folder+f'Waves_{index}_{tmintag}_{tmaxtag}.h5'
    rw.write_h5(filesave,data)
    print(f'Writing S2 file for phone {phone}: ')
    print(filesave)

def save_W2(data,tmin,tmax,index):
    tmintag = tmin.replace(':','_')
    tmaxtag = tmax.replace(':','_')

    folder = summary_folder(date)
    filesave = folder+f'Waves_{index}_{tmintag}_{tmaxtag}.h5'
    rw.write_h5(filesave,data)
    print('Writing W2 file : ')
    print(filesave)

def filtering(y,fc=0.01,flow=0.0001):
    #correspond to 4Hz and 0.04Hz at fs = 400Hz

    [b1,a1] = sig.butter(4,fc,'high')
    y_high =  sig.filtfilt(b1,a1,y)
    
    [b2,a2] = sig.butter(4,fc,'low')
    y_low =  sig.filtfilt(b2,a2,y)

    [b3,a3] = sig.butter(4,flow,'high')    
    y_wave  =  sig.filtfilt(b3,a3,y_low)

    [b4,a4] = sig.butter(4,flow,'low')
    y_trend  =  sig.filtfilt(b4,a4,y_low)

    sigma = np.std(y)
    err = np.std(y-y_high-y_wave-y_trend)/sigma
    return y_high,y_wave,y_trend,err

def summary_folder(date):
    #base = f'/media/turbots/BlueDisk/Shack25_local/Data/'
    base = df.find_path(year='2025')#
    print(base)
    folder = base +f'{date}/Phone/Summary/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder



