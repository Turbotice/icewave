
import scipy.interpolate as interp

import icewave.phone.pipeline as pl

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df

import icewave.field.multi_instruments as multi
import icewave.field.time as timest
import scipy.signal as sig

import numpy as np
import os


def smooth(data):
    for key in data.keys():
        pass

def load_data(date,phonelist,nums,tmin,tmax,orientation,dt=0.02):
    tmin = multi.convert_time(tmin)
    tmax = multi.convert_time(tmax)

    data = {}
    keylist = []
    for phone in phonelist:
        data[phone]=load_data_phone(date,phone,nums,tmin,tmax,orientation,dt=dt)
        #save_W2_phone(data,date,phone,tmin,tmax,nums[0])
    return data

def load_data_phone(date,phone,nums,tmin,tmax,orientation,dt=0.02):
    print(f'Load data N1 for phone {phone}')
    variables = ['a','g','m']
    names = orientation #z : vertical direction, y: transverse direction to the line, x: longitudinal direction of the line

    data={}
    buffer = {}
    for var in variables:
        buffer['t'+var]=[]
        for c,name in zip(orientation,names):
            buffer[var+name]=[]
        
    for num in nums:
        hf = pl.load_lvl_1(date,phone,num)
        if hf is not None:
            t = list(hf['t_sync'])
            t = timest.today_time(t)
            tlag = hf['t_sync'][0]-hf['ta'][0]
                
            for var in variables:
                tvar = np.asarray(hf['t'+var])+tlag
                tinit = tvar[0]
                tend = tvar[-1]
                ti = np.arange(tinit,tend,dt)
                buffer['t'+var] = buffer['t'+var] + list(ti)
                for c,name in zip(orientation,names):
                    y_high,y_wave,y_trend,err = filtering(hf[var+c])    
                    f = interp.interp1d(tvar,y_wave)
                    yi = f(ti)
                    buffer[var+name] = list(buffer[var+name]) + list(yi)
    for var in variables:
        tsi = np.asarray(buffer['t'+var])
        tsi = np.asarray(timest.today_time(tsi))
        print(tsi)
        indices = np.logical_and(tsi>=tmin,tsi<tmax)
        data['ti'+var] =  tsi[indices]
        for c,name in zip(orientation,names):
            data[var+'i_'+name] = np.asarray(buffer[var+name])[indices]
    return data

def save_W2_phone(data,date,phone,tmin,tmax,index):
    tmin,tmax = timest.display_time([tmin,tmax])
    tmintag = tmin.replace(':','_')
    tmaxtag = tmax.replace(':','_')

    folder = pl.get_folder(date)+'Waves/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    filesave = folder+f'Waves_{phone}_{index}_{tmintag}_{tmaxtag}.h5'
    
    print(f'Writing S2 file for phone {phone}: ')
    print(filesave)
    rw.write_h5(filesave,data)


def save_W2(data,date,tmin,tmax,index):
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



