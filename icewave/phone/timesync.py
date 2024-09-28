
import scipy.signal as sig
import numpy as np


def timesync(data,ref=25,prominence=2):
    data[ref],tmax = process_ref(data[ref])
    yiref = data[ref]['yi']
    for key in data.keys():
        data[key] = get_dt(data[key],yiref,tmax,prominence=prominence)
        #print(data[key]['dt'])
    return data

def get_comb(d,tmin,tmax,distance=50,dt=1/2000,prominence=2):
    t = np.asarray(d['ta'])
    y = np.asarray(np.sqrt(np.asarray(d['az'])**2+np.asarray(d['ay'])**2+np.asarray(d['ax'])**2))    
    peaks = sig.find_peaks(y,distance=distance,prominence=prominence)[0]
        
    #d['tpks'] = np.asarray(d['ta'])[peaks]
    #d['vpks'] = np.asarray(d['az'])[peaks]
    yf = y*0
    yf[peaks]=y[peaks]/y[peaks]#keep the magnitude of the peaks.
    #Might not work perfectly for phones that are noiser

    ti = np.arange(tmin,tmax,dt)
    yfi = np.interp(ti,t,yf)#interpolate to give each peak some width
    return ti,yfi

def process_ref(d):
    d['dt']=0
    tmax = np.max(d['ta'])
    tiref,yiref = get_comb(d,0,tmax)

    d['ti'] = tiref
    d['yi'] = yiref
    return d,tmax

def get_dt(d,yiref,tmax,prominence=2):
    ti,yi = get_comb(d,0,tmax,prominence=prominence)
    C = sig.correlate(yi,yiref,mode="same")
    dt = ti-np.mean(ti)
    i = np.argmax(C)
    good = np.nanmax(C)/len(yi)/np.nanstd(yi)/np.nanstd(yiref)
    print(good)
    Dt = dt[i]
    d['dt'] = Dt
    d['ti'] = ti
    d['yi']= yi
    return d

def get_timetable(data,ref=25):
    res = {}
    for phone in data.keys():
        res[phone]={}
        res[phone]['dt'] = data[phone]['dt']
        res[phone]['t0'] = data[phone]['syst'][0]
        res[phone]['dt_tot'] = data[phone]['syst'][0]-data[ref]['syst'][0]+data[phone]['dt']
    return res
    
