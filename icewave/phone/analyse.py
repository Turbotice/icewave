


#keys_a = ['ta','ax','ay','az']
#keys_g = ['tg','gx','gy','gz']
#keys_m = ['tm','mx','my','mz']
import rw_pyphone as rw
import timesync
import load as load
import argparse

import numpy as np
import glob
import pylab as plt
from pprint import pprint

global BicWin2024_datas
#BicWin2024_datas={'0226':{'folder':'Bic24*/','subfolders':'0001*/','sync':}//
#                  '0223':{'folder':'Bic24*/','subfolders':'0001*/'}//
#                  '0221':{'folder':'Bic24*/','subfolders':'0001*/'}}



def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate smartphone data")
    parser.add_argument('-date', dest='date', type=str,default='0226',help='select date to process data')
#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args

def main(date):
    base = '/media/turbots/Hublot24/Share_hublot/Data/'
    #date = '0221'
    datafolder = base+date+'/Telephones/'    
    folder = glob.glob(datafolder+'Bic24*/')[0]
    print(folder)
    phonefolders = glob.glob(folder+'0001*/')
    pprint(phonefolders)

    #phonefolders=[phonefolder[:-1] for phonefolder in phonefolders]
    
    savefolder = folder+'Results/'

    for phonefolder in phonefolders:
        data = load.load(phonefolder)
        data = load.sort(data)
        data = find_measure_interval(data)
        data = cut(data)
        rw.save_data_single_phone(data,phonefolder)
#    testfolder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_4_U400cms/'
#    datafolder =  '/Volumes/labshared2/Telephones/Soufflerie_dec23/Results/'
#    basefolder = 'Telephones/Soufflerie_dec23/131223/Telephones/'
    #    basefolder = rw.find_path(testfolder)

#    allfolders = glob.glob(rw.find_path(basefolder)+'*/')
#    print(allfolders)
#    savefolder = rw.find_path('Telephones/Soufflerie_dec23/Results/')

#    folder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_11h39_T1/'
    #time_sync(rw.find_path(folder),savefolder)
#   folders = glob.glob(rw.find_path(basefolder)+'*_T*/')
    
#    for folder in folders:#allfolders[11:]:
#        name = folder.split('_')[-1][:-1]
#        print(name,folder)
#        time_sync(folder,savefolder,name)

        #print("Unzip files ...")
        #rw.extract_all(folder)
        #print("Done")
#        wind = folder.split('_')[-1].split('/')[0]
#        print(wind)
#        run(folder,savefolder,name=wind)

def find_measure_interval(data,var='a',Dt=5,S0=1,display=False):
    #work for same sampling frequency for all sensors
    print(data.keys())
    t = data['t'+var]
    y = np.sqrt(np.sum([data[var+coord]**2 for coord in data['coords']],axis=0))

    n = int(np.round(np.sum(t<Dt)/2)*2) #windows of Dt seconds.
    print(f"Number of points per bin :{n}")
    N = int(np.floor(len(t)/n))

    Y = np.reshape(y[:N*n],(N,n))# cut the signal in pieces of n points (Dt time slots)
    T = np.reshape(t[:N*n],(N,n))
    S = np.std(Y,axis=1) # compute the std on these slots
    Smean = np.mean(Y,axis=1)
        
    indices = np.where(S<S0)[0]
    data[var+'i0']=indices[1]*n
    data[var+'i1']=indices[-2]*n
    data[var+'t0']=t[indices[1]*n]
    data[var+'t1']=t[indices[-2]*n]

    if display:
        fig,axs = plt.subplots(figsize=(10,6),nrows=2,sharex=True)#,nrows=3,sharex=True)
        axs[0].plot(t,y,'k')
        axs[0].errorbar(T[:,int(n/2)],Smean,S,marker='o',color='r')
        axs[0].set_ylim([5,20])
        axs[0].vlines(t[indices[1]*n],0,20)
        axs[0].vlines(t[indices[-2]*n],0,20)

        text = t[data[var+'i0']:data[var+'i1']]
        yext = y[data[var+'i0']:data[var+'i1']]
        axs[1].plot(text,yext,'k')
        axs[1].set_ylim([8,12])
        fig.subplots_adjust(hspace=0)

    return data

def cut(data,v='a'):
    variables = ['a','g','m']
    t0 = data[v+'t0']
    t1 = data[v+'t1']
    for var in variables:
        t = data['t'+var]
        indices=np.logical_and(t>=t0,t<=t1)
        data['t'+var]=data['t'+var][indices]
        for coord in data['coords']:
            data[var+coord]=np.asarray(data[var+coord])[indices]

    # same for GPS data
    if 'loc' in data.keys():
        t = data['loc']['t']
        indices = np.logical_and(t>=t0,t<=t1)
        for key in ['lat','elev','lon','t']:
            data['loc'][key]=data['loc'][key][indices]
    Dt = np.round(t1-t0,decimals=1)
    print(f"Duration of recording : {Dt}")
    return data

def spatiotemp(basefolder,savefolder,datafolder):
    folderlist = rw.get_phone_list(basefolder)
    data={}

    tmin = 60
    tmax = 300
    ti = np.arange(tmin,tmax,0.02)
    n = len(ti)
    nphone = len(folderlist)
    M = np.zeros((n,nphone))

    timefile = glob.glob(datafolder+'timetable*.csv')[0]
    times = load_timetable(timefile)

    ref = 1
    for folder in folderlist:
        phone = rw.get_phone_num(folder)
        print(phone)
        d = rw.load_data(folder)
        ti,yi = extract_angle(d,times[phone],tmin,tmax)        
        M[:,phone-ref] = yi

    
def extract_angle(d,dt,tmin,tmax,tstep=0.02,fcut=0.001,order=4):
    t = np.asarray(d['tg'])
    y = np.asarray(d['gx'])
    
    y = y-np.mean(y)
    [b,a] = sig.butter(order,Wn=[fcut],btype='high')
    yfilt = sig.filtfilt(b, a, y)

    n = len(y)
    y= np.asarray(yfilt)
    
    t = t-dt#times[key]
    
    indices = np.logical_and(t>tmin-1,t<tmax+1)
    
    y = np.cumsum(y)
    y = y-np.mean(y)
    
    t = t[indices]
    y = y[indices]

    ti = np.arange(tmin,tmax,tstep)
    yi = np.interp(ti,t,y)
    return ti,yi

def load_timetable(timefile):
    import csv
    times={}
    with open(timefile, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')    
        for i,row in enumerate(spamreader):
            #print(row)
            if i>0:
                times[int(row[0])]=float(row[3])#old version : row[1]
#        print(', '.join(row))
    return times#pprint(times)


#for oscillating phones : az, gx, mz
def time_sync(basefolder,savefolder,name,ref=25,write=True):
    folderlist = rw.get_phone_list(basefolder)
    print(folderlist)
    data={}
    for folder in folderlist:
        phone = rw.get_phone_num(folder)
        print(phone)
        d = rw.load_data(folder)
        data[phone]=d
    data = timesync.timesync(data,ref=ref)
    results = timesync.get_timetable(data,ref=ref)

    if write:
        rw.write_csv(results,savefolder,title='timetable_'+name+'_')
    return data,results

def run(basefolder,savefolder,name=''):
    folderlist = rw.get_phone_list(basefolder)
    #print(folderlist)
    results={}
    for folder in folderlist:
        phone = rw.get_phone_num(folder)
        #print(phone)
        data = rw.load_data(folder)
        results[phone]= analyse(data)
        
    rw.write_csv(results,savefolder,title=name+'_v2_')
    

def analyse(data):
    results={}

    #compute averages
    #res = averages(data)
    #results.update(res)

    #compute damping coefficients
    res = damping(data)
    results.update(res)
    #print(results)
    return results

def averages(data,keys='all'):
    xyz = ['x','y','z']
    variables = ['a','g','m']

    results={}
    for var in variables:
        for l in xyz:
            key = var+l
            y = data[key]
            #print(key,np.mean(y),np.std(y))
            for fun in [np.mean,np.std]:
                results[key+'_'+fun.__name__] = fun(y)
    return results


import scipy.signal as sig

def damping(data):
    results={}
    var = 'g'
    l = 'x'
    key = var+l
    
    y = data[key]
    dt = np.mean(np.diff(data['tg']))
    Lmb,Lmb_err=Lambda(y,dt)

    results[key+'_Lambda']=Lmb
    results[key+'_Lambda_err']=Lmb_err

    return results
            
def Lambda(y,dt,twin=8,dist=50,fcut=0.001):
    #    twin : window length in seconds
    [b,a] = sig.butter(4,Wn=[fcut],btype='high') #high pass filter at 0.4Hz
    yfilt = sig.filtfilt(b, a, y)
    y = yfilt-np.mean(yfilt) #remove average
#    y = np.asarray(y)#[int(n/6):]
    y = np.cumsum(y)#compute position over time
    n = len(y)

    Nb = int(twin/dt)
    c = int(np.floor(n/Nb))

    Ctab = np.zeros((Nb*2-1,c-2))#do not count the two extremities for safety
    for i in range(1,c-1):
        Ctab[:,i-1] = sig.correlate(y[i*Nb:(i+1)*Nb],y[i*Nb:(i+1)*Nb])
        tC = dt*np.arange(-Nb+1,Nb,1)#*np.linspace(-Nb+1,Nb-1,Nb*2+1)
    Cmoy=np.mean(Ctab,axis=1)#compute the averaged correlation

    [b,a] = sig.butter(4,Wn=[fcut],btype='high') #filter again the low frequency variations
    Cfilt = sig.filtfilt(b, a, Cmoy)

    peaks = sig.find_peaks(np.abs(Cfilt),distance=dist)[0]
    indices = [tC[peaks]>=0][0]#keep only the positive events

#    print(indices)
    nmax = np.min([10,len(indices)])
#    print(nmax)
    x = tC[peaks][indices][:nmax]
    y = np.abs(Cfilt)[peaks][indices][:nmax]
    p = np.polyfit(x,np.log(y),1)
#    plt.plot(x,np.exp(np.polyval(p,x)),'r--')
    print(p[0])
    return p[0] #damping coefficient in s$^{-1}$

def time_spectrum(data,key):
    TFt = np.fft.rfft(data[key])


if __name__ =='__main__':
    args = gen_parser()
    main(args)
