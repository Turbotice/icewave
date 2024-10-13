


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
    parser.add_argument('-step', dest='step', type=int,default=3,help='select Step to be performed')

#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args

def step4(folders):
    #step4: load .csv file
    #       draw a map of the phone location
    pass

def step3(folder):
    import pickle
    #step3: load .pkl data
    #       compute for each phone
    #           starting and end time (using cellphone values, no sync at this step)
    #           duration, mean location, total distance traveled
    #           ai_rms, vi_rms, mi_rms
    #           main frequency f0
    #           amplitude in a frequency band [f0/2,2*f0]
    #       agregate the results for each folder, save a summary csv file with the quantity of interest
    phonefiles = glob.glob(folder+'000*/*.pkl')
    pprint(phonefiles)

    results={}
    for filename in phonefiles:
        with open(filename, 'rb') as handle:
            data = pickle.load(handle)
        phone = int(filename.split('/')[-2].split('_')[1])
        print(data.keys(),phone)
        result = averages(data)
        results[phone] = result

    rw.write_csv(results,folder,title='averages')

    #save results in .csv format
    
def step2(folder):
    #step 2 :   load data in .csv format
    #           make a dictionnary data
    #           find the measurement interval,
    #               using a criteria on the modulus of the acceleration,
    #               on time interval Dt (default 5s)
    #           cut the temporal serie
    #           save the dictionnary data in a .pkl format (one for each phone)    
        phonefolders = glob.glob(folder+'000*/')
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

def step1(folders):
    #step 1 : find all zip files, and extract the data
    for folder in folders:#allfolders[11:]:
#        name = folder.split('_')[-1][:-1]
#        print(name,folder)
#        time_sync(folder,savefolder,name)

        print("Unzip files ...")
        rw.extract_all(folder)
        print("Done")
#        wind = folder.split('_')[-1].split('/')[0]
#        print(wind)
#        run(folder,savefolder,name=wind)

def find_measure_interval(data,var='a',Dt=5,S0=1,display=False):
    #work for same sampling frequency for all sensors
    print(data.keys())
    t = data['t'+var]
    y = np.sqrt(np.sum([data[var+coord]**2 for coord in data['coords']],axis=0))

    n = int(np.round(np.sum(t<Dt)/2)*2) #windows of Dt seconds. may not work for fs = 400Hz (values are designed for fs = 50Hz)
    print(f"Number of points per bin :{n}")
    print(f"Number of bins : {N}")
    N = int(np.floor(len(t)/n))

    Y = np.reshape(y[:N*n],(N,n))# cut the signal in pieces of n points (Dt time slots)
    T = np.reshape(t[:N*n],(N,n))
    S = np.std(Y,axis=1) # compute the std on these slots
    Smean = np.mean(Y,axis=1)
        
    indices = np.where(S<S0)[0]
    if len(indices)<2:
        indices = np.where(S<S0*10)[0]
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
    variables = ['a','g','m']
    funlist = [np.mean,np.std]
    results={}
    for var in variables:
        for coord in data['coords']:
            key = var+coord
            y = data[key]
            #print(key,np.mean(y),np.std(y))
            y_high,y_wave,y_trend,err = filtering(y)
            results[key+'_err'] = err
            for fun in funlist:
                if var=='g' and fun==np.mean:
                    continue
                results[key+'_w_'+fun.__name__] = fun(y_wave)
                results[key+'_high_'+fun.__name__] = fun(y_high)
                results[key+'_trend_'+fun.__name__] = fun(y_trend)

    var = 'a'
    coord = 'z'
    key = var+coord
    t = data['t'+var]
    y = data[key]
    y_high,y_wave,y_trend,err = filtering(y)
    f,TFmoy,fmax = time_spectrum(t,y_wave)
    results[key+'_w_freq']=fmax
    
    #location
    if 'loc' in data.keys():
        for key in ['elev','lat','lon']:
            y = data['loc'][key]
            for fun in funlist:
                    results[key+'_'+fun.__name__] = fun(y)
    else:
        for key in ['elev','lat','lon']:
            for fun in funlist:
                    results[key+'_'+fun.__name__] = np.nan
    #time
    if 'time' in data.keys():
        results['time_start']=data['time']['system_START']
        results['time_end']=data['time']['system_PAUSE']
    else:
        results['time_start']='date 00:00:00.000'#data['time']['system_START']
        results['time_end']='date 00:00:00.000'#data['time']['system_PAUSE']
    return results

import scipy.signal as sig

def filtering(y,fc=0.1,flow=0.002):
    #correspond to 5Hz and 0.1Hz at fs = 50Hz
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

def time_spectrum(t,y):
#    t = np.asarray(data['t'+var])
#    y = np.asarray(data[var+coord])
    y = y-np.mean(y)

    n = y.shape[0]
    nt = 50*60
    N = int(np.floor(n/nt))
#    nt = int(np.floor(n/N))
    print("Number of samples : "+str(N))

    Y = np.reshape(y[:nt*N],(N,nt))
    win = np.transpose(np.repeat(np.reshape(np.hanning(nt),(nt,1)),N,axis=1))
    Y = Y*win

    Ypad = np.zeros((N,nt*9))
    Ypad[:,4*nt:5*nt]=Y
    Nt = 9*nt

    dtmean = np.mean(np.diff(t))
    fe = 1/dtmean
    f = np.linspace(0,fe/2,int(Nt/2))
    TF = np.abs(np.fft.fft(Ypad,axis=1))
    df = f[1]-f[0]

    TF = TF[:,:int(Nt/2)]/np.sqrt(df)/nt  #normalisation de la transformée de Fourier
    TFmoy = np.mean(TF,axis=0)#/np.sqrt(N)

    #remove first 10 points to find the maximu
    Amax = np.max(TFmoy[10:])
    i = np.argmax(TFmoy[10:])+10
    fmax = f[i]
    return f,TFmoy,fmax

def main(args):
    date = args.date
    base = '/media/turbots/Hublot24/Share_hublot/Data/'
    #date = '0221'
    datafolder = base+date+'/Telephones/'
    folders1 = glob.glob(datafolder+'Bic24*/')
    folders2 = glob.glob(datafolder+'Sag24*/')
    folders = folders1+folders2
    for folder in folders:
        #func = locals()['step'+str(args.step)]
        #func(folder)
        step2(folder)
if __name__ =='__main__':
    args = gen_parser()
    main(args)
