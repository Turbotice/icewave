


#keys_a = ['ta','ax','ay','az']
#keys_g = ['tg','gx','gy','gz']
#keys_m = ['tm','mx','my','mz']

import numpy as np
import rw_pyphone as rw
import glob
import timesync

def main():
    testfolder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_4_U400cms/'
    datafolder =  '/Volumes/labshared2/Telephones/Soufflerie_dec23/Results/'
    basefolder = 'Telephones/Soufflerie_dec23/131223/Telephones/'
    #    basefolder = rw.find_path(testfolder)

    allfolders = glob.glob(rw.find_path(basefolder)+'*/')
    print(allfolders)
    savefolder = rw.find_path('Telephones/Soufflerie_dec23/Results/')

    folder = 'Telephones/Soufflerie_dec23/131223/Telephones/121223_11h39_T1/'
    #time_sync(rw.find_path(folder),savefolder)

    folders = glob.glob(rw.find_path(basefolder)+'*_T*/')
    
    
    for folder in folders:#allfolders[11:]:
        name = folder.split('_')[-1][:-1]
        print(name,folder)
        time_sync(folder,savefolder,name)

        #print("Unzip files ...")
        #rw.extract_all(folder)
        #print("Done")
#        wind = folder.split('_')[-1].split('/')[0]
#        print(wind)
#        run(folder,savefolder,name=wind)

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
    main()
