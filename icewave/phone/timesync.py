
import scipy.signal as sig
import numpy as np

import pylab as plt
import icewave.display.graphes as graphes

def timesync(data,ref=25,prominence=2):
    data[ref],tmax = process_ref(data[ref],prominence=prominence)
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
    yf[peaks]=y[peaks]/y[peaks]#remove the magnitude of the peaks.
    #Might not work perfectly for phones that are noiser

    ti = np.arange(tmin,tmax,dt)
    yfi = np.interp(ti,t,yf)#interpolate to give each peak some width
    return ti,yfi

def process_ref(d,prominence=2):
    d['dt']=0
    tmax = np.max(d['ta'])
    tiref,yiref = get_comb(d,0,tmax,prominence=prominence)
    
    d['ti'] = tiref
    d['yi'] = yiref
    return d,tmax

def get_dt(d,yiref,tmax,prominence=2):
    ti,yi = get_comb(d,0,tmax,prominence=prominence)

    #print(np.sum(np.isnan(yi)))
    #print(np.sum(np.isnan(yiref)))

    C = sig.correlate(yi,yiref,mode="same")
    dt = ti-np.mean(ti)
    i = np.argmax(C)
        
    good = np.max(C)/len(yi)/np.nanstd(yi)/np.nanstd(yiref)
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

def compute_table_sync(data_sync,pair=(0,1)):
    phonelist = np.sort(list(data_sync[0].keys()))#[0,4,6,9,11,16,17,18]#np.sort(list(res[1].keys()))#[0,4,6,9,11,13,16,17,18,19]

    n = len(phonelist)

    tab_dt=np.zeros((n,n))
    tab_dt_tot={}
    for i in data_sync.keys():
        tab_dt_tot[i] = np.zeros((n,n))

    res = {}
    for j1,ref in enumerate(phonelist):
        for i in data_sync.keys():
            data_sync[i] = timesync(data_sync[i],ref=ref,prominence=0.1)
            res[i] = get_timetable(data_sync[i],ref=ref)

        print('Compare time scales for ref'+str(ref))
        for j2,key in enumerate(phonelist):
            dt = res[pair[1]][key]['dt_tot']-res[pair[0]][key]['dt_tot']
            tab_dt[j1,j2]=dt
            for i in data_sync.keys():
                tab_dt_tot[i][j1,j2]=res[i][key]['dt_tot']
    return phonelist,tab_dt_tot,tab_dt
            
def check_timesync_accuracy(data_sync,ref,prominence=0.1,pair=(0,1)):
    res = {}
    for i in data_sync.keys():
        data_sync[i] = timesync(data_sync[i],ref=ref,prominence=prominence)
        res[i] = get_timetable(data_sync[i],ref=ref)

    print('Compare time scales')
    notsyncs = []

    dt_final={}
    for i in data_sync.keys():
        dt_final[i]={}

    if len(list(data_sync.keys()))<2:
        print("Need at least two time sync recordings to check timesync accuracy")
    else:
        for key in res[0].keys():
            dt = res[pair[1]][key]['dt_tot']-res[pair[0]][key]['dt_tot']
            #print(key,dt)
            if np.abs(dt)<0.1:
                print(key,str(np.round(dt*1000,decimals=1))+' ms')
                for i in res.keys():
                    dt_final[i][key]={}
                    dt_final[i][key]['dt']=res[i][key]['dt']
                    dt_final[i][key]['t0']=res[i][key]['t0']
                    dt_final[i][key]['dt_tot']=res[i][key]['dt_tot']
                    dt_final[i][key]['ref']=ref
            else:
                for i in res.keys():
                    dt_final[i][key]={}
                    dt_final[i][key]['dt']=np.nan
                    dt_final[i][key]['t0']=res[i][key]['t0']
                    dt_final[i][key]['dt_tot']=np.nan#res[i][key]['dt_tot']
                    dt_final[i][key]['ref']=np.nan
                    notsyncs.append(key)
                    print(str(key)+' not sync')
    return dt_final,notsyncs


def correct_timesync(phonelist,notsyncs,dt_final,tab_dt,tab_dt_tot):
    phonelist = list(phonelist)
    phonesyncs=np.ones(len(phonelist),dtype=bool)
    for j,phone in enumerate(phonelist):
        if not phone in notsyncs:
            phonesyncs[j]=True
        else:
            phonesyncs[j]=False

    for phone in notsyncs:
        j1 = np.where(np.asarray(phonelist)==phone)[0]
        #print(tab_dt[:,j2])
        #find the ref phone :
        #smallest error among the telephones that are already synced
        j2 = np.argmin(np.abs(tab_dt[j1,phonesyncs]))
        #get the corresponding error
        dt = tab_dt[j1,phonesyncs][j2]
        phoneref = np.asarray(phonelist)[phonesyncs][j2]
        print(phone,phoneref,dt)
        
        key = phone
        for i in dt_final.keys():
#            dt_final[i][key]={}
            dt_final[i][key]['dt']=dt
            dt_final[i][key]['dt_tot']=tab_dt_tot[i][j1,phonesyncs][j2]
            dt_final[i][key]['ref']=phoneref
    return dt_final

def merge_timesync_references(dt_final,ref):
    for i in dt_final.keys():
        for key in dt_final[i].keys():
            if dt_final[i][key]['ref']!=ref:
                print(dt_final[i][key]['ref'])
                key_int = dt_final[i][key]['ref']
                dt_final[i][key]['dt'] = dt_final[i][key]['dt']+dt_final[i][key_int]['dt']
                dt_final[i][key]['dt_tot'] = dt_final[i][key]['dt_tot']+dt_final[i][key_int]['dt_tot']
                dt_final[i][key]['ref'] = ref
    return dt_final


def display_timesync_table(tab_dt,phonelist,vmax=200):
    n = len(phonelist)
    cm = plt.colormaps['Oranges']

    fig,ax = plt.subplots(figsize=(5,4))#,sharex=True)

    sc = ax.pcolormesh(np.abs(tab_dt)*1000,vmin=0,vmax=vmax,cmap=cm)
    plt.xticks(np.arange(n)+0.5,labels=phonelist)
    plt.yticks(np.arange(n)+0.5,labels=phonelist)
    #plt.colorbar()
    plt.axis('equal')
    #cbar_ax = fig.add_axes([0.87, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(sc)#, cax=cbar_ax,ticks=[1,3,10,30])

    #cbar = plt.colorbar(sc)#,
    #cbar.ax.set_yticklabels(['1','3','10','30'])#,'$10^1$'])  
    cbar.set_label(r'$\Delta t$ (ms)', rotation=270,fontsize=18,labelpad=20)

    figs = graphes.legende(r'Reference phone #',r'Phone #','',cplot=True)
    #graphes.save_figs(figs,savedir=savefolder,prefix=f'Time_table_{date}_',overwrite=True)
    return figs

def exemple():
    
    phonelist,tab_dt_tot,tab_dt = compute_table_sync(data_sync,pair=(0,1))
    display_timesync_table(tab_dt,phonelist,vmax=200)
    ref = 1
    dt_final,notsyncs = check_timesync_accuracy(data_sync,ref,prominence=0.1,pair=(0,1))
    dt_final = correct_timesync(phonelist,notsyncs,dt_final,tab_dt,tab_dt_tot)
    dt_final = merge_timesync_references(dt_final,ref)



    
