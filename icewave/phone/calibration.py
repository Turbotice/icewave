import numpy as np
import pylab as plt
import glob
import os

import scipy.signal as sig

import stephane.display.graphes as graphes
from pprint import pprint

import icewave.phone.load as load
import icewave.phone.analyse as analyse

import icewave.tools.rw_data as rw


def load_param_table(filename):
    param = {}

    raw = rw.read_csv(filename)
    data = rw.csv2dict(raw)
    for k in data.keys():
        key = k.split(' ')[1]
        param[key] = data[k]

    return param


def load_data(param,datafolder,key='rate'):
    datafiles = glob.glob(datafolder+'*.csv')
    datas={}
    for i,fichier in enumerate(param['fichier']):
        files=[]
        try:
            k = int(param[key][i])#.keys()
        except:
            k = np.round(float(param[key][i]),decimals=1)#.keys()

        for datafile in datafiles:
            file = os.path.basename(datafile).split('_D')[1].split('-android')[0]
            if fichier==file:
                #print(file,os.path.basename(datafile))
                files.append(datafile)

        if len(files)>0:
            datas[k]={}
            for filename in files:
                datas[k].update(load.load_gobfile(filename))

    return datas

def display_data(datas,var='g'):
    n = len(datas.keys())
    fig,axs = plt.subplots(figsize=(15,2*n),nrows=n,sharex=True)

    for key,ax in zip(datas.keys(),axs):
        data = datas[key]
        Vars = [f'${var}_{c}$ (rad/s)' for c in data['coords']]

        t = data['t'+var]
        t = (t-t[0])/1000000 #convert in seconds
        for c in data['coords']:
            y = data[var+c]
    #        if var+c =='a'
            y = y# -np.mean(y)
            nb = int(len(y)/2)
            ax.plot(t[nb:],y[nb:])
        ax.legend(Vars)
        figs = graphes.legende('',str(key),'',ax=ax)

    figs = graphes.legende('$t$ (s)',str(key),'',ax=axs[-1])
    fig.subplots_adjust(hspace=0)

    return figs,axs

def compute_ar(datas,c='x'):
    n = len(datas.keys())
    var = 'a'
    for key in datas.keys():
        data = datas[key]
        rate = int(key)
        y = data[var+c]
        nb = int(len(y)/2)#cut the first half to avoid transient regime
        ax = np.mean(y[nb:])
        datas[key]['Acc_r']=ax
    return datas

def display_ar_a(datas,axis='x'):
    fig,ax = plt.subplots(figsize=(4.78,4.78))#,nrows=n,sharex=True)

    Xs = [int(key) for key in datas.keys()]
    ar = [datas[key]['Acc_r'] for key in datas.keys()]
    
    p = np.polyfit(Xs,ar,1)
    Xth = np.linspace(np.min(Xs),np.max(Xs),100)

    X0 = -p[1]/p[0]
    ax.plot(Xs,ar,'ko')
    ax.plot(Xth,np.polyval(p,Xth),'r--')
    ax.plot(X0,0,'r*',markersize=14)

    title = f'${axis}_S =$ '+str(np.round(X0,decimals=2))+' mm'
    #ax.set_xlim([0,110])
    #ax.set_ylim([0,10])
#    plt.legend(['$g_y$','$c$ ='+str(np.round(p[0],decimals=4))])
    figs = graphes.legende(f'${axis}$ (mm)','$a_r$ (m/s$^{2}$)',title,ax=ax)
    #fig.subplots_adjust(hspace=0)

    return figs,ax


def compute_Omega_g(datas,c='y'):
    n = len(datas.keys())
    var = 'g'
    for key in datas.keys():
        data = datas[key]
        rate = int(key)
        y = data[var+c]
        nb = int(len(y)/2)#cut the first half to avoid transient regime
        gy = np.mean(y[nb:])
        datas[key]['Omega_g']=np.abs(gy)        
    return datas

def display_Omega_g(datas):
    fig,ax = plt.subplots(figsize=(4.78,4.78))#,nrows=n,sharex=True)

    Rates = [int(key) for key in datas.keys()]
    w_g = [datas[key]['Omega_g'] for key in datas.keys()]
    
    p = np.polyfit(Rates,w_g,1)
    Oth = np.linspace(0,np.max(Rates),100)

    ax.plot(Oth,Oth*p[0],'r--')
    ax.plot(Rates,w_g,'ko')

    #ax.set_xlim([0,110])
    #ax.set_ylim([0,10])
    plt.legend(['$g_y$','$c$ ='+str(np.round(p[0],decimals=4))])
    figs = graphes.legende('Command value (a.u.)','$\Omega$ (rad/s)','',ax=ax)
    #fig.subplots_adjust(hspace=0)

    return figs,ax

def compute_Omega_m(datas,c='x'):
    for key in datas.keys():
        data = datas[key]
        tm = data['tm']
        y = data['m'+c]

        nt = len(y)

        f,TFmoy,fmax,Amax = analyse.time_spectrum(tm/10**6,y,nt=nt)
        datas[key]['Omega_m']=fmax*2*np.pi
    return datas

def compare_Omegas(datas):
    Omega_g = [datas[key]['Omega_g'] for key in datas.keys()]
    Omega_m = [datas[key]['Omega_m'] for key in datas.keys()]

    fig,ax = plt.subplots(figsize=(4.78,4))#,nrows=n,sharex=True)
    ax.plot(Omega_m,Omega_g,'ko')
    Oth = np.linspace(0,np.max(Omega_m),100)
    ax.plot(Oth,Oth,'r--')

    figs = graphes.legende('$\Omega_m$ (rad/s)','$\Omega_g$ (rad/s)','',ax=ax)
    return figs,ax

def display_all(data):
    Vars = ['a','g','m']
        
    for var,ax in zip(Vars,axs):
        #print(data.keys())
        t = data['t'+var]
        t = (t-t[0])/1000000
        for c in data['coords']:
            y = data[var+c]
    #        if var+c =='a'
            y = y# -np.mean(y)
            nb = int(len(y)/2)
            ax.plot(t[nb:],y[nb:])
            

def exemple():
    base = '/Volumes/Fabien_2024/Telephones/Gyro_FP_Rotation/'
    savefolder = base+'Results/'
    datafolder = base+'Data/'

    tablefiles = glob.glob(base+'*.txt')

    pprint(list(tablefiles))

    headers ={}
    for filename in tablefiles:
        key = os.path.basename(filename).split('.')[0]
        headers[key] = {}
        print(key)
        headers[key] = calib.load_param_table(filename)

    name = 'axe_X_omega'
    param = headers[name]

    datas = calib.load_data(param,datafolder)

    fig,axs = calib.display_data(datas)    

    datas = calib.compute_Omega_g(datas,c='y')
    figs,ax = calib.display_Omega_g(datas)

    datas = calib.compute_Omega_m(datas)

    figs,ax = calib.compare_Omegas(datas)
    graphes.save_figs(figs,savedir=savefolder,prefix='calibration_gyro_',overwrite=True)

    pass
