import numpy as np
import pylab as plt
import glob
import os

import scipy
import scipy.signal as sig
import stephane.display.graphes as graphes



def get_channels(a):
    """
    Correspondance between columns and channels
    """
    d={}
    d['Z1'] = a[:,0]
    d['L1'] = a[:,1]
    d['T1'] = a[:,2]
    d['Z2'] = a[:,3]
    d['L2'] = a[:,4]
    d['T2'] = a[:,5]
    d['Z3'] = a[:,6]
    d['L3'] = a[:,7]
    d['T3'] = a[:,8]
    return d

def compute_all(d):
    Corr = {}
    
    for k in ['Z','L','T']:
        for i in ['1','2','3']:
            for j in ['1','2','3']:
                Corr[k+i+j],dt = corr(d[k+i],d[k+j])

    for k1,k2 in [('Z','L'),('L','T'),('T','Z')]:
        for i in ['1','2','3']:
            Corr[k1+k2+i+i],dt = corr(d[k1+i],d[k2+i])
            
    return d,Corr,dt

def corr(y1,y2,facq=1000):
    n = len(y1)
    C = sig.correlate(y1-np.mean(y1),y2-np.mean(y2))/(np.std(y1)*np.std(y2))/n
    dt = np.arange(-n+1,n,1)/facq
    
    return C,dt

def get_max(C,dt):
    indmin = np.argmin(C)
    indmax = np.argmax(C)

    vmin = C[indmin]
    vmax = C[indmax]
    
    return indmin,indmax,vmin,vmax

def extract(C,dt,i0=0,b=500):
    N = len(dt)
    n = N//2-1
    indices = range(n-b+i0,n+b+i0)
    Ce = C[indices]-np.mean(C[indices])
    return Ce,dt[indices]

def display_corr(C,dt,b=500,ax=None,name='',figs=None):
    Ce,dte = extract(C,dt,b=b)
    indmin,indmax,vmin,vmax = get_max(Ce,dte)
        
    if ax is None:
        plt.plot(dte,Ce)
        plt.plot(dte[indmin],vmin,'rv')
        plt.plot(dte[indmax],vmax,'r^')     
    else:
        ax.plot(dte,Ce)
        ax.plot(dte[indmin],vmin,'rv')
        ax.plot(dte[indmax],vmax,'r^') 
        
    if figs is None:
        figs={}
    figs.update(graphes.legende(r'$\Delta t$ (s)',r'$C_{'+name+'}$','',ax=ax))
    return figs

def display_corrs(Z,dt,b=500,ax=None,name='',figs=None):
    for C in Z:
        Ce,dte = extract(C,dt,b=b)
        if ax is None:
            plt.plot(dte,Ce)  
        else:
            ax.plot(dte,Ce)

    
    for C in Z:
        indmin,indmax,vmin,vmax = get_max(Ce,dte)
        
        if ax is None:
            plt.plot(dte[indmin],vmin,'rv')
            plt.plot(dte[indmax],vmax,'r^')     
        else:
            ax.plot(dte[indmin],vmin,'rv')
            ax.plot(dte[indmax],vmax,'r^') 
    if figs is None:
        figs={}
    figs.update(graphes.legende(r'$\Delta t$ (s)',r'$C_{'+name+'}$','',ax=ax))
    return figs


def display_map(Z,dt,name='',b1=5000,b2=500):
    fig,axs = plt.subplots(ncols=2,nrows=1,figsize=(10,5))
    figs = {}

    figs = display_corrs(Z,dt,b=b1,ax=axs[0],name=name,figs=figs)
    figs = display_corrs(Z,dt,b=b2,ax=axs[1],name=name,figs=figs)

    return figs


def env(X,t):
    A = X[0]
    phi = X[1]
    w = X[2]
    t0 = 16#X[3]

    n = len(t)
    modul = A*np.max([np.zeros(n),(t0-np.abs(t))/t0],axis=0)
    Cth = np.cos(w*t+phi)*modul
    return Cth

def error(X,t,C,fun=env):
    Cth = fun(X,t)
    return np.sum((C-Cth)**2)


def fit(t,C):
    x0 = [0.2,np.pi,1.5,15]
    xf = scipy.optimize.fmin(error,x0,args=(t,C))
    Cth = env(xf,t)
    return xf,Cth
