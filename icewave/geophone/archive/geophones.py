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
    Z : composante verticale 
    L : composante 'N', normalement longitudinal (si geophone en ligne)
    T : composante 'E', normalement transversale (si geophone en ligne) 
    1,2,3 : géophones, normalement numéroté G4,G6,G12. Voir cahier de manip pour détail
    """
    print(a.shape)
    nt,nc = a.shape
    
    global geokeys

    d={}
    d['keys'] = []
    if nc >=3:
        d['Z1'] = a[:,0]
        d['L1'] = a[:,1]
        d['T1'] = a[:,2]
        d['keys'].append('1')
    if nc >=6:
    	d['Z2'] = a[:,3]
    	d['L2'] = a[:,4]
    	d['T2'] = a[:,5]
    	d['keys'].append('2')
    if nc >=9:
    	d['Z3'] = a[:,6]
    	d['L3'] = a[:,7]
    	d['T3'] = a[:,8]
    	d['keys'].append('3')
    return d
    
def compute_all(d):
    Corr = {}
    
    for k in ['Z','L','T']:
        for i in d['keys']:
            for j in d['keys']:
                Corr[k+i+j],dt = corr(d[k+i],d[k+j])

    for k1,k2 in [('Z','L'),('L','T'),('T','Z')]:
        for i in d['keys']:
            Corr[k1+k2+i+i],dt = corr(d[k1+i],d[k2+i])
            
    return d,Corr,dt

def corr(y1,y2,facq=1000,detrend=False):
    n = len(y1)
    if detrend:
        y1 = sig.detrend(y1, axis=-1, type='linear', bp=0, overwrite_data=False)
        y2 = sig.detrend(y2, axis=-1, type='linear', bp=0, overwrite_data=False)
        
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

def display_corr(C,dt,b=500,ax=None,name='',title='',figs=None):
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
    figs.update(graphes.legende(r'$\Delta t$ (s)',r'$C_{'+name+'}$',title,ax=ax))
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
    """
    Theoretical form of temporal correlation, A*cos(wt+phi)*Triangle[(t0-t)/t0)]
    """
    A = X[0]
    phi = X[1]
    w = X[2]
    t0 = X[3]

    n = len(t)
    modul = A*np.max([np.zeros(n),(t0-np.abs(t))/t0],axis=0)
    Cth = np.cos(w*t+phi)*modul
    
    param = {}
    if X[0]<0:
        X[0] = -X[0]
        X[1] = X[1]+np.pi
        
    param['A'] = X[0]
    param['phi'] = X[1]
    param['f'] = X[2]/(2*np.pi)
    param['tc'] = X[3]
    
    return Cth,param

def error(X,t,C,fun=env):
    Cth,param = fun(X,t)
    return np.sum((C-Cth)**2)

def fit(t,C,x0 = [0.5,np.pi,1,15]):
    xf = scipy.optimize.fmin(error,x0,args=(t,C))
    Cth,param = env(xf,t)
    return xf,Cth,param
