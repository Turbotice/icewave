import matplotlib.pyplot as plt
import glob
import numpy as np
import os
import scipy.signal as sig
import pickle

import stephane.display.graphes as graphes
import icewave.geophone.geophones as geo
import icewave.tools.browse as browse


import socket
import platform

#Global variables
global osname,ostype
ostype = platform.platform().split('-')[0]
osname = socket.gethostname()

print(ostype,osname)

def demod(y,f0,facq=26.5):
    n = len(y)
    dt = 1/facq
    t = np.arange(0,n*dt,dt)
    
    return np.mean(y*np.exp(1j * 2 * np.pi * t * f0),axis=0,dtype='complex64')

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def compute_simple():
    frange = np.arange(0,10,0.1)

def compute_dephasage(Xfilt,Zfilt,facq=26.5):
    fig,axs = plt.subplots(nrows=1,ncols=3,figsize=(12,6))
    #fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(12,6))
    N = len(Xfilt)
    print(N)
    Np = int(2**np.floor(np.log2(N)))
    n=int(Np/2**5)
    nb=int(Np/n*2-1)

    print(n)
    dt = 1/facq
    
    frange = np.arange(0,10,0.05)
    C = np.zeros((n,len(frange)),dtype='complex64')
    Cx = np.zeros((n,len(frange)))#,dtype='complex64')
    Cz = np.zeros((n,len(frange)))#,dtype='complex64')
    
    for j,f0 in enumerate(frange):
        for i in range(n):
            X = Xfilt[i*int(Np/n):(i+1)*int(Np/n)]
            cX = demod(X,f0)
            
        
            Z = Zfilt[i*int(Np/n):(i+1)*int(Np/n)]
            cZ = demod(Z,f0)
            #print(cX)
            C[i,j]=cX*np.conj(cZ)#sig.correlate(X,Z)/(sigX*sigZ)
            Cx[i,j] = np.abs(cX)
            Cz[i,j] = np.abs(cZ)
            
    [R,Theta] = cart2pol(np.real(C),np.imag(C))

    Theta = np.unwrap(Theta,axis=1)
    axs[0].pcolormesh(frange,range(n),np.log10(R))
    axs[1].pcolormesh(frange,range(n),np.log10(Cx))
    axs[2].pcolormesh(frange,range(n),np.log10(Cz))

    #plt.colorbar()
    #for j,f0 in enumerate(frange):
    #    ax.plot(np.real(C[:,j]),np.imag(C[:,j]),'ko')


    A_moy = np.mean(R,axis=0)
    Phi_moy = np.unwrap(np.mean(Theta,axis=0))
    
    ind = np.argmax(A_moy)

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.plot(frange,A_moy,'ko')

    plt.show()

    #axs[0].plot(R[:,ind],'ko')
    #graphes.legende('#',r'|$c_x$ $c_z$*|','',ax=axs[0])
    #axs[1].plot(Theta[:,ind],'ko')
    #graphes.legende('#',r'arg($c_x$ $c_z$*)','',ax=axs[1])
    
    #ax.plot(frange,Phi_moy,'+')
    return A_moy
    #return T,Cmoy,param,figs,ax
    
def compute_corr(Xfilt,Zfilt,facq=26.5,f0=1,A0=1,display=True):
    if display:
        fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))

    N = len(Xfilt)
    Np = int(2**np.floor(np.log2(N)))
    n=int(Np/2**5)
    nb=int(Np/n*2-1)
    C = np.zeros((n,nb))

    dt = 1/facq
    for i in range(n):
        X = Xfilt[i*int(Np/n):(i+1)*int(Np/n)]
        X = X-np.mean(X)
        sigX = np.sqrt(np.sum(X**2))
    
        Z = Zfilt[i*int(Np/n):(i+1)*int(Np/n)]
        Z = Z-np.mean(Z)
        sigZ = np.sqrt(np.sum(Z**2))
    
        C[i,:]=sig.correlate(X,Z)/(sigX*sigZ)

        
        T = np.linspace(-nb/2*dt,nb/2*dt,nb)

        if display:
            ax.plot(T,C[i,:])
    
    Cmoy = np.median(C,axis=0)

    ind = np.argmax(Cmoy)
    vmax = Cmoy[ind]
    phi0 = T[ind]*f0*2*np.pi

    xlist = [vmax,phi0,f0*2*np.pi,6/f0]
    print(xlist)  #A, phi, f, tau
    
    xf,Cth,param = geo.fit(T,Cmoy,x0=xlist) #
    print(param)

    if display:
        figs,ax = display_corr(T,Cmoy,Cth,C)
    else:
        figs,ax=None,None
        
    err = error(Cmoy,Cth)
    return T,Cmoy,Cth,C,err,param,figs,ax

def display_corr(T,Cmoy,Cth,C):
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    [n,nb] = C.shape
    for i in range(n):
        ax.plot(T,C[i,:])
        
    ax.plot(T,Cmoy,'k--')
    ax.plot(T,Cth,'ro--')
    
#    ax.plot(T[ind],vmax,'bs')
    figs={}
    figs.update(graphes.legende(r'$\tau$ (frames)','C',r'$\langle X(t) Z(t+\tau) \rangle$'))
    return figs,ax

def error(Y1,Y2):
    var = np.sum((Y1-Y2)**2)/np.sqrt(np.sum(Y1**2)*np.sum(Y2**2))
    return err

    
def filtered(data,fc=1.5,facq=26.5):
    if len(fc)==1:
        [b,a] = sig.butter(4,fc[0]/facq,btype='high')
    if len(fc)==2:
        [b,a] = sig.butter(4,np.asarray(fc)/facq,btype='bandpass')
    
    X = sig.filtfilt(b,a,data['x'])
    Z = sig.filtfilt(b,a,data['z'])
    return X,Z

def load(filename,facq=26.5):
    with open(filename, 'rb') as f:
        data = pickle.load(f)

    print(len(data['x']),len(data['z']))
    n = len(data['x'])
    dt = 1/facq
    t = np.linspace(0,n*dt,n)
    data['t'] = t
    data['facq'] = facq
    
    return data

def display_signal(data,facq=26.5,title=''):
    fig,axs = plt.subplots(nrows=3,ncols=1,figsize=(12,8))
    axs[0].plot(data['t'],data['x'])
    figs = graphes.legende(r'$t$ (s)',r'$X$ (mm)',title,ax=axs[0])
    axs[1].plot(data['t'],data['z'])
    figs = graphes.legende(r'$t$ (s)',r'$Z$ (mm)','',ax=axs[1])
    axs[2].plot(data['x'],data['z'])
    graphes.legende(r'$X$ (mm)',r'$Z$ (mm)','',ax=axs[2])
    
    return figs,axs
    
#plt.plot(Zfilt)
def error(Y,Yth):
    return np.sum((Y-Yth)**2)/np.sum(Y**2)
    

def process(filename,savefolder,title=''):
    basename = os.path.basename(filename).split('.')[0]

    
    data = load(filename)

    X,Z = filtered(data,fc=[0.05],facq=26.5)#1.5#high-pass filter
    data['x'] = X
    data['z'] = Z
    
    figs,axs = display_signal(data,title=title)
    #Fmin = np.linspace(0.2,5,10)
    #Fmax = np.linspace(0.5,6,10)
    #res = {}

    
    #for (fmin,fmax) in zip(Fmin,Fmax):
    #    X,Z = filtered(data,fc=[fmin,fmax],facq=26.5)#1.5
    #figs,axs = display_signal(data)
    #compute_dephasage(X,Z)
    return figs,axs

def process_fuck(filename,savefolder):
    basename = os.path.basename(filename).split('.')[0]
    data = load(filename)
    figs,axs = display_signal(data)

    Fmin = np.linspace(0.2,5,10)
    Fmax = np.linspace(0.5,6,10)
    res = {}

    for (fmin,fmax) in zip(Fmin,Fmax):
        fmoy = fmin#(fmin+fmax)/2
        
        T,Cmoy,Cth,C,err,param,figs,ax = compute_corr(X,Z,facq=26.5,f0=fmoy,A0=0.5,display=False)
        print(err)
        if err<0.1:
            res[(fmin,fmax)] = {}
            for key in param.keys():
                res[(fmin,fmax)][key] = param [key]        
        #plt.show()
            #data = res[(fmin,fmax)]
            #fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8))
            #ax.plot(data['f'],data['phi'],'ko')
            print(param['phi'])
            #display_corr(T,Cmoy,Cth,C)
            #plt.show()
    #T,C,xf,figs,ax = compute_corr(X,Z)
    #plt.show()
    #    graphes.save_figs(figs,savedir=savefolder,overwrite=True,suffix=basename)
    #param = {}
    return res
    #graphes.save_figs(figs,savedir=savefolder)


def main2():
    base = 'Banquise/Sebastien/Experiences/'
    #exemple
    base = browse.find_path(base,disk='labshared2')
    folders = base+'Test_frequency*/Datas/'

    filelist = glob.glob(folders+'data*.pkl')
    print(len(filelist))

    for filename in filelist:
        folder = os.path.dirname(filename)
        savefolder = folder+'/Figures/'
        
        if not os.path.isdir(savefolder):
            os.makedirs(savefolder)

        title = "_".join(savefolder.split('/')[6].split('_')[2:5])
        toto
        figs,axs = process(filename,savefolder,title=title)
        #plt.show()
        graphes.save_figs(figs,savedir=savefolder,prefix=title+'_')


def main():
    base = 'Banquise/Sebastien/Experiences/Test_frequency_h10mm_16_05_2023/white_polypropylene_10mm_d4cm_fps26p5Hz_Tchang120ms_v70/'
    folder = browse.find_path(base,disk='labshared2')    
    filelist = glob.glob(folder+'*/results/data*.pkl')
    savefolder = folder+'Results/'

    print(len(filelist))
    #toto

    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)

    #fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8))
#    F0,Phi = [],[]      
    for i,filename in enumerate(filelist):
        print(i,filename)
        #res = process(filename,savefolder)
        process(filename,savefolder)
        #Phi.append(param['phi'])
        #F0.append(param['f'])

        #F,Phi=[],[]
        #for key in res.keys():
        #    F.append(np.abs(res[key]['f']))
        #    Phi.append(res[key]['phi'])
        
        #Phi = np.unwrap(Phi)
        #ax.plot(F,Phi,'o')
        #plt.show()

#    ax.plot(F0,Phi,'ko')
    figs = {}
    figs.update(graphes.legende(r'$f$ (Hz)',r'$\Phi$ (rad)',''))
    plt.axis([0,4.5,0,16])
    plt.show()
#    graphes.save_figs(figs,savedir=savefolder,overwrite=True,prefix='Summary')

if __name__=='__main__':
    main2()
