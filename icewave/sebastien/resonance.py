import matplotlib.pyplot as plt
import glob
import numpy as np
import os
import scipy.signal as sig
import pickle

import stephane.display.graphes as graphes
import stephane.rimouski.geophones as geo


def compute_corr(Xfilt,Zfilt,facq=26.5):
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
        ax.plot(T,C[i,:])
    
    Cmoy = np.mean(C,axis=0)

    xf,Cth = geo.fit(T,Cmoy,x0=[0.8,np.pi/2,15,3])
    print(xf)



    ax.plot(T,Cmoy,'k--')
    ax.plot(T,Cth,'ro--')
    figs={}
    figs.update(graphes.legende('$\tau$ (frames)','C',r'$\langle X(t) Z(t+\tau) \rangle$'))

    return T,Cmoy,xf,figs,ax

def filtered(data,fc=1.5,facq=26.5):
    [b,a] = sig.butter(4,fc/facq,btype='high')
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

def display_signal(data):
    fig,axs = plt.subplots(nrows=3,ncols=1,figsize=(12,8))
    axs[0].plot(data['x'],data['z'])
    graphes.legende(r'$X$ (mm)',r'$Z$ (mm)','',ax=axs[0])
    axs[1].plot(data['t'],data['x'])
    axs[1].plot(data['t'],data['z'])
    figs = graphes.legende(r'$t$ (s)',r'$X,Z$ (mm)','',ax=axs[1])
    return figs,axs
    
#plt.plot(Zfilt)
def process(filename,savefolder):
    basename = os.path.basename(filename).split('.')[0]

    data = load(filename)
    figs,axs = display_signal(data)

    X,Z = filtered(data,fc=1.5,facq=26.5)
    axs[2].plot(data['t'],X,'o-')
    axs[2].plot(data['t'],Z,'o-')
    #plt.show()
    #graphes.save_figs(figs,savedir=savefolder,overwrite=True,suffix=basename)

    T,C,xf,figs,ax = compute_corr(X,Z)
    #plt.show()
#    graphes.save_figs(figs,savedir=savefolder,overwrite=True,suffix=basename)

    return xf
    #graphes.save_figs(figs,savedir=savefolder)

def main():
    base = '/Volumes/labshared2/'
    folder = base+'Banquise/Sebastien/Experiences/Test_frequency_h10mm_16_05_2023/'
    
    filelist = glob.glob(folder+'*/results/data*.pkl')
    savefolder = folder+'Results/'
    if not os.path.isdir(savefolder):
        os.makedirs(savefolder)

    F0,Phi = [],[]      
    for i,filename in enumerate(filelist):
        print(i,filename)
        xf = process(filename,savefolder)
        Phi.append(xf[1])
        F0.append(xf[2]/(2*np.pi))

    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(8,8))

    ax.plot(F0,Phi,'ko')
    figs = {}
    figs.update(graphes.legende(r'$f$ (Hz)',r'$\Phi$ (rad)',''))
    plt.axis([0,4.5,0,np.pi])
    #plt.show()
    graphes.save_figs(figs,savedir=savefolder,overwrite=True,prefix='Summary')

if __name__=='__main__':
    main()
