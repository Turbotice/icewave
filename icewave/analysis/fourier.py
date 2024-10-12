

def spectrum(y,N=200,win=False):
    n = y.shape[0]
    N = 200
    nt = int(np.floor(n/N))

    Y = np.reshape(y[:nt*N],(N,nt))

    #win = np.repeat(np.reshape(np.hanning(nt),(nt,1)),N,axis=1)
    fe = 1/dtmean
    f = np.linspace(0,fe/2,int(nt/2))

print(Y.shape)
print(nt)
TF = np.abs(np.fft.fft(Y,axis=1))
#n = int(len(TF)/2)
df = f[1]-f[0]

TF = TF[:,:int(nt/2)]/np.sqrt(df)/nt  #normalisation de la transformée de Fourier
TFmoy = np.mean(TF,axis=0)#/np.sqrt(N)

for i in range(N):
    plt.loglog(f,TF[i,:],'k')
plt.loglog(f,np.abs(TFmoy),'r-')
figs = graphes.legende('$f$ (Hz)',r'$\hat a_z~~$(m s$^{-2}/\sqrt{Hz}$)','')
#plt.loglog([100,100],[0.001,10],'r--')

figs[1]['fignum'] = 'calib_az'

plt.xlim(0.1,200)
plt.ylim(10**(-6),1)
#plt.figure()
#plt.plot(np.sum(TF**2,axis=0)*df,np.mean(Y**2,axis=0),'ko')
graphes.save_figs(figs,savedir=savefolder,suffix='_'+angle,overwrite=True)
