#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import re
from scipy.signal import savgol_filter
#matplotlib.use('TkAgg')
# %% def fonctions
def savedict(filename,dico):
    with open(filename+'.pickle', 'wb') as handle:
        pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
def loaddict(filename):
    with open(filename+'.pickle', 'rb') as handle:
        dico = pickle.load(handle)
    return dico
def extract_data_mat_file(path_mat_file):
    with h5py.File(path_mat_file, 'r') as file:
        data = file['data'][:]  # Accessing 'data' variable
    return data
#    return data

def load_complex_field(path_mat_file):
    matrix = extract_data_mat_file(path_mat_file)
    # Initialize an empty matrix to store complex numbers
    complex_matrix = np.empty(matrix.shape, dtype=complex)
    
    # Iterate through each element in the matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            # Get the tuple from the original matrix
            a, b = matrix[i, j]
            
            # Create a complex number from the tuple and assign it to the corresponding position
            complex_matrix[i, j] = complex(a, b)
    return complex_matrix    
#%% import test data

savgoltest = False

date = '20240604'
#f_exc =500
#freq_acq = 124.69
f_exc =300
freq_acq = 149.25

W = 32
Dt = 8

general_folder = 'X:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
matfile_path = folder_video_demod + '/figdata_complex.mat'
print('data loading...')
Data_demod = load_complex_field(matfile_path)
print('data loaded!')
complex_field = np.copy(Data_demod.data)
#%% montrer le champ démodulé test
#print(complex_field)
plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field)))
#for i in range(48):
#    img = np.real(complex_field*np.exp(2*np.pi*1j*i/48))/np.max(np.abs(complex_field))
#    plt.imshow(img)
#    plt.pause(0.01)
plt.show()
#%% process test
nb_periods = 4
nb_pts_per_period = 256
index_profile_line = 4
m = nb_periods*nb_pts_per_period
matrix = np.zeros((m,complex_field.shape[1]))
for i in range(m):
    img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))#/np.max(np.abs(complex_field))
    img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
    if savgoltest==True:
        img = savgol_filter(img,21,3)
    #img = np.real(complex_field*np.exp(2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
    p = img[index_profile_line,:]
    #plt.imshow(img)
    #plt.pause(0.01)
    matrix[i,:] = p

X_MAXS = np.arange(matrix.shape[1])
Y_MAXS = []
for i in range(matrix.shape[1]):
    indices = find_peaks(matrix[:,i])[0]
    if len(indices)==3:
        indices=np.hstack((indices,np.array([np.inf])))
#        i0 = indices[0]
#        i1 = indices[1]
#        i2 = indices[2]
    #i3 = indices[3]        
#        tab_max_along_column = np.array([i0,i1,i2])

    if len(Y_MAXS)==0:
        Y_MAXS.append(i1)
    else:
        index_closest = np.argmin(abs(Y_MAXS[-1] - indices))
        Y_MAXS.append(indices[index_closest])
"""X_MAXS = np.arange(matrix.shape[1])
Y_MAXS = []

for i in range(matrix.shape[1]):
    indices = find_peaks(matrix[:,i])[0]
    i0 = indices[0]
    i1 = indices[1]
    if len(Y_MAXS)==0:
        Y_MAXS.append(i0)
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=nb_pts_per_period/3)):
    #    Y_MAXS.append(i0)
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>nb_pts_per_period/3)):
    #    Y_MAXS.append(i1)
    elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<abs(Y_MAXS[-1]-i1))):
        Y_MAXS.append(i0)
    elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>abs(Y_MAXS[-1]-i1))):
        Y_MAXS.append(i1)
    else:
        Y_MAXS.append(np.nan)
"""        
#%% plot the line to fit on top of the matrix imshow (test)
plt.figure()
plt.imshow(matrix,aspect='auto')
plt.plot(Y_MAXS,'r')
plt.colorbar()
plt.xlabel('position along profile [PIV box units]',fontsize=15)
plt.ylabel('time [256 px = T = 1/$f_{exc}$]',fontsize=15)
#plt.ylim(np.min(Y_MAXS),np.max(Y_MAXS))
plt.show()

# %% fit the line (test)
xlim_fit = (10,len(Y_MAXS)) # à modifier en fonciton de ce qu'on a!
#xlim_fit = (0,50) # à modifier en fonciton de ce qu'on a!
#W = 64
dcm = 10
dpx = 501
x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
x_meters = x_px * W/2 * (dcm*1e-2)/dpx

x_meters = np.unwrap(x_meters)

t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
print(t_px)
t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


def linear(x,a,b):
    return a*x+b
from scipy.optimize import curve_fit
opt,pcov = curve_fit(linear,t_sec,x_meters)
plt.plot(t_sec,x_meters,'o')
t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
plt.xlabel('Time (seconds)',fontsize=15)
plt.ylabel('Position of the wave front (meters)',fontsize=15)
plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(pcov[0][0],2))+' m/s',fontsize=15)
plt.show()
print('phase velocity = '+str(opt[0])+' +- '+str(pcov[0][0])+' m/s')
#fu.linfitxy(t_sec,x_meters,plot=True)

# %% maintenant on veut implementer ca dans une sorte de routine

#%% Dictionary business

date = '20240604'
#date = '20240531'

dataset = 0 # si dataset=0 alors on va dans les dossier matData sans numéro à la fin

general_folder = 'X:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
#general_folder = 'I:/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
pathdict = general_folder + 'DICO' # la focntion save va rajouter .pickle après
if os.path.exists(pathdict+'.pickle'):
    DICO = loaddict(pathdict)
else:
    DICO = {}
#test
#DICO['hello'] = 'wooooorld'
#savedict(namedict,DICO)
# rentrer les tableaux de fréquences dans le dico
if 'dataset'+str(dataset) in DICO:
    pass
else:
    DICO['dataset'+str(dataset)] = {}
    ll = os.listdir(general_folder)
    ld = []
    list_f_exc = []
    list_freq_acq = []
    for l in ll:
        if ('Hz_' in l)&(l[-2:]=='Hz'):
            ld.append(l)
            print(l)
            list_f_exc.append(int(l[:l.find('Hz')]))
            if '.' in l[l.find('_')+1:-2]:
                list_freq_acq.append(float(l[l.find('_')+1:-2]))
            elif ('.' in l[l.find('_')+1:-2])==False:
                list_freq_acq.append(int(l[l.find('_')+1:-2]))
tab_f_exc = np.array(list_f_exc)
tab_freq_acq = np.array(list_freq_acq,dtype='object')
indsort = np.argsort(tab_f_exc)
tab_f_exc = tab_f_exc[indsort]
tab_freq_acq = tab_freq_acq[indsort]

dico_dataset = DICO['dataset'+str(dataset)]
#dico_dataset['tab_f_exc'] = tab_f_exc
#dico_dataset['tab_freq_acq'] = tab_freq_acq

#dico_dataset['Complex_Fields'] = {}
#DCF = dico_dataset['Complex_Fields'] # DCF pour "Dico_Complex_Fields"

#dico_dataset['W_values'] = {}
#DWV = dico_dataset['W_values'] # Dico W values

#dico_dataset['Dt_values'] = {}
#DDtV = dico_dataset['Dt_values'] # Dico Dt values


#for i in range(len(tab_f_exc)):
#        f_exc = tab_f_exc[i]
#        freq_acq = tab_freq_acq[i]
#        if f_exc<200:
#            Dt = 8
#        elif (f_exc>=200)&(f_exc<500):
#            Dt = 4
#        elif (f_exc>=500):
#            Dt = 2

        #matfile_path = folder_video_demod + '/figdata_complex.mat'
        #print('data loading... '+str(f_exc))
        #Data_demod = load_complex_field(matfile_path)
        #print('data loaded!')
        #complex_field = Data_demod.data
        #DCF[str(f_exc)] = complex_field
        #DWV[str(f_exc)] = W
        #DDtV[str(f_exc)] = Dt

#%% routine "à moitié à la mano"

#index_data_in_dataset = 10
#f_exc = tab_f_exc[index_data_in_dataset]
#freq_acq = tab_freq_acq[index_data_in_dataset]
#W = DWV[str(f_exc)]
#Dt = DDtV[str(f_exc)]

do_savgol = True

remove_left_part = False
uncertainty_pixels=False
#f_exc = 250
#freq_acq = 124.38

tab_Dt = np.array([16,16,16,16,16,16,16,16,16,16,16,8,8,8,4])

for i in range(len(tab_f_exc)):

    f_exc = tab_f_exc[i]
    freq_acq = tab_freq_acq[i]

    W = 64
    #Dt = 4
    Dt = tab_Dt[i]
    if not('v_phase' in dico_dataset):
        dico_dataset['v_phase'] = []
        dico_dataset['v_phase_err'] = []
        dico_dataset['xlim_fit_px'] = []
        dico_dataset['tab_f_exc'] = []
        dico_dataset['tab_freq_acq'] = []

    if dataset==0:
        folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
    elif dataset!=0:
        folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData'+str(dataset)+'/video_demod_W'+ str(W) +'_Dt' +str(Dt)

    matfile_path = folder_video_demod + '/figdata_complex.mat'
    print('data loading... '+str(f_exc))
    Data_demod = load_complex_field(matfile_path)
    print('data loaded!')
    complex_field = np.copy(Data_demod.data)
    #complex_field = DCF[str(f_exc)] # affect the corresponding complex field to the variable complex_field
    print('imshow of the real part of complex field for f_exc = '+str(f_exc)+' Hz')
    plt.figure()
    plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field)))

    plt.show()

    nb_periods = 4
    nb_pts_per_period = 256
    index_profile_line = 4
    m = nb_periods*nb_pts_per_period
    matrix = np.zeros((m,complex_field.shape[1]))
    if remove_left_part==True:
        matrix = np.zeros((m,complex_field.shape[1]-5))
    for i in range(m):
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.max(np.abs(complex_field))
        #img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        if remove_left_part==True:
            img = img[:,5:]
        if do_savgol==True:
            img = savgol_filter(img,21,3)
        p = img[index_profile_line,:]
        #plt.imshow(img)
        #plt.pause(0.01)
        matrix[i,:] = p

    print('profile vs time at line '+str(index_profile_line)+' computed.')


    X_MAXS = np.arange(matrix.shape[1])
    Y_MAXS = []
    for i in range(matrix.shape[1]):
        indices = find_peaks(matrix[:,i])[0]
        if len(indices)==3:
            indices=np.hstack((indices,np.array([np.inf])))
#        i0 = indices[0]
#        i1 = indices[1]
#        i2 = indices[2]
        #i3 = indices[3]        
#        tab_max_along_column = np.array([i0,i1,i2])

        if len(Y_MAXS)==0:
            Y_MAXS.append(i1)
        else:
            index_closest = np.argmin(abs(Y_MAXS[-1] - indices))
            Y_MAXS.append(indices[index_closest])
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=nb_pts_per_period/3)):
    #    Y_MAXS.append(i0)
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>nb_pts_per_period/3)):
    #    Y_MAXS.append(i1)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<abs(Y_MAXS[-1]-i1))):
        #    Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>abs(Y_MAXS[-1]-i1))):
        #    Y_MAXS.append(i1)
        #else:
        #    Y_MAXS.append(np.nan)
    print('indices maxima (using find_peaks) found.')
    print('plot peaks on top of imshow of profile vs time : ')

    plt.figure()
    plt.imshow(matrix,aspect='auto')
    plt.plot(Y_MAXS,'r')
    plt.colorbar()
    #plt.ylim(np.min(Y_MAXS),np.max(Y_MAXS))
    plt.show()
    



    #  fitter la partie d'intérêt
    
    
    
    xlim_fit = (10,33) # à modifier en fonciton de ce qu'on a!
    #xlim_fit = (1,20)
    W = 64
#    dcm = 11
#    dpx =  978
    dcm = 10
    dpx =  501

    x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
    x_meters = x_px * (W/2) * (dcm*1e-2)/dpx
    t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
    print(t_px)
    t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


    def linear(x,a,b):
        return a*x+b
    from scipy.optimize import curve_fit
    if uncertainty_pixels==True:
        opt,pcov = curve_fit(linear,t_sec,x_meters,sigma=np.ones(len(x_meters))*(48/2)* W/2 * (dcm*1e-2)/dpx)
    else:
        opt,pcov = curve_fit(linear,t_sec,x_meters)
    plt.plot(t_sec,x_meters,'o')
    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Position of the wave front (meters)',fontsize=15)
    plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(pcov[0][0],2))+' m/s',fontsize=15)

    #plt.errorbar(t_sec,x_meters,np.ones(len(x_meters))*(48/2)* W/2 * (dcm*1e-2)/dpx)
    t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
    plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
    plt.show()
    print('phase velocity = '+str(opt[0])+' +- '+str(pcov[0][0])+' m/s')
    
    
    
    #  affecter valeurs fittées dans dictionnaire
    
    
    
    #dico_dataset['v_phase'][str(f_exc)] = opt[0]
    #dico_dataset['v_phase_err'][str(f_exc)] = pcov[0][0]
    #dico_dataset['xlim_fit_px'][str(f_exc)] = xlim_fit

    dico_dataset['tab_f_exc'].append(f_exc)
    dico_dataset['tab_freq_acq'].append(freq_acq)
    dico_dataset['v_phase'].append(opt[0])
    dico_dataset['v_phase_err'].append(pcov[0][0])
    dico_dataset['xlim_fit_px'].append(xlim_fit)
# %% A LA FIN : enregistrer le dico

savedict(filename=pathdict,dico=DICO)

# %% relation dispersion
def compute_omega(k,D):
    return np.sqrt(((D/rho)*k**5 + (T/rho)*k**3 + g*k)*np.tanh(k*H))

def compute_D(e):
    return (E*(e**3)/(12*(1-nu**2)))

e = 1.5e-3 
H = 2.5e-2 - e 
E = 9e9
nu = 0.4
D_value = compute_D(e)
rho = 1e3
T = 0
lambdamin = 0.5e-1
lambdamax = 1
g = 10
tab_lambda = np.linspace(lambdamin,lambdamax,10000)

mask = (np.array(dico_dataset['v_phase_err'])<(1/10)*np.array(dico_dataset['v_phase']))&(np.array(dico_dataset['v_phase'])>0)

plt.figure()
plt.plot((1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value),tab_lambda * (1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value))
plt.errorbar(np.array(dico_dataset['tab_f_exc'])[mask],np.array(dico_dataset['v_phase'])[mask],np.array(dico_dataset['v_phase_err'])[mask],marker='o')
plt.ylim(0,50)
plt.xlim(0,500)
plt.xlabel('frequency (Hz)',fontsize=15)
plt.ylabel('$v_{\phi}$',fontsize=15)
plt.errorbar([60],[10.47],[0.15],marker='o') # obtenu sur une autre ligne de profil (6 au lieu de 4)
plt.errorbar([150],[18.40],[0.56],marker='o') # obtenu sur une autre ligne de profil (6 au lieu de 4)
plt.errorbar([250],[22.5],[0.51],marker='o')
plt.show()
# %%
