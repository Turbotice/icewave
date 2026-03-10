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
from scipy.optimize import curve_fit
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

savgoltest = True

#date = '20241125'
#f_exc =500
#freq_acq = 124.69
f_exc = 10
freq_acq = 9.901
#f_exc = 20
#freq_acq = 19.802

date = '20241126'



acq_num = 1
camera_SN = '40300722'

#W = 32
#Dt = 20
W = 64
Dt = 50

#general_folder = 'F:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
#general_folder = f'D:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_1/camera_{camera_SN}/'

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
index_profile_line = 17
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


#%%


X_MAXS = np.arange(matrix.shape[1])
Y_MAXS = []

matrix = np.flip(matrix,axis=0)

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
xlim_fit = (0,30) # à modifier en fonciton de ce qu'on a!
#xlim_fit = (0,50) # à modifier en fonciton de ce qu'on a!
#W = 64
dcm = 3.2
dpx = 269
x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
x_meters = x_px * W/2 * (dcm*1e-2)/dpx

x_meters = np.unwrap(x_meters)

t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
print(t_px)
t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


def linear(x,a,b):
    return a*x+b
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

def routine_vitesse_phase(f_exc=10, freq_acq=9.901, date='20241125',acq_num=1,camera_SN='40300722', dosavgol=True,W=32,Dt=20,plot=True,index_profile_line=34, xlim_fit = (75,90), dcm = 3.2, dpx = 269):
    ## load data
    #general_folder = 'F:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
    general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
    #general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

    folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
    matfile_path = folder_video_demod + '/figdata_complex.mat'
    print('data loading...')
    Data_demod = load_complex_field(matfile_path)
    print('data loaded!')
    complex_field = np.copy(Data_demod.data)
    if plot:
        plt.figure()
        plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field)))
        plt.show()
    ## process
    nb_periods = 4
    nb_pts_per_period = 256
    #index_profile_line = 34 # pour piv du 25/11 c'est 34 (W=32), pour W=64 on met donc 17
    m = nb_periods*nb_pts_per_period
    matrix = np.zeros((m,complex_field.shape[1]))
    for i in range(m):
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))#/np.max(np.abs(complex_field))
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        if dosavgol==True:
            img = savgol_filter(img,21,3)
        #img = np.real(complex_field*np.exp(2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        p = img[index_profile_line,:]
        #plt.imshow(img)
        #plt.pause(0.01)
        matrix[i,:] = p


    # calcul du spatio-temporel du champ démodulé en temps


    X_MAXS = np.arange(matrix.shape[1])
    Y_MAXS = []

    matrix = np.flip(matrix,axis=0)

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

            
    ## plot the line to fit on top of the matrix imshow (test)
    if plot:
        plt.figure()
        plt.imshow(matrix,aspect='auto')
        plt.plot(Y_MAXS,'r')
        plt.colorbar()
        plt.xlabel('position along profile [PIV box units]',fontsize=15)
        plt.ylabel('time [256 px = T = 1/$f_{exc}$]',fontsize=15)
        #plt.ylim(np.min(Y_MAXS),np.max(Y_MAXS))
        plt.show()

    ##  fit the line on the spatio-temporal diagram

    x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
    x_meters = x_px * W/2 * (dcm*1e-2)/dpx

    x_meters = np.unwrap(x_meters)

    t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
    print(t_px)
    t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


    def linear(x,a,b):
        return a*x+b
    opt,pcov = curve_fit(linear,t_sec,x_meters)
    plt.plot(t_sec,x_meters,'o')
    t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
    plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Position of the wave front (meters)',fontsize=15)
    plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(pcov[0][0],2))+' m/s',fontsize=15)
    plt.show()
    print('phase velocity = '+str(opt[0])+' +- '+str(pcov[0][0])+' m/s')

    return opt[0],pcov[0][0]

#%% changer generalfolder si besoin
date = '20241126'

acq_num = 2
camera_SN = '40300722'

#W = 32
#Dt = 20
W = 64
Dt = 50

general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_1/camera_{camera_SN}/'

#%% lister toutes les frequences

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

#%% lancer routine
tab_v_phase = np.zeros(len(tab_f_exc))
tab_v_phase_err = np.zeros(len(tab_f_exc))
for i in range(len(tab_f_exc)):
    tab_v_phase[i],tab_v_phase_err[i] = routine_vitesse_phase(f_exc=tab_f_exc[i],freq_acq=tab_freq_acq[i], date='20241126',acq_num=2,W=64,Dt=50,index_profile_line=17,xlim_fit=(0,30),camera_SN='40300722')#'40437120')

#%%

#  relation dispersion
def compute_omega(k,D):
    return np.sqrt(((D/rho)*k**5 + (T/rho)*k**3 + g*k)*np.tanh(k*H))

def compute_D(e):
    return (E*(e**3)/(12*(1-nu**2)))

e = 4e-3 
H = 14.5e-2 - e 
E = 9e9
nu = 0.4
D_value = compute_D(e)
rho = 1e3
T = 0
lambdamin = 0.5e-1
lambdamax = 1
g = 10
tab_lambda = np.linspace(lambdamin,lambdamax,10000)

#mask = (np.array(dico_dataset['v_phase_err'])<(1/10)*np.array(dico_dataset['v_phase']))&(np.array(dico_dataset['v_phase'])>0)

plt.figure()
plt.plot((1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value),tab_lambda * (1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value),label='E=9GPa')
E = 4e9
D_value = compute_D(e)
plt.plot((1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value),tab_lambda * (1/(2*np.pi))*compute_omega(2*np.pi/tab_lambda,D_value),label='E=4GPa')

#plt.errorbar(np.array(dico_dataset['tab_f_exc'])[mask],np.array(dico_dataset['v_phase'])[mask],np.array(dico_dataset['v_phase_err'])[mask],marker='o')
plt.errorbar(tab_f_exc,tab_v_phase,tab_v_phase_err,marker='o',linestyle='')
plt.ylim(0,80)
plt.xlim(0,300)
plt.xlabel('frequency (Hz)',fontsize=15)
plt.ylabel('$v_{\phi}$',fontsize=15)
plt.legend()
plt.show()

# %%
