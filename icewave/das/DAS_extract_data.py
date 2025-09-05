# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:41:26 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py 
import glob
import os 
from time import strftime, localtime
from datetime import datetime
import pytz
import time 
import scipy.signal as signal
from scipy.fftpack import fft,ifft 
from scipy.linalg import svd
import pickle

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

font_size_medium = 20
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

#%%


def epoch2datetime(t,timezone = None,format_date = '%Y-%m-%d %H:%M:%S.f'):
    """ Convert time since epoch to UTC time 
    Inputs: - t, array or list of times since epoch
            - timezone, pytz timezone in which we want datetime 
            - format_date, string, format in which datetime will be written """  
    
    t_datetime = []
    for t_epoch in t:
        current_string = strftime(format_date, localtime(t_epoch))
        current_time = datetime.strptime(current_string, format_date)
        if timezone is None:
            adjusted_time = current_time
        else:
            adjusted_time = current_time.astimezone(timezone)
            
        t_datetime.append(adjusted_time)
    
    t_datetime = np.array(t_datetime)
    t_datetime = np.reshape(t_datetime,(t.shape))
    return t_datetime

#-------------------------------------------------------------------------------------------------------


def generate_datetime_txt(UTC_0):
    """ Generate a string from a datetime UTC_0"""
    
    txt = f'{UTC_0.year}_{UTC_0.month}_{UTC_0.day}_{UTC_0.hour}_{UTC_0.minute}_{UTC_0.second}'
    return txt


#-------------------------------------------------------------------------------------------------------

def extents(f):
    """ Computes the extents of an array, returns extremities to be used with plt.imshow """
    delta = f[1] - f[0]
    return [f[0] - delta/2, f[-1] + delta/2]

#--------------------------------------------------------------------------------------------------------

def time_stacking(strain_rate,nb_seconds,fiber_length):
    """ Stack strain rate measurements over time
        Inputs : - strain_rate, 3D numpy array, #dim 0 : each second of acquisition, #dim 1 : time sampled at fs within 
        the corresponding second, #dim 2 : space
                 - nb_seconds , int type, number of seconds to stack
                 - fiber_length, float type, fiber length
                 
        Outputs : - spatio_long, 2D numpy array, #dim 0 : time in seconds, #dim 1 : space
                  - s, 1D numpy array, curvilinear coordinate
                  - t_sec, 1D numpy array, time in seconds """
                  
    time_length = np.shape(strain_rate)[0]
    fs = np.shape(strain_rate)[1] # time frequency sampling 
    fx = np.shape(strain_rate)[2]/fiber_length # space frequency sampling
    
    # curvilinear coordinate 
    s = np.arange(0,fiber_length,1/fx)
    
    for i in range(nb_seconds):
        if i == 0 :    
            spatio_long = np.vstack((strain_rate[i,:,:],strain_rate[i + 1,:,:]))
        elif i > 1:
            spatio_long = np.vstack((spatio_long,strain_rate[i,:,:]))
            
    print(f'Spatio-temp computed for a total number of seconds : {nb_seconds}')
    t_sec = np.arange(0,nb_seconds,1/fs)  # time in seconds
    
    return spatio_long,s,t_sec

#---------------------------------------------------------------------------------------------------

def new_stack_strain_rate(strain_rate,t,fiber_length,Nb_minutes = 1):
    """ Structure strain rate array so that several minutes of spatio-temporal recording are stacked in a 3D array
    Inputs: - strain_rate, 3D numpy array, #dim 0 : each second of acquisition, #dim 1 : time sampled at fs within 
            the corresponding second, #dim 2 : space
            - t, numpy array, time since epoch for each seconds of recording
            - fiber_length, float, length of the fiber
            - Nb_minutes, float, duration of each created stack
    Outputs: - stack_strain, 3D numpy array, #dim 0 : each stack of duration Nb_minutes, #dim 1 : time sampled at fs within 
            the corresponding stack, #dim 2 : space
            - stack_time, 2D numpy array # dim 0: each stack of duration Nb_minutes, # dim 1: time in seconds, taking 
            beginning of recording as a reference
            - s, numpy array, curvilinear coordinate """

    fs = np.shape(strain_rate)[1] # time frequency sampling 
    fx = np.shape(strain_rate)[2]/fiber_length # space frequency sampling
    
    # create an array # dim 0: stacks of duration Nb_minutes, 
    # dim 1: time samples within Nb_minutes, dim 2: space samples
    Nb_stacks = int(np.shape(strain_rate)[0]/60/Nb_minutes) 
    print(f'Created {Nb_stacks} stacks of duration {Nb_minutes} minutes')  
    Nb_seconds = int(Nb_minutes*60) # nb_seconds for a single section 
    
    stack_strain = np.zeros((Nb_stacks,Nb_seconds*fs,np.shape(strain_rate)[2]))
    stack_time = np.zeros((Nb_stacks,Nb_seconds*fs))
    
    for n in range(Nb_stacks):
        idx_min = int(n*Nb_seconds)
        idx_max = int((n + 1)*Nb_seconds)
        
        current_section = strain_rate[idx_min : idx_max, : ,:]
        spatio_long,s,t_sec = time_stacking(current_section,Nb_seconds,
                                            fiber_length) # stack for a given duration Nb_minutes
        stack_strain[n,:,:] = spatio_long
    
        # create array of time
        current_time = np.arange(n*Nb_seconds, (n + 1)*Nb_seconds, 1/fs)
        stack_time[n,:] = current_time 
        

    return stack_strain,stack_time,s 
    
#----------------------------------------------------------------------------------------------------

# function for plotting spatio temporal
def plot_spatio_temp(spatio,t,s,fiber_length):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - t, numpy array or list, time array 
            - s, numpy array or list, curvilinear coordinate array
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsho, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,9))
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
              interpolation = 'gaussian', extent = extents(t) + extents(s))
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

#----------------------------------------------------------------------------------------------------

def wavenumbers_hydro(freq,rho_w,rho_ice,E,h,nu,H,c_w,equation):
    """ Compute wavenumbers associated to hydroelastic waves """ 
    
    g = 9.81
    D = E*pow(h,3)/(12*(1-nu**2)) # flexural modulus

    k = np.linspace(1e-4,10,200000)
    
    idx_zero = np.zeros(len(freq)) 
    flag = 0
    for i in range(len(freq)):
        omeg = 2*np.pi*freq[i]
        if omeg == 0:
            flag = 1
            idx_flag = i
        else:
            
            if equation == 'Squire_deep':
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + 1) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Squire_shallow':
                coth = 1/np.tanh(H*k)
                func = pow(omeg, 2) * (k * h * rho_ice / rho_w + coth) - D * pow(k, 5) / rho_w - g * k
            elif equation == 'Stein':
                cph = omeg/k
                func = rho_w/D*(g-omeg/np.lib.scimath.sqrt((1/cph)**2 - (1/c_w)**2  )) - h*omeg**2*rho_w/D + pow(omeg/cph,4)
            else : 
                print('inappropriate equation name, choose between : Squire_deep / Squire_shallow / Stein')

            func[func.imag != 0] = -1
            func = func.real # keep only real part 
            print(np.where(np.diff(np.signbit(func)))[0])
            idx_zero[i] = (np.where(np.diff(np.signbit(func)))[0]) # index of the array k at which func(k) = 0
            
    idx_zero = idx_zero.astype(int)        
    k_QS =  k[idx_zero] # wave vector associated to flexural mode 
    if flag:
        k_QS[idx_flag] = 0

    return k_QS  
    

def svd_DAS(signals,fs,xs,rang,*varargin):
    """ MODIFICATIONS REQUIRED ! Computes spatio-temporal Fourier Transform using SVD decomposition.
    Arguments : 
        - signals : matrix #dim 1 : nb_recepters, #dim 2 : time, #dim 3 : nb_emitters (sources)
        - fs : sampling frequency in time 
        - xs : distance between cellphones (inverse of spatial sampling)
        - rang : row used to perform decomposition using SVD method 
        (number of singulavr values used to perform decomposition) 
        - varargin : optionnal argument
        
    Returns : 
        
        - f : array of frequencies, scaled
        - k : array of wavenumber, scaled 
        - FK : amplitude of FFT of the signal in time and space, using SVD method 
        #dim 1 : k
        #dim 2 : f """
        
    if varargin:
        if varargin[0] == 'threshold':
            threshold_user = varargin[1]
        elif varargin[0] == 'rang':
            rang = varargin[1]
        else:
            print('varargin(1) unknown')
            return
    else:
        print('varargin empty')

    Nreceiv, Nt, Nemit = signals.shape # (16, 1000, 3) (geophones, nb de valeurs, sources) # Nt not used

    # Time domain fft
    Nf = 2048
    Nk = 2048
    f = (np.arange(Nf) / Nf) * (fs if fs else 1)
    # f_axename = 'f/fs' if not fs else 'f'

    SIGNALS = fft(signals - np.mean(signals,axis = 1,keepdims = True), Nf, axis=1)
    SIGNALS = SIGNALS[:, :Nf + 1, :]

    # svd
    # Creation matrice U S V,  D???
    U = np.zeros((Nreceiv, Nf, Nemit), dtype=complex) # 16, 2048, 3
    S = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    V = np.zeros((Nemit, Nf, Nemit), dtype=complex)
    D = np.zeros((Nemit, Nf), dtype=complex)

    fig, ax = plt.subplots()

    for ii in range(Nf):
        U[:, ii, :], S[:, ii, :], V[:, ii, :] = svd(SIGNALS[:, ii, :], full_matrices=False)
        D[:, ii] = np.diag(S[:, ii, :])
    for ne in range(Nemit):
        titi = 20 * np.log10(np.abs(D[ne, :]) / np.max(np.abs(D[0, :])))
        #plt.plot(f, titi, label=f'Slice {ne}')
        ax.plot(f[1:-1],titi[1:-1],label = f'Slice {ne}')
        ax.set_xlabel('frequency (Hz)')
        ax.set_ylabel('Singular values (in dB of peak value)')

    if threshold_user is None: 
        nb_points = 5
        print(f'Select {nb_points} on the graph to define noise minimum')
        [fcut, sigmacutsup] = plt.ginput(nb_points)
        fcut = [f[0]] + fcut.tolist() + [f[-1]]
        sigmacutsup = [sigmacutsup[0]] + sigmacutsup.tolist() + [sigmacutsup[-1]]
        sigmacutsup = np.interp(f, fcut, sigmacutsup)
    else:
        sigmacutsup = np.full(int(Nf/2+1), threshold_user)
        sigmacutsup = threshold_user

    for ne in range(Nemit):
        titi = 20 * np.log10(D[ne, :] / np.max(D[0, :]))
        idx = np.where(titi <= sigmacutsup)[0]
        U[:, idx, ne] = 0

    # projection onto each singular vector

    k = (np.arange(Nk) / Nk) * (2 * np.pi / xs)  if xs else np.arange(Nk + 1)
    k_axename = 'k/ks' if not xs else 'k'

    projections = ifft(U, Nk, axis=0)#np.fft.fftshift(ifft(U, Nk, axis=0), axes=0)
    projections_sum = np.zeros((Nf, Nk, Nemit))

    for kemit in rang:
        for ii in range(Nf):
            max_value = 1  # np.max(np.abs(projections[:, ii, kemit]))
            projections_sum[ii, :, kemit] = np.abs(projections[:, ii, kemit]/max_value) ** 2

    FK = np.abs(np.mean(projections_sum, axis=2))
    
    return f,k,FK
#%% Load data 

date = '0210'

# Load parameters for DAS
path2DAS_param = 'U:/Data/parameters_Febus_2025.pkl'
with open(path2DAS_param,'rb') as pf:
    param = pickle.load(pf)
print('Parameters file loaded')

# Set parameters
# fs = np.shape(strain_rate)[1] # time sampling frequency 
# facq_x = np.shape(strain_rate)[2]/fiber_length # spatial sampling frequency
fs = param[date]['fs']
fiber_length = param[date]['fiber_length'] # fiber length in meters (set on DAS)
fiber_length = 700
facq_x = param[date]['facq_x'] 

format_date = '%Y-%m-%d %H:%M:%S.f'
local_timearea = pytz.timezone('America/Montreal')
UTC_timearea = pytz.timezone('UTC')

# path2data = f'E:/Data/{date}/DAS_h5/'
# path2data = 'F:/20250211/'
path2data = f'U:/Data/{date}/DAS/'

# Create folder for saving graphs
fig_folder = f'{path2data}Figures/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

filelist = glob.glob(path2data + '*.h5')

idx_file = 2
file2load = filelist[idx_file]
with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
    
    a_group_key = list(f.keys())[0]
    data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])
    
# shape strain rate : 
# dim 0 : time in second / dim 1 : time sampled at fs / dim 2 : space 
strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
t = data['Source1']['time'] #time since epoch


#%% Stack strain rate into sections of duration Nb_minutes
Nb_minutes = 1
stack_strain,stack_time,s = new_stack_strain_rate(strain_rate, t, fiber_length,Nb_minutes)

stack_epoch = t[0] + stack_time
UTC_stack = np.empty(stack_time.shape,dtype = 'object')
for i in range(stack_epoch.shape[0]):
    current_UTC = epoch2datetime(stack_epoch[i,:],timezone = UTC_timearea)
    UTC_stack[i,:] = current_UTC
    
#%% Show stacks 

UTC_0 = UTC_stack[0,0]

txt_UTC = generate_datetime_txt(UTC_0)
subfig_folder = f'{fig_folder}idx_file_{idx_file}_UTC_{txt_UTC}/'
if not os.path.isdir(subfig_folder):
    os.mkdir(subfig_folder)

#%%

# loop over all stacks
# for idx in range(stack_strain.shape[0]):
idx = 0
spatio = stack_strain[idx,:,:]
curent_t = stack_time[idx,:]
current_UTC = UTC_stack[idx,:]

fig,ax,imsh,cbar = plot_spatio_temp(spatio, current_UTC, s, fiber_length)
imsh.set(clim = [1, 0.3e4])
ax.set_title(f'{current_UTC[0]} - {current_UTC[-1]}')
ax.set_xlabel(r'UTC',labelpad = 5)

# initial_UTC = generate_datetime_txt(current_UTC[0])
# figname = f'{subfig_folder}spatiotemporal_duration_{Nb_minutes}_UTC_{initial_UTC}'

# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

    # plt.close('all')























###########################################################################################################
###########################################################################################################
# %% PREVIOUS STUDIES 
###########################################################################################################
###########################################################################################################

#%% Concatenate several seconds 
Nb_minutes = 1
Nb_seconds = int(Nb_minutes*60)
spatio_long,s,t_sec = time_stacking(strain_rate,Nb_seconds,fiber_length)

normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,9))
imsh = ax.imshow(spatio_long.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
          interpolation = 'gaussian', extent = extents(t_sec) + extents(s),vmin = 1, vmax = 0.6e4)
ax.set_ylim([0,fiber_length])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)


# figname = f'{fig_folder}spatio_{date}_norm_{normalization}_Nbmin_{Nb_minutes}'
# plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
# plt.savefig(f'{figname}.png',bbox_inches = 'tight')






#%% Perform 2D Fourier transform

FFT_2D,omega,k = FT.fft_2D(spatio_long - np.mean(spatio_long),[fs,facq_x],add_pow2 = [0,1])

normalization = 'linear'
fig,ax = plt.subplots()
imsh = ax.imshow(abs(FFT_2D),origin = 'lower',aspect = 'auto',norm = normalization,interpolation = 'gaussian',
          extent = extents(k) + extents(omega))
cbar = plt.colorbar(imsh)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')

ax.set_xlim([0,0.5])
ax.set_ylim([0,6])


#%% Try to perform 2D Fourier transform on two different space domains

# define space domains #[s_min,s_max]
space_domains = {'1' : {'bounds' : [0,240]}, '2' : {'bounds' : [240,700]}}

for domain in space_domains.keys():
    s_min = space_domains[domain]['bounds'][0]
    s_max = space_domains[domain]['bounds'][1]
    idx_min = np.argmin(abs(s - s_min))
    idx_max = np.argmin(abs(s - s_max))
    
    FFT_2D,omega,k = FT.fft_2D(spatio_long[:,idx_min:idx_max] - np.mean(spatio_long[:,idx_min:idx_max]),[fs,facq_x],add_pow2 = [0,1])

    space_domains[domain]['FFT_2D'] = FFT_2D
    space_domains[domain]['omega'] = omega
    space_domains[domain]['k'] = k

#%%
dom = '2'
FFT_2D = space_domains[dom]['FFT_2D']
omega = space_domains[dom]['omega']
k = space_domains[dom]['k']

# Compute theoretic curve
rho_w = 1027 # water density 
rho_ice = 917 # ice density 
E = 3.0e9 # Young modulus
nu = 0.3 # Poisson coefficient
H = 5.0 # water depth 
c_w = 1450 # sound celerity in water  

h = 0.5 # ice thickness
freq = np.linspace(0,1,100)
k_QS = wavenumbers_hydro(freq, rho_w, rho_ice, E, h, nu,H,c_w,'Squire_shallow')

# Create figure
fig, ax = plt.subplots(figsize = (12,9))
imsh = ax.imshow(abs(FFT_2D),origin = 'lower',aspect = 'auto',norm = 'log',interpolation = 'gaussian',
          extent = extents(k) + extents(omega),cmap = parula_map,vmin = 2,vmax = 2e2)
cbar = plt.colorbar(imsh)

label_theory = f'E = {E:.1e} Pa, h = {h:.2f} m, H = {H:.2f} m'
ax.plot(k_QS,2*np.pi*freq,'r-',label = label_theory)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')
ax.legend()
ax.set_xlim([0,0.5])
ax.set_ylim([0,6])

figname = f'{fig_folder}FK_{date}_spacedom_{dom}'
plt.savefig(f'{figname}.pdf',bbox_inches = 'tight')
plt.savefig(f'{figname}.png',bbox_inches = 'tight')



########################################################
############# Create stacks of strain rate #############
#%% Create new stack of strain rate

Nb_minutes = 1
Nb_seconds = Nb_minutes*60
Nb_stack = int(np.shape(strain_rate)[0]/60/Nb_minutes) # number of stacks
print(Nb_stack)

new_strain = np.zeros((Nb_stack,Nb_seconds*fs,np.shape(strain_rate)[2]))
new_time = np.zeros((Nb_stack,Nb_seconds))
for i in range(Nb_stack):
    current_strain,s,_ = time_stacking(strain_rate[i*Nb_seconds:(i+1)*Nb_seconds,:,:],Nb_seconds,fiber_length)
    new_strain[i,:,:] = current_strain


#%% 

normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,9))
current_spatio = new_strain[0,:,:]
imsh = ax.imshow(curent_spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = parula_map,
          extent = extents(t_sec) + extents(s),vmin = 1, vmax = 0.6e4)
ax.set_ylim([0,fiber_length])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$s \; \mathrm{(m)}$',labelpad = 5)








################################################
#%% Perform FFT 2D for each section
add_pow2 = [0,1]
padding= [2**(FT.nextpow2(np.shape(new_strain)[d]) + add_pow2[d - 1]) for d in range(1,3)]
FFT_stack = np.zeros((np.shape(new_strain)[0],padding[0],padding[1]),dtype = complex)

for i in range(Nb_stack):
    current_FFT2,omega,k = FT.fft_2D(new_strain[i,:,:] - np.mean(new_strain[i,:,:]),[fs,facq_x],add_pow2 = add_pow2)
    FFT_stack[i,:,:] = current_FFT2

#%%

fig, ax = plt.subplots()
mean_FFT = np.mean(abs(FFT_stack),axis = 0)
normalization = 'linear'
imsh = ax.imshow(mean_FFT,origin = 'lower',aspect = 'auto',norm = normalization,interpolation = 'gaussian',
          extent = extents(k) + extents(omega))
cbar = plt.colorbar(imsh)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$')
ax.set_ylabel(r'$\omega \; \mathrm{(rad.s^{-1})}$')

ax.set_xlim([0,0.5])
ax.set_ylim([0,6])


#%% Try SVD on several time segments of recording 

# transpose 
signals = np.transpose(new_strain,axes = (2,1,0))
rang = [0,1,2]
f,k,FK = svd_DAS(signals[:,:,:3],fs,1/facq_x,rang,'threshold', -30)
F, K = np.meshgrid(f, k)

fig, ax = plt.subplots()
imsh = ax.imshow(FK,origin = 'lower',aspect = 'auto',norm = 'linear',interpolation = 'gaussian', cmap = 'gnuplot2',
          extent = extents(k) + extents(f))






#%%%%%%%%

#%% Build low pass filter 

wn = 4 # critical frequency of low pass filter 
order_filter = 4
b,a = signal.butter(order_filter,wn,'low',fs = fs)
f_array,h = signal.freqz(b,a,fs = fs)

fig,ax = plt.subplots()
ax.semilogx(f_array,20*np.log10(abs(h)))
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'Amplitude (dB)')
ax.grid(which = 'both',axis = 'both')
ax.axvline(wn,color = 'green')

# Apply filter to signal 

filtered_spatio = signal.filtfilt(b,a,spatio_long,axis = 0)

fig,ax = plt.subplots()
imsh = ax.imshow(spatio_long.T,origin = 'lower',aspect = 'auto',norm = 'linear',
          extent = extents(t_sec) + extents(s),vmin = 1e-5, vmax = 1e5)
ax.set_ylim([450,700])
cbar = plt.colorbar(imsh)

fig,ax = plt.subplots()
imsh = ax.imshow(filtered_spatio.T,origin = 'lower',aspect = 'auto',norm = 'linear',
          extent = extents(t_sec) + extents(s),vmin = 1e-5, vmax = 1e5)
ax.set_ylim([450,700])
cbar = plt.colorbar(imsh)

#%% Look at a given position on the fiber 

s_check = 500 # position (in meter) at which we look the fiber signal 
idx = np.argmin(abs(s - s_check))

profile = spatio_long[:,idx]
profile_filtered = filtered_spatio[:,idx]


fig,ax = plt.subplots()
ax.plot(t_sec,profile)
ax.plot(t_sec,profile_filtered)


#%%

t0 = 150
idx = np.argmin(abs(t_sec - t0))
profile = spatio_long[idx,:]

fig,ax = plt.subplots()
ax.plot(s,profile)
# ax.plot(t_sec,profile_filtered)
