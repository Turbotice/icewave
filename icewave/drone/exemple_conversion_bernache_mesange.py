# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 14:42:07 2025

@author: sebas
"""


import numpy as np
import os
import glob
import pandas as pd 
import h5py
import pickle
import cv2 as cv

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import matplotlib.cm as cm

import scipy.signal as signal
from scipy.interpolate import LinearNDInterpolator
from datetime import datetime,timedelta
import pytz

import icewave.tools.matlab2python as mat2py
import icewave.sebastien.set_graphs as set_graphs
import icewave.drone.drone_projection as dp


#%% FUNCTION SECTION 

def transpose_PIVmat_fields(m):
    """ Change dimensions of different fields computed using PIVlab """
    key_fields = ['Vx','Vy','Vz','X','Y']

    for key in key_fields:
        m[key] = np.transpose(m[key])
        
    return m

def synchronize_with_all_drones(m,istart = 0, iend = -1):
    """ Keep only data that is common to all drones """
    
    fields_with_time = ['Vx','Vy','Vz']
    for key in fields_with_time:
        m[key] = m[key][:,:,istart:iend]
    
    m['t'] = m['t'][istart:iend]
    
    return m

def interpolate_drone_Vz_on_trajectory(Vz,X,Y,xobj,yobj,bool_evolution_XY = 0):
    """ Interpolate a vertical velocity measured by a drone along a given tajectory along time
    Inputs : - Vz, numpy array [nx,ny,nt], time dimension should be last dimension
             - X,Y, numpy array [nx,ny] or [nx,ny,nt] if bool_evolution_XY = 1, meshgrid of (X,Y) coordinates over which Vz is known
             - xobj,yobj, numpy array (nt,), (X,Y) coordinates of the tracked object along time
             - bool_evolution_XY, boolean, equals 1 if the (X,Y) coordinates of the drone are varying with time
    Outputs : - interp_Vz, numpy array (nt,), interpolated Vz along the trajectory """
    
    # stack coordinates of meshgrid
    if not bool_evolution_XY :
        points = np.column_stack((X.ravel(),Y.ravel())) # shape (M*N,2)
        interp_Vz = np.zeros(Vz.shape[2])
        for i in range(Vz.shape[2]):
            print(f'Iteration :{i}')
            # create interpolator
            interpolator = LinearNDInterpolator(points, Vz[:,:,i].ravel())
            # interpolate Vz at time i on object position
            interp_Vz[i] = interpolator(xobj[i],yobj[i])
            
    else :
        interp_Vz = np.zeros(Vz.shape[2])
        for i in range(Vz.shape[2]):
            print(f'Iteration :{i}')
            # create interpolator
            points = np.column_stack((X[:,:,i].ravel(),Y[:,:,i].ravel())) # shape (M*N,2)
            interpolator = LinearNDInterpolator(points, Vz[:,:,i].ravel())
            # interpolate Vz at time i on object position
            interp_Vz[i] = interpolator(xobj[i],yobj[i])
        
    return interp_Vz


#%% Load PIV data from both drones 

path2data = 'Y:/Banquise/Baptiste/Resultats/Analyse_1102/Data/'
drones = ['bernache','mesange']
data = {}
for i in range(len(drones)):
    filelist = glob.glob(f'{path2data}*{drones[i]}.mat')    
    file2load = filelist[0]
    with h5py.File(file2load, 'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
    
        data[drones[i]] = mat2py.mat_to_dict(fmat['m'],fmat['m'])

#%% Load/ save synchro dictionnary
synchro = {'fulmar' : 0, 'mesange' : 520, 'bernache' : 520}

path_save = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path_save}synchro_3drones.pkl'

with open(filename,'wb') as pf:
    pickle.dump(synchro,pf)
    
# Synchronize data and structure them
drone_keys = ['bernache','mesange']
for key in drone_keys :
    data[key] = transpose_PIVmat_fields(data[key])
    istart = synchro[key]
    data[key] = synchronize_with_all_drones(data[key],istart = istart)
    
#%% Load all procruste operations 
path2procruste = 'W:/SagWin2024/Data/0211/Drones/'
procruste_file = f'{path2procruste}procruste_operation_i0_520_Nbframes_4000_any_drones.pkl'

with open(procruste_file,'rb') as pf:
    all_op = pickle.load(pf)

# keep the given set of procruste operations 
ref_drone = 'mesange' # drone chosen as reference
projected_drone = 'bernache'
key_procruste = f'{projected_drone}_2_{ref_drone}'
procruste_op = all_op[key_procruste]

# get drone GPS position and orientation 
path2parameters = 'W:/SagWin2024/Data/0211/Drones/'
file2parameters = f'{path2parameters}parameters_3drones.pkl'
with open(file2parameters,'rb') as pf :
    param = pickle.load(pf)

#%% Apply procruste_operations to bernache vertical wave field for a given frame 

Vz = data[projected_drone]['Vz']

# select frame 
i = 500

# get procruste operations
R = procruste_op['rot'][:,:,i]
translation = procruste_op['translat'][:,i]
scaling = procruste_op['scaling'][i]

# get X,Y positions of the whole field
X = data[projected_drone]['X']
Y = data[projected_drone]['Y']

# compute X,Y in the new reference frame 
X_proj,Y_proj = dp.change_XY_reference_system(X, Y, param[projected_drone], param[ref_drone], R, translation, scaling)

#%% Plot wave field 
Vz_ref = data[ref_drone]['Vz']
X_ref = data[ref_drone]['X']
Y_ref = data[ref_drone]['Y']

fig, ax = plt.subplots()
c = ax.pcolormesh(X_ref,Y_ref,Vz_ref[:,:,i],shading = 'gouraud', cmap = 'viridis',vmin = -4, vmax = 4)
c.set_rasterized(True)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_title(f'{ref_drone} in frame work of {ref_drone}')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])
fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}Vz_i{i}_{ref_drone}_framework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


fig, ax = plt.subplots()
c = ax.pcolormesh(X_proj,Y_proj,Vz[:,:,i],shading = 'gouraud', cmap = 'viridis',vmin = -4, vmax = 4)
c.set_rasterized(True)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$V_z \; \mathrm{(m.s^{-1})}$',labelpad = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$')
ax.set_ylabel(r'$Y \; \mathrm{(m)}$')
ax.set_title(f'{projected_drone} in frame work of {ref_drone}')

ax.set_aspect(1)
ax.set_xlim([-100,110])
ax.set_ylim([-60,80])
fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}Vz_i{i}_{projected_drone}_framework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Check evolution of rescaling factor with time
fig, ax = plt.subplots()
ax.plot(procruste_op['scaling'])


#%% Load buoys position dictionnary 

path2buoys = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path2buoys}structure_buoys_tracking_3drones.pkl'

with open(filename,'rb') as pf:
    Z = pickle.load(pf)

    
#%% Load buoy data 
path2buoys = 'F:/Rimouski_2024/Data/0211/BoueeVague/'
filelist = glob.glob(f'{path2buoys}B*')
buoy_files = []
for buoy_folder in filelist:
    folder = f'{buoy_folder}/mat/'
    current_file = glob.glob(f'{folder}*2000.mat')
    buoy_files.append(current_file[0])


buoy_ref = {'0':'B5','1':'B1','2':'B4'}
data_buoy = {}
for i in range(len(buoy_files)):
    file2load = buoy_files[i]
    with h5py.File(file2load, 'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
        
        key_buoy = buoy_ref[str(i)]
        data_buoy[key_buoy] = mat2py.mat_to_dict(fmat['IMU'],fmat['IMU'])


#%% Create dictionnary of buoys and Plot all buoys data

UTC_buoy = {}
Vz_buoy = {}

# build a filter for buoys
fs = 50 # sampling frequency
fc = 0.1 # cutoff frequency
order_filter = 4
b,a = signal.butter(order_filter,fc,'high',fs = fs)

for key in data_buoy.keys():
    
    UTC_buoy[key] = []
    # Build UTC time for buoys
    Y = int(data_buoy[key]['UTC_TIME']['YEAR'][0])
    M = int(data_buoy[key]['UTC_TIME']['MONTH'][0])
    D = int(data_buoy[key]['UTC_TIME']['DAY'][0])
    H = int(data_buoy[key]['UTC_TIME']['HOUR'][0])
    MIN = int(data_buoy[key]['UTC_TIME']['MIN'][0])
    SEC = int(data_buoy[key]['UTC_TIME']['SEC'][0])

    # initial UTC time for buoys
    UTC0 = datetime(Y,M,D,H,MIN,SEC)
    t0_epoch = UTC0.timestamp()
    t_epoch = t0_epoch + data_buoy[key]['IMU_DATA']['t']

    for t in t_epoch:
        if not np.isnan(t):
            UTC_buoy[key].append(datetime.fromtimestamp(t).replace(tzinfo = pytz.utc))

    mask = ~np.isnan(t_epoch)
    az = data_buoy[key]['IMU_DATA']['ACCEL_Z'][mask]
    # Apply filter to signal 
    filtered_az = signal.filtfilt(b,a,az)
    # Compute vertical velocity
    Vz_buoy[key] = np.cumsum(filtered_az - np.mean(filtered_az))/fs
    
#%% Plot all buoys signal

fig, ax = plt.subplots()
for key in Vz_buoy.keys():
    label_buoy = key
    ax.plot(UTC_buoy[key],Vz_buoy[key],label = label_buoy)
    
ax.legend()

#%% INTERPOLATION ALONG BUOY TRAJECTORY
# Select a buoy and its trajectory in frame work of ref_drone
Vz_buoy_traj = {}
Vz = {}
X = {}
Y = {}
i0 = synchro[ref_drone]

for key_drone in [ref_drone,projected_drone]: 
    Vz[key_drone] = data[key_drone]['Vz']
    X[key_drone] = data[key_drone]['X']
    Y[key_drone] = data[key_drone]['Y']

for key in data_buoy.keys():
    buoy_idx = int(key[-1])-1
    print(f'{key} : buoy_idx = {buoy_idx}')
    Vz_buoy_traj[key] = {}
    
    # Interpolate Vz from ref_drone along trajectory of buoy
    interp_Vz = {}
    for key_drone in [ref_drone,projected_drone]:
        Xbuoy = Z[key_drone]['real'][buoy_idx,0,i0:]
        Ybuoy = Z[key_drone]['real'][buoy_idx,1,i0:]
        interp_Vz[key_drone] = interpolate_drone_Vz_on_trajectory(Vz[key_drone],X[key_drone],Y[key_drone],Xbuoy,Ybuoy)
        print(f'Interpolated Vz computed for buoy {key} and drone {key_drone}')
        
    Vz_buoy_traj[key]['raw'] = interp_Vz

#%% Save dictionnary of interpolated velocities
path2save = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path2save}Interpolated_Vz_buoys_trajectories.pkl'
with open(filename,'wb') as pf:
    pickle.dump(Vz_buoy_traj,pf)

#%% Filter interpolated velocities

path2save = 'W:/SagWin2024/Data/0211/Drones/'
filename = f'{path2save}Interpolated_Vz_buoys_trajectories_ref_{ref_drone}_projected_{projected_drone}.pkl'
with open(filename,'rb') as pf:
    Vz_buoy_traj = pickle.load(pf)

# filter drone signal
fs = 30 # sampling frequency
fc = 1 # cutoff frequency
order_filter = 4
b,a = signal.butter(order_filter,fc,'low',fs = fs)

# Apply filter to signal 
for key_buoy in Vz_buoy_traj.keys():
    filtered_Vz = {}
    for key in Vz_buoy_traj[key_buoy]['raw'].keys():
        filtered_Vz[key] = signal.filtfilt(b,a,Vz_buoy_traj[key_buoy]['raw'][key])
    Vz_buoy_traj[key_buoy]['filtered'] = filtered_Vz

#%% Apply time shift to drone signals 

# UTC0 of videos
UTC0_drone = {'bernache' : datetime(2024,2,11,20,35,11,690000),
              'mesange' : datetime(2024,2,11,20,35,10,589000)}
fps = 30 # frame rate of each drone
# create new UTC0 based on drone synchronization with fulmar

for key_drone in UTC0_drone.keys():
    delta_t = synchro[key_drone]/fps
    t0 = UTC0_drone[key_drone].timestamp()
    new_t0 = datetime.fromtimestamp(t0 + delta_t).replace(tzinfo = pytz.utc)
    
    UTC0_drone[key_drone] = new_t0
    
# compute UTC_t for drones, based on UTC0 of bernache
ref_time = 'bernache'
t_drone = np.arange(0,Vz.shape[2])/fps
UTC_drone = []
for t in t_drone:
    t_epoch = UTC0_drone['bernache'].timestamp() + t
    UTC_drone.append(datetime.utcfromtimestamp(t_epoch).replace(tzinfo = pytz.utc))

UTC_drone = np.array(UTC_drone)

# apply time shift to synchronize drones and buoys
path2syncrho = 'W:/SagWin2024/Data/0211/Drones/'
file2save = f'{path_save}Buoy_drone_synchro_time_delta_i0_520.pkl'

with open(file2save,'rb') as pf:
    synchro_buoy_drone = pickle.load(pf)

UTC_drone = UTC_drone + synchro_buoy_drone['shift_timedelta']

#%% Plot superposition interpolated Vz and Vz measured by buoy

key_buoy = 'B5'
type_signal = 'filtered'

ylim = {'B5' : [-0.6,0.6],'B4': [-0.9,0.9], 'B1': [-1.5,1.5]}

for key_buoy in Vz_buoy_traj.keys():
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots(figsize = (12,9))
    for key_drone in Vz_buoy_traj[key_buoy][type_signal].keys():
        label_drone = key_drone
        ax.plot(UTC_drone,Vz_buoy_traj[key_buoy][type_signal][key_drone],'-',label = key_drone)
    ax.plot(UTC_buoy[key_buoy],Vz_buoy[key_buoy],'k-',label = key_buoy)

    ax.set_xlabel(r'UTC')
    ax.set_ylabel(r'$V_z \; \mathrm{(m.s^{-1})}$')
    ax.legend()
    
    ax.set_xlim([UTC_drone[0] , UTC_drone[-1]])
    ax.set_ylim(ylim[key_buoy])

    fig_folder = 'W:/SagWin2024/Data/0211/Drones/Buoy_Drone_superposition/Procruste_and_trajectory/'
    figname = f'{fig_folder}{key_buoy}_interpolated_Vz_ref_{ref_drone}_projected_{projected_drone}'
    plt.savefig(f'{figname}.png', bbox_inches = 'tight')
    plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')


# # stack coordinates of grid
# points_ref = np.column_stack((X_ref.ravel(),Y_ref.ravel())) # shape (M*N,2)
# points_proj = np.column_stack((X_proj.ravel(),Y_proj.ravel())) # shape (M*N,2)

# interp_Vz[ref_drone] = np.zeros((Vz_ref.shape[2]))
# interp_Vz[projected_drone] = np.zeros((Vz_ref.shape[2]))
# for i in range(Vz_ref.shape[2]):
#     print(i)
#     # create interpolator
#     interpolator_ref = LinearNDInterpolator(points_ref, Vz_ref[:,:,i].ravel())
#     interp_Vz[ref_drone][i] = interpolator_ref(Xbuoy[i],Ybuoy[i])

#     interpolator_proj = LinearNDInterpolator(points_proj, Vz[:,:,i].ravel())
#     interp_Vz[projected_drone][i] = interpolator_proj(Xbuoy[i],Ybuoy[i])








































###########################################################################
#%% ############## Syncrhonization between B5 and drones ##################
###########################################################################

key_buoy = 'B5'

# Build UTC time for buoys
Y = int(data_buoy[key_buoy]['UTC_TIME']['YEAR'][0])
M = int(data_buoy[key_buoy]['UTC_TIME']['MONTH'][0])
D = int(data_buoy[key_buoy]['UTC_TIME']['DAY'][0])
H = int(data_buoy[key_buoy]['UTC_TIME']['HOUR'][0])
MIN = int(data_buoy[key_buoy]['UTC_TIME']['MIN'][0])
SEC = int(data_buoy[key_buoy]['UTC_TIME']['SEC'][0])

# initial UTC time for buoys
UTC0 = datetime(Y,M,D,H,MIN,SEC)
t0_epoch = UTC0.timestamp()
t_epoch = t0_epoch + data_buoy[key_buoy]['IMU_DATA']['t']

UTC_t = {'buoy':[],'drone':[]}
for t in t_epoch:
    if not np.isnan(t):
        UTC_t['buoy'].append(datetime.fromtimestamp(t).replace(tzinfo = pytz.utc))

mask = ~np.isnan(t_epoch)
az = data_buoy[key_buoy]['IMU_DATA']['ACCEL_Z'][mask]

UTC0_drone = {'bernache' : datetime(2024,2,11,20,35,11,690000),
              'mesange' : datetime(2024,2,11,20,35,10,589000)}

fps = 30 # frame rate of each drone
for key_drone in UTC0_drone.keys():
    delta_t = synchro[key_drone]/fps
    t0 = UTC0_drone[key_drone].timestamp()
    new_t0 = datetime.fromtimestamp(t0 + delta_t).replace(tzinfo = pytz.utc)
    
    UTC0_drone[key_drone] = new_t0
    
# compute UTC_t for drones, based on UTC0 of bernache
ref_time = 'bernache'
fps = 30
t_drone = np.arange(0,Vz.shape[2])/fps
UTC_t['drone'] = []
for t in t_drone:
    t_epoch = UTC0_drone['bernache'].timestamp() + t
    UTC_t['drone'].append(datetime.utcfromtimestamp(t_epoch).replace(tzinfo = pytz.utc))
    
#%% Compute vertical velocity from vertical acceleration 
fs = 50 # sampling frequency
fc = 0.1 # cutoff frequency

order_filter = 4
b,a = signal.butter(order_filter,fc,'high',fs = fs)
f_array,h = signal.freqz(b,a,fs = fs)

fig,ax = plt.subplots()
ax.semilogx(f_array,20*np.log10(abs(h)))
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'Amplitude (dB)')
ax.grid(which = 'both',axis = 'both')
ax.axvline(fc,color = 'green')

# Apply filter to signal 
filtered_az = signal.filtfilt(b,a,az)

# Compute vertical velocity
Vz_buoy = np.cumsum(filtered_az - np.mean(filtered_az))/fs

#%% Superpose vertical velocities

# filter drone signal
fs = 30 # sampling frequency
fc = 1 # cutoff frequency

order_filter = 4
b,a = signal.butter(order_filter,fc,'low',fs = fs)
f_array,h = signal.freqz(b,a,fs = fs)

fig,ax = plt.subplots()
ax.semilogx(f_array,20*np.log10(abs(h)))
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'Amplitude (dB)')
ax.grid(which = 'both',axis = 'both')
ax.axvline(fc,color = 'green')

# Apply filter to signal 
filtered_Vz = {}
for key in interp_Vz.keys():
    filtered_Vz[key] = signal.filtfilt(b,a,interp_Vz[key])


#%%
fig, ax = plt.subplots()

ax.plot(UTC_t['drone'],filtered_Vz[ref_drone])
ax.plot(UTC_t['drone'],filtered_Vz[projected_drone])
ax.plot(UTC_t['buoy'],Vz_buoy,'k')

#%%

# select ROI of buoy signal 
tstart = datetime(2024,2,11,20,34,50).replace(tzinfo = pytz.utc)
tend = datetime(2024,2,11,20,37,26).replace(tzinfo = pytz.utc)

mask = [np.logical_and(t > tstart,t < tend) for t in UTC_t['buoy']]
Vz_ROI = Vz_buoy[mask]
UTC_ROI = np.array(UTC_t['buoy'])[mask]

delta_t0 = UTC_t['drone'][0] - UTC_ROI[0]
add_dt = timedelta(seconds = 18)
shifted_UTC_drone = np.array(UTC_t['drone']) - delta_t0 + add_dt

fig, ax = plt.subplots()
ax.plot(shifted_UTC_drone,filtered_Vz[ref_drone])
ax.plot(shifted_UTC_drone,filtered_Vz[projected_drone])
ax.plot(UTC_ROI,Vz_ROI,'k')

#%%

tmin = datetime(2024,2,11,20,35,46).replace(tzinfo = pytz.utc)
tmax = datetime(2024,2,11,20,36,10).replace(tzinfo = pytz.utc)

# keep buoy signal within interval
short_UTC = {}
short_Vz = {}
mask_buoy = [np.logical_and(t > tmin,t < tmax) for t in UTC_ROI]
short_UTC['buoy'] = UTC_ROI[mask_buoy]
short_Vz['buoy'] = Vz_ROI[mask_buoy]

shift = np.arange(-3000,3000,10)
dist = {ref_drone : np.zeros(shift.shape), projected_drone : np.zeros(shift.shape)}
for i in range(len(shift)):
    x = shifted_UTC_drone + timedelta(milliseconds = int(shift[i]))
    mask_drone = [np.logical_and(t > tmin,t < tmax) for t in x]
    short_UTC['drone'] = x[mask_drone]
    for key in filtered_Vz.keys():
        short_Vz[key] = filtered_Vz[key][mask_drone]
    
    x = np.array([t.timestamp() for t in short_UTC['drone']])
    xp = np.array([t.timestamp() for t in short_UTC['buoy']])
    interp_buoy = np.interp(x,xp,short_Vz['buoy'])
        
    # Compute mean square distance between each signals 

    for key in [ref_drone,projected_drone]:
        
        dist[key][i] = np.mean(np.abs(short_Vz[key] - interp_buoy))

#%% 
idx_min = {}
for key in dist.keys():
    idx_min[key] = np.argmin(dist[key])
    print(f'Minimal shift for drone {key} is : {shift[idx_min[key]]} ms')

time_precision = np.diff([shift[idx_min[key]] for key in idx_min.keys()])[0]
x = shifted_UTC_drone + timedelta(milliseconds = int(time_precision/2))
mask_drone = [np.logical_and(t > tmin,t < tmax) for t in x]
short_UTC['drone'] = x[mask_drone]
for key in filtered_Vz.keys():
    short_Vz[key] = filtered_Vz[key][mask_drone]

x = np.array([t.timestamp() for t in short_UTC['drone']])
xp = np.array([t.timestamp() for t in short_UTC['buoy']])
interp_buoy = np.interp(x,xp,short_Vz['buoy'])
        
fig, ax = plt.subplots(figsize = (12,9))
for key in [ref_drone,projected_drone]:
    label_drone = key
    ax.plot(short_UTC['drone'],short_Vz[key],label = label_drone)
ax.plot(short_UTC['drone'],interp_buoy,'k',label = 'buoy')

ax.legend()
ax.set_xlabel(r'UTC time')
ax.set_ylabel(r'$V_z \; \mathrm{(m.s^{-1})}$')

fig_folder = 'W:/SagWin2024/Data/0211/Drones/Buoy_Drone_superposition/Procruste_and_trajectory/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
figname = f'{fig_folder}Buoy_drone_synchro_ref_{ref_drone}_projected_{projected_drone}_B5_precision_{time_precision}ms'

plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')

#%% Compute total shift between two instruments and apply it on raw signals of drones 
total_shift =  add_dt + timedelta(milliseconds = int(time_precision/2)) - delta_t0

fig, ax = plt.subplots()
new_UTC_drone = np.array(UTC_t['drone']) + total_shift
ax.plot(new_UTC_drone,filtered_Vz[ref_drone])
ax.plot(new_UTC_drone,filtered_Vz[projected_drone])
ax.plot(UTC_t['buoy'],Vz_buoy,'k')

synchro_buoy_drone = {'shift_timedelta':total_shift,'shift_seconds':total_shift.total_seconds(),
                      'to_be_applied':'drone','initial_frame_drone':520}


path_save = 'W:/SagWin2024/Data/0211/Drones/'
file2save = f'{path_save}Buoy_drone_synchro_time_delta_i0_520.pkl'

with open(file2save,'wb') as pf:
    pickle.dump(synchro_buoy_drone,pf)

#%%


#%% INTERPOLATION ALONG BUOY TRAJECTORY
# Select a buoy and its trajectory in frame work of ref_drone
Vz_buoy_traj = {}
i0 = synchro[ref_drone]

Vz_ref = data[ref_drone]['Vz']
X_ref = data[ref_drone]['X']
Y_ref = data[ref_drone]['Y']

Vz_proj = data[projected_drone]['Vz']
X = data[projected_drone]['X']
Y = data[projected_drone]['Y']

################## This part is not really needed for buoys as we know there exact position on each drone ####################
X_proj = np.zeros(np.shape(Vz_proj))
Y_proj = np.zeros(np.shape(Vz_proj))
# Compute (X_proj,Y_proj) along time 
for i in range(Vz_proj.shape[2]):

    # get procruste operations
    R = procruste_op['rot'][:,:,i]
    translation = procruste_op['translat'][:,i]
    scaling = procruste_op['scaling'][i]

    # compute X,Y in the new reference frame 
    X_proj[:,:,i],Y_proj[:,:,i] = dp.change_XY_reference_system(X, Y, param[projected_drone], param[ref_drone], 
                                                                R, translation, scaling)
##############################################################################################################################

for key in data_buoy.keys():
    buoy_idx = int(key[-1])-1
    print(f'{key} : buoy_idx = {buoy_idx}')
    Vz_buoy_traj[key] = {}
    
    # Interpolate Vz from ref_drone along trajectory of buoy
    interp_Vz = {}
    Xbuoy = Z[ref_drone]['real'][buoy_idx,0,i0:]
    Ybuoy = Z[ref_drone]['real'][buoy_idx,1,i0:]
    interp_Vz[ref_drone] = interpolate_drone_Vz_on_trajectory(Vz_ref,X_ref,Y_ref,Xbuoy,Ybuoy)
    print(f'Interpolated Vz computed for buoy {key} and drone {ref_drone}')
    
    interp_Vz[projected_drone] = interpolate_drone_Vz_on_trajectory(Vz_proj,X_proj,Y_proj,Xbuoy,Ybuoy,bool_evolution_XY = 1)
    print(f'Interpolated Vz computed for buoy {key} and drone {projected_drone}')
    Vz_buoy_traj[key]['raw'] = interp_Vz

























#%% Look at Vz for a given point along time

# Coordinates of the chosen point
Xp = 0
Yp = 0

# Find indices of this point in (X,Y) ref_drone
idx = np.argmin(np.sqrt((X_ref - Xp)**2 + (Y_ref - Yp)**2))
row,col = np.unravel_index(idx, X_ref.shape)
print(X_ref[row,col],Y_ref[row,col])

# vertical velocity at the designated point
signal_ref = Vz_ref[row,col,:]

# Compute X_proj and Y_proj for all time steps
X_proj = np.zeros((np.shape(Vz)))
Y_proj = np.zeros((np.shape(Vz)))
signal_proj = np.zeros((Vz.shape[2],))

for i in range(Vz.shape[2]):
    # get procruste operations
    R = procruste_op['rot'][:,:,i]
    translation = procruste_op['translat'][:,i]
    scaling = procruste_op['scaling'][i]

    # compute X,Y in the new reference frame 
    current_X,current_Y = dp.change_XY_reference_system(X, Y, param[projected_drone], param[ref_drone], R, translation, scaling)

    # Find indices of this point in (X,Y) ref_drone
    idx = np.argmin(np.sqrt((current_X - Xp)**2 + (current_Y - Yp)**2))
    row,col = np.unravel_index(idx, current_X.shape)

    # Compute Vz 
    signal_proj[i] = Vz[row,col,i]

#%% Compare two Vz

fig, ax = plt.subplots()
ax.plot(data[ref_drone]['t'],signal_ref,label = f'{ref_drone}')
ax.plot(data[projected_drone]['t'],signal_proj,label = f'{projected_drone}')

ax.set_xlabel(r'$t \; \mathrm{(s)}$')
ax.set_ylabel(r'$V_z \; \mathrm{(m.s^{-1})}$')
ax.legend()

fig_folder = 'W:/SagWin2024/Data/0211/Drones/Resultats_f2700/'
figname = f'{fig_folder}Vz_X{Xp}_Y{Yp}_framework_{ref_drone}'
plt.savefig(f'{figname}.png', bbox_inches = 'tight')
plt.savefig(f'{figname}.pdf', bbox_inches = 'tight')






