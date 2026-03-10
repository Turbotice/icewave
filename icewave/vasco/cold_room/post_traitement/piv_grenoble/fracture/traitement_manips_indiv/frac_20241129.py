#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import csv
import re
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator
#matplotlib.use('TkAgg')
# %%
#%matplotlib widget
%matplotlib qt
#%%


# %% load data
W = 32
#Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 21
N = 1200
refimg = 21

date = '20241129'
name_frac_file = 'Acquisition_7'
#camera_SN = '22458101'

f_exc = 0.94
freq_acq = 15

frame_frac = 767
#computer = 'adour'
ypix_surf = 480
alpha0_deg = 17.1

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/PIV_results/{date}_frac_{name_frac_file}/'


path2data = general_folder

matfile = f'{path2data}PIV_processed_4passages_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

#matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'


from scipy.io import loadmat

mat_dict = loadmat(matfile)


u_original = mat_dict['u_original'][:,0]
v_original = mat_dict['v_original'][:,0]

def reshape_array(arr):
    array_new = np.zeros((len(arr),arr[0].shape[0],arr[0].shape[1]))
    for i in range(len(arr)):
        array_new[i,:,:] = arr[i]
    return array_new

u = reshape_array(u_original)
v = reshape_array(v_original)

u = -u
v = -v


xpix = mat_dict['x'][0][0][0]
ypix = mat_dict['y'][0][0][:,0]




#%%

yind = 10 # il vaut mieux utiliser les coordonnées en pixels

frame_plot = frame_frac - 5



plt.figure()
plt.imshow(v[frame_plot-i0],vmin=-5,vmax=5)
plt.colorbar()
plt.show()

plt.figure()
for i in range(20):
    plt.plot(v[i,yind,:])
plt.show()


plt.figure()
plt.title('frame frac = '+str(frame_frac))
for i in range(frame_plot-i0,frame_plot+10-i0): # fenetre temporelle où il y a la fracture
    plt.plot(v[i,yind,:],label=str(i+i0))
plt.legend()
plt.show()


# %% courbure de la plaque juste avant la fracture ?

dcm_sur_dpx = 0.07
dcm_sur_dpx_err = 0.005

xvals = np.arange(v.shape[2]) * dcm_sur_dpx * W/2 * 1e-2
v_converted_meters = v * dcm_sur_dpx * 1e-2

# 0.075 est la valeur approximative de dcm/dpx, mais une utilisation de la variation de dcm/dpx 
# est requise pour une meilleure estimation de kappa_c

if W==32:
    ind_inf_fit = 55
    ind_sup_fit = 85
else:
    ind_inf_fit = 20
    ind_sup_fit = 45

print(xvals)
fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
print(fit_params)

plt.figure()
plt.title('profile at yind='+str(yind))
plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='frame '+str(frame_frac))
plt.plot(xvals[ind_inf_fit:ind_sup_fit],(fit_params[0]*xvals**2 + fit_params[1]*xvals + fit_params[2])[ind_inf_fit:ind_sup_fit],label='fit : $\kappa$ = '+str(np.round(2*fit_params[0],3)))
plt.legend()
plt.show()

kappa_c = 2*fit_params[0]


#%% calcul echelles avec photos regle
tab_dcm = np.array([23,30,30])
tab_dpx = np.array([259,384,497])

delta_y = 1100-200 # car photos regles pas memes dimensions que photos fracture...
tab_y = np.array([1100,1185,1350]) - delta_y
plt.figure()
plt.plot(tab_y,tab_dcm/tab_dpx,'o')
popt_echelle,_ = curve_fit(lambda x,a,b:a*x + b,tab_y,tab_dcm/tab_dpx)
plt.plot(tab_y,popt_echelle[0]*tab_y+popt_echelle[1])
plt.show()

def compute_aspect_ratio_simple(y,popt=popt_echelle):
    return popt[0]*y+popt[1]

#%% calcul de tab_dcmsurdpx pour toutes les positions de cases piv (x,y) (avec x et y en pixels)

XPIX,YPIX = np.meshgrid(xpix,ypix)

DCM_SUR_DPX_2 = np.zeros((v.shape[1],v.shape[2]))
for j in range(len(ypix)):
    DCM_SUR_DPX_2[j,:] = compute_aspect_ratio_simple(ypix[j],popt=popt_echelle)

plt.figure()
plt.imshow(DCM_SUR_DPX_2)
plt.colorbar()
plt.show()

DcmSurDpx_3dim = np.tile(DCM_SUR_DPX_2,(u.shape[0],1,1))

v_cm = DcmSurDpx_3dim * v
v_converted_meters = v_cm * 1e-2



xvals_px = np.arange(v.shape[2]) * W/2
xvals_test =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
v_converted_meters_1 = v * 0.07 * 1e-2


plt.figure()
#plt.plot(xvals,v[frame_frac-i0,yind,:]*0.075,label='frame '+str(frame_frac))
plt.plot(xvals_test,v_cm[frame_frac-i0,yind,:],label='frame '+str(frame_frac))

plt.show()


# %% affichage de differents profils avec correction angulaire vs y 
# et variation echelle horizontale vs y
"""
sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/cold_room/post_traitement/piv_grenoble/fracture/python_functions/')

from spatial_scale import *

array_alphas = np.zeros(len(ypix))

for i in range(len(ypix)):
    array_alphas[i] = compute_angle(ypix[i],ypix_surf, dcm_sur_dpx_vertical=47/250, H_cam= 45, d_cam_vitre=87)
array_alphas = np.reshape(array_alphas,(len(array_alphas),1))
ArrAlph = np.tile(array_alphas, (v.shape[0],1,v.shape[2]))
#print(ArrAlph[0,0,:])
#print(ArrAlph[0,:,0])
v_angle_corrected = v_converted_meters/np.cos(ArrAlph)
"""
v_angle_corrected = (1/np.cos(np.radians(alpha0_deg))) * v_converted_meters


kappa_c_vals = []

yindices = np.array([6,9,10,11,12,13,14])

#ind_inf_fit = 20
#ind_sup_fit = 42

plt.figure()
for yind in yindices:
    #print(xvals)
    xvals =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
    xvals_fit = xvals[ind_inf_fit:ind_sup_fit]
    fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_angle_corrected[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
    #print(fit_params)
#    fit_params_2 = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
#    plt.figure()
    plt.title('profiles at frame '+str(frame_frac))
#    plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals,v_angle_corrected[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals_fit,fit_params[0]*xvals_fit**2+fit_params[1]*xvals_fit+fit_params[2],color='red',alpha=0.5)
    plt.plot()
    plt.legend()
#    plt.show()

    kappa_c = 2*fit_params[0]
    print('kappac',kappa_c)
    kappa_c_vals.append(kappa_c)

plt.show()

kappa_c_vals = np.array(kappa_c_vals)
print(np.mean(np.array(kappa_c_vals)))
print(np.std(np.array(kappa_c_vals)))



# %% enregistrement des données
file_dict_results = 'R:/Gre25/Summary/fracture_postprocessing/resultats/fracture_results.pkl'
"""
if os.path.exists(file_dict_results):
    with open(file_dict_results, 'rb') as f:
        dict_results = pickle.load(f)
else:
    dict_results = {}



if date in dict_results:
    dict_results[date][name_frac_file] = {'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf, 'f_exc':f_exc}
else:
    dict_results[date] = {name_frac_file: {'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf,'f_exc':f_exc}}


"""

if os.path.exists(file_dict_results):
    with open(file_dict_results, 'rb') as f:
        dict_results = pickle.load(f)
else:
    dict_results = {}



if date in dict_results:
    dict_results[date][name_frac_file] = {'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf,'f_exc':f_exc}
else:
    dict_results[date] = {name_frac_file: {'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf,'f_exc':f_exc}}


ee = input('Are you sure you want to erase last results ?(y/n)')
if ee=='y':
    with open(file_dict_results, 'wb') as handle:
        pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    pass

# %%
