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
i0 = 0
N = 5540
Dt = 1

date = '0514'
name_frac_file = 'img_seq1'
#camera_SN = '22458101'

f_exc = 0.9042
freq_acq = 20

frame_frac = 5095
#computer = 'adour'
ypix_surf = 435
alpha0_deg = 15 # pas mesuré exactement

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

matfile = f'{path2data}PIV_processed_4passages_i0{i0}_N{N}_W{W}_Dt{Dt}.mat'

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

yind = 3 # il vaut mieux utiliser les coordonnées en pixels

frame_plot = frame_frac - 10

plt.figure()
plt.imshow(v[frame_plot-i0],vmin=-5,vmax=5)
plt.colorbar()
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


ind_inf_fit = 20
ind_sup_fit = 40
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

#%%load file echelles fracture pour avoir en x et y les variations de dpx/dcm
#if computer=='Leyre':
#    file_echelles_fracture = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/echelles/echelles_fracture.txt'
if system_loc=='windows_server':
    file_echelles_fracture = 'R:/Gre25/Data/0512/cameras/ref_matin/echelles.txt'
elif system_loc=='linux_server':
    file_echelles_fracture = '/media/turbots/GreDisk/Gre25/Data/0512/cameras/ref_matin/echelles.txt'
data_ech_frac = np.loadtxt(file_echelles_fracture,skiprows=1,usecols=range(5))

d = {}
"""
d['xmoy'] = data_ech_frac[:,1]
d['ymoy'] = data_ech_frac[:,2]
d['dcm'] = data_ech_frac[:,3]
d['dpx'] = data_ech_frac[:,4]
d['delta_y'] = data_ech_frac[:,5]

d['tab_ymoy_refmanip'] = d['ymoy'] - d['delta_y']
"""
d['xmoy'] = data_ech_frac[:,0]
d['ymoy'] = data_ech_frac[:,1]
d['dcm'] = data_ech_frac[:,2]
d['dpx'] = data_ech_frac[:,3]
d['ypix_surf_ref'] = data_ech_frac[:,4]

delta_y = d['ypix_surf_ref'] - ypix_surf

d['tab_ymoy_refmanip'] = d['ymoy'] - delta_y



rng = np.random.default_rng()
x = d['xmoy']
y = d['tab_ymoy_refmanip']
z = d['dcm']/d['dpx']
X = np.linspace(min(x), max(x))
Y = np.linspace(min(y), max(y))
X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
interp_function = LinearNDInterpolator(list(zip(x, y)), z)

d['interp_function'] = interp_function

Z = interp_function(X, Y)
plt.pcolormesh(X, Y, Z, shading='auto')
plt.plot(x, y, "ok", label="input point")
plt.legend()
plt.colorbar()
plt.axis("equal")
plt.xlabel('xmoy')
plt.ylabel('ymoy')
plt.xlim(np.min(x)-50,np.max(x)+50)
plt.show()

def linfitinterp(x,d=d,plot=False):
    interp_function = d['interp_function']
    tab_y = np.linspace(np.min(d['tab_ymoy_refmanip']), np.max(d['tab_ymoy_refmanip']))
    X,Y = np.meshgrid(x,tab_y)
    Z = interp_function(X,Y)
    if np.sum(np.isnan(Z)==False)==0:
        return np.zeros(2)*np.nan,np.zeros((2,2))*np.nan
    popt,pcov = curve_fit(lambda x,a,b:a*x + b,Y[np.isnan(Z)==False],Z[np.isnan(Z)==False])
    if plot:
        plt.figure()
        plt.plot(Y.flatten(),Z.flatten(),'o')
        plt.plot(Y.flatten(),popt[0]*Y.flatten()+popt[1])
        plt.show()
    return popt,pcov


def compute_aspect_ratio(x,y,d=d,plot=False):
    popt,_ = linfitinterp(x,d=d,plot=plot)
    dcm_sur_dpx = popt[0]*y+popt[1]
    return dcm_sur_dpx    

#%% calcul de tab_dcmsurdpx pour toutes les positions de cases piv (x,y) (avec x et y en pixels)

XPIX,YPIX = np.meshgrid(xpix,ypix)

DCM_SUR_DPX = interp_function(XPIX,YPIX)

plt.figure()
plt.imshow(DCM_SUR_DPX)
plt.show()

DCM_SUR_DPX_2 = np.zeros((DCM_SUR_DPX.shape[0],DCM_SUR_DPX.shape[1]))
for j in range(len(ypix)):
    for i in range(len(xpix)):
        try:
            DCM_SUR_DPX_2[j,i] = compute_aspect_ratio(xpix[i],ypix[j],d=d)
        except:
            DCM_SUR_DPX_2[j,i] = np.nan

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

#%% hidden cell
# 0.075 est la valeur approximative de dcm/dpx, mais une utilisation de la variation de dcm/dpx 
# est requise pour une meilleure estimation de kappa_c
"""
plt.figure()

for yind in [3,6]:
    print(xvals)
    if yind==3:
        factor=0.9
    else:
        factor=1
    fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],factor*v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
    print(fit_params)

#    plt.figure()
    plt.title('profile at yind='+str(yind))
    plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='frame '+str(frame_frac),color='tab:blue')
    plt.plot(xvals,v_converted_meters_1[frame_frac-i0,yind,:],label='frame '+str(frame_frac),color='tab:orange')
    plt.plot(xvals[ind_inf_fit:ind_sup_fit],(fit_params[0]*xvals**2 + fit_params[1]*xvals + fit_params[2])[ind_inf_fit:ind_sup_fit],label='fit : $\kappa$ = '+str(np.round(2*fit_params[0],3)),color='tab:green')
    plt.legend()
#    plt.show()

    kappa_c = 2*fit_params[0]
    print(kappa_c)
"""
# %% affichage de differents profils avec correction angulaire vs y 
# et variation echelle horizontale vs y
"""
sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/cold_room/post_traitement/piv_grenoble/fracture/python_functions/')

from spatial_scale import *

array_alphas = np.zeros(len(ypix))

for i in range(len(ypix)):
    array_alphas[i] = compute_angle(ypix[i],ypix_surf)
array_alphas = np.reshape(array_alphas,(len(array_alphas),1))
ArrAlph = np.tile(array_alphas, (v.shape[0],1,v.shape[2]))
#print(ArrAlph[0,0,:])
#print(ArrAlph[0,:,0])
v_angle_corrected = v_converted_meters/np.cos(ArrAlph)

"""
v_angle_corrected = (1/np.cos(np.radians(alpha0_deg))) * v_converted_meters


kappa_c_vals = []

yindices = np.array([3,4,5,6])

ind_inf_fit = 20
ind_sup_fit = 42

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
    plt.plot(xvals_fit,fit_params[0]*xvals_fit**2+fit_params[1]*xvals_fit+fit_params[2],label='fit',color='red',alpha=0.5)
    plt.plot()
    plt.legend()
#    plt.show()

    kappa_c = 2*fit_params[0]
    print('kappac',kappa_c)
    kappa_c_vals.append(kappa_c)

plt.show()

print(np.mean(np.array(kappa_c_vals)))
print(np.std(np.array(kappa_c_vals)))



# %% enregistrement des données
"""
file_dict_results = 'R:/Gre25/Summary/fracture_postprocessing/resultats/fracture_results.pkl'

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
"""
# %%
