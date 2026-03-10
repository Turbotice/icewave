#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
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

import sys

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/cold_room/post_traitement/piv_grenoble/fracture/python_functions/')

from spatial_scale import *

#matplotlib.use('TkAgg')
# %%

def mat_to_dict(mat_object,ref_matobj):
    """
    Recursively convert a MATLAB structure (HDF5 group or dataset) to a Python dictionary.
    

    INPUTS : - mat_object : matlab object extracted from a .mat file using h5py
             - ref_matobj : matlabo object of reference, main root of the matlab structure, used to dereference some values, 
                     needed for cell_arrays for instance 
                     
    OUTPUT : - whatever was in the .mat file : structure, substructures, cell_array, strings etc..
    
    """
    if isinstance(mat_object, h5py.Dataset):  # If it's a dataset, return its value
        data = mat_object[()]
                  
        # Handle MATLAB strings (stored as bytes)
        if data.dtype == 'uint16':  # Check if it's a string
            # Convert uint16 array to Python string (decode as Unicode characters)
            return ''.join(chr(code_point[0]) for code_point in data)

        # Handle case of cell array
        if data.dtype == 'O':
            new_data = np.empty(data.shape,dtype = object).ravel()
            for i,pointer in enumerate(data.flat):
                new_mat_object = ref_matobj[pointer]
                new_data[i] = mat_to_dict(new_mat_object,ref_matobj)
                
            new_data = new_data.reshape(data.shape)
            return new_data
    
    
        if isinstance(data, np.ndarray):
            data = np.squeeze(data) 
            if data.size == 1:  # If the array contains only one element, convert to float
                return float(data)
            else :
                return data
            
    
    elif isinstance(mat_object, h5py.Group):  # If it's a group (structure), create a dictionary
        result_dict = {}
        for key, item in mat_object.items():
            result_dict[key] = mat_to_dict(item,ref_matobj)  # Recursively call for each element
        return result_dict
    
    else:
        raise TypeError(f"Unsupported type {type(mat_object)}")

#%matplotlib widget
%matplotlib qt

# %% load data
W = 64
Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 0
N = 0
#refimg = 0

date = '0513'
name_frac_file = 'img_seq5'
#camera_SN = '22458101'

f_exc = 0.85
freq_acq = 20

frame_frac = 3209
#computer = 'adour'
ypix_surf = 486
alpha0_deg = 15 # pas mesuré exactement

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/{date}/cameras/frac/img_seq5/'


#path2data = general_folder
#matfile = f'{path2data}PIV_processed_3passages_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

path2data = f'{general_folder}/matData/'
matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'




#matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

with h5py.File(matfile, 'r') as fmat:
    mat_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    mat_dict = mat_to_dict(fmat['m'],fmat['m'])





Vx = mat_dict['Vx']
Vy = mat_dict['Vy']

u = -Vx * (1/freq_acq)
v = -Vy * (1/freq_acq)


xpix = mat_dict['xpix']
ypix = mat_dict['ypix']



#%%

yind = 9 # il vaut mieux utiliser les coordonnées en pixels

frame_plot = frame_frac - 8

plt.figure()
plt.title('vertical displacement map (px/frame)')
plt.imshow(v[frame_plot-i0],vmin=-2.5,vmax=2.5)
plt.colorbar()
plt.show()


plt.figure()
plt.title('displacement profiles (px/frame). frame frac = '+str(frame_frac))
for i in range(frame_plot-i0,frame_plot+10-i0): # fenetre temporelle où il y a la fracture
    plt.plot(v[i,yind,:],label=str(i+i0))
plt.legend()
plt.show()


# %% amplitude max de la vitesse juste avant que ça casse et conversion en amplitude en metres
frame_Vy_max = 3204

xvals = np.arange(v.shape[2]) * 0.075 * W/2 * 1e-2
v_mpersec = v * 0.075 * 1e-2 * freq_acq


ind_inf_fit = 17
ind_sup_fit = 35
#print(xvals)
fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_mpersec[frame_Vy_max-i0,yind,ind_inf_fit:ind_sup_fit],2)

# on peut simplement prendre le max :
max_val = np.max(v_mpersec[frame_Vy_max-i0,yind,ind_inf_fit:ind_sup_fit])
print("amplitude max (m/s) : ", max_val)

print("amplitude convertie en élévation (m) : ", max_val/(2*np.pi*f_exc))
plt.figure()
plt.plot(xvals,v_mpersec[frame_Vy_max-i0,yind,:],label='frame '+str(frame_Vy_max))
plt.show()
"""
# on peut fitter proche du max avec un polydeg2
fit_params


plt.figure()
plt.title('velocity profile (px) at yind='+str(yind))
plt.plot(xvals,v_mpersec[frame_Vy_max-i0,yind,:],label='frame '+str(frame_Vy_max))
plt.plot(xvals[ind_inf_fit:ind_sup_fit],(fit_params[0]*xvals**2 + fit_params[1]*xvals + fit_params[2])[ind_inf_fit:ind_sup_fit],label='fit : $\kappa$ = '+str(np.round(2*fit_params[0],3)))
plt.legend()
plt.show()
"""

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
#%% calcul de tab_dcmsurdpx pour toutes les positions de cases piv (x,y) (avec x et y en pixels)

if len(xpix)!=v.shape[2]:
    xpix = xpix[:-(len(xpix)-v.shape[2])]
else:
    pass

if len(ypix)!=v.shape[2]:
    ypix = ypix[:-(len(ypix)-v.shape[1])]
else:
    pass


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



xvals_px = np.arange(v.shape[2]) * W/2
xvals_test =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
v_converted_meters = v_cm * 1e-2


######### si  on s'intéresse au xhamp de vitesses ####################
v_converted_meters = v_converted_meters * freq_acq # converti en m/s 
######################################################################

#v_converted_meters_1 = v * 0.07 * 1e-2

#%%
plt.figure()
#plt.plot(xvals,v[frame_frac-i0,yind,:]*0.075,label='frame '+str(frame_frac))
plt.plot(xvals_test,v_cm[frame_Vy_max-i0,yind,:],label='frame '+str(frame_frac))

plt.show()

# %% affichage de differents profils avec correction angulaire vs y 
# et variation echelle horizontale vs y
"""
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



kappa_c_velocity_vals = []

yindices = np.array([9,10,11])

ind_inf_fit = 20
ind_sup_fit = 30

plt.figure()
for yind in yindices:
    #print(xvals)
    xvals =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
    xvals_fit = xvals[ind_inf_fit:ind_sup_fit]
    fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_angle_corrected[frame_Vy_max-i0,yind,ind_inf_fit:ind_sup_fit],2)
    #print(fit_params)
#    fit_params_2 = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
#    plt.figure()
    plt.title('profiles at frame '+str(frame_Vy_max))
#    plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals,v_angle_corrected[frame_Vy_max-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals_fit,fit_params[0]*xvals_fit**2+fit_params[1]*xvals_fit+fit_params[2],label='fit',color='red',alpha=0.5)
    plt.plot()
    plt.legend()
#    plt.show()

    kappa_c_velocity = 2*fit_params[0]
    print('max curvature of velocity :',kappa_c_velocity)
    kappa_c_velocity_vals.append(kappa_c_velocity)

plt.show()
kappa_c_velocity_vals = np.array(kappa_c_velocity_vals)


kappa_c_vals = kappa_c_velocity_vals/(2*np.pi*f_exc)


print("converti en courbure spatiale :")
print("kappa_c=",kappa_c_vals)


# %% enregistrement des données
file_dict_results = 'R:/Gre25/Summary/fracture_postprocessing/resultats/fracture_results.pkl'

if os.path.exists(file_dict_results):
    with open(file_dict_results, 'rb') as f:
        dict_results = pickle.load(f)
else:
    dict_results = {}



if date in dict_results:
    dict_results[date][name_frac_file] = {'kappa_c_velocity_vals':kappa_c_velocity_vals,'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf,'f_exc':f_exc}
else:
    dict_results[date] = {name_frac_file: {'kappa_c_velocity_vals':kappa_c_velocity_vals,'kappa_c_vals':kappa_c_vals, 'yindices':yindices , 'ypix':ypix[yindices], 'ypix_surf':ypix_surf,'f_exc':f_exc}}








ee = input('Are you sure you want to erase precedent results ?(y/n)')
if ee=='y':
    with open(file_dict_results, 'wb') as handle:
        pickle.dump(dict_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
else:
    pass

# %%
