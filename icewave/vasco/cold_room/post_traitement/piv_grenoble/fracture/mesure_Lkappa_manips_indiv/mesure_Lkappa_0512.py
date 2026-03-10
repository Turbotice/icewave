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
import sys
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import find_peaks
#matplotlib.use('TkAgg')
# %%
#%matplotlib widget
%matplotlib qt
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

def matcell2dict_PIV(matcell,dim_keys = 0):
    """ Create a dictionnary for a 2xN matlab cell array 
    
    INPUT : - matcell, cell array converted as an array using h5py and function mat_to_dict, 
                row 0 -> keys and row 1 -> values 
                
    OUTPUT : - a python dictionnary whose keys correspond to keys stored in the first dimension 
    
    """
    my_dict = {}
    keys = matcell[dim_keys]
    for i,key in enumerate(keys) :
        my_dict[key] = matcell[1,i]
        
    return my_dict




# %% load data
W = 64
Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 0
N = 0

date = '0512'
name_frac_file = 'img_seq2'
#camera_SN = '22458101'

f_exc = 0.89
freq_acq = 20

frame_frac = 2572
#computer = 'adour'
ypix_surf = 490

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/{date}/cameras/frac/{name_frac_file}/'


path2data = f'{general_folder}/{f_exc}Hz_{freq_acq}Hz/matData/'
matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'


#general_folder = f'R:/Gre25/{date}/cameras/frac/img_seq_3/'

#general_folder = f'R:/Gre25/{date}/cameras/frac/image_sequence/'

path2data = general_folder + 'matData/'
matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'

with h5py.File(matfile, 'r') as fmat:
    mat_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    mat_dict = mat_to_dict(fmat['m'],fmat['m'])
# %% visualize piv data
ymin = 100
ymax = 400 # bords de l'eau sur images

t_plot = 408 - i0 # numero de frame par rapport à la premiere frame considerée dans la piv (i0)

Vy = mat_dict['Vy']
Vx = mat_dict['Vx']
xpix = mat_dict['xpix']
ypix = mat_dict['ypix']




Vx = mat_dict['Vx']
Vy = mat_dict['Vy']

u = -Vx * (1/freq_acq)
v = -Vy * (1/freq_acq)


xpix = mat_dict['xpix']
ypix = mat_dict['ypix']


y_indices = np.where((ypix>ymin)&(ypix<ymax))[0]

Vx_converted_px = Vx/freq_acq
Vy_converted_px = Vy/freq_acq



plt.figure()
plt.imshow(Vy_converted_px[t_plot,y_indices,:],extent=[np.min(xpix),np.max(xpix),np.max(ypix[y_indices]),np.min(ypix[y_indices])],vmin=-np.max(Vy_converted_px)/5,vmax=np.max(Vy_converted_px)/5)
plt.show()

plt.figure()
plt.imshow(Vy[t_plot,y_indices,:])
plt.show()

#%%

yind = 6
plt.figure()
plt.title('profile at yind='+str(yind)+' when fracture occurs')
for i in range(frame_frac - 10, frame_frac+2):
    plt.plot(v[i,yind,:],label=str(i))
plt.legend()
plt.show()




#%%load file echelles fracture pour avoir en x et y les variations de dpx/dcm
#if computer=='Leyre':
#    file_echelles_fracture = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/echelles/echelles_fracture.txt'
if system_loc=='windows_server':
    file_echelles_fracture = 'R:/Gre25/Data/0507/cameras/ref_matin/echelles.txt'
elif system_loc=='linux_server':
    file_echelles_fracture = '/media/turbots/GreDisk/Gre25/Data/0507/cameras/ref_matin/echelles.txt'
data_ech_frac = np.loadtxt(file_echelles_fracture,skiprows=1,usecols=range(6))

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
#d['delta_y'] = data_ech_frac[:,5]

d['tab_ymoy_refmanip'] = d['ymoy']# - d['delta_y']



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
v_converted_meters = v_cm * 1e-2

v_converted_meters = v_converted_meters * freq_acq # meters per sec

xvals_px = np.arange(v.shape[2]) * W/2

xvals_test =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
#v_converted_meters_1 = v * 0.07 * 1e-2

#%%
plt.figure()
#plt.plot(xvals,v[frame_frac-i0,yind,:]*0.075,label='frame '+str(frame_frac))
plt.plot(xvals_test,v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:],label='frame '+str(frame_frac))
peaks = find_peaks(-(v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:]))[0]
print(peaks)
plt.vlines(xvals_test[peaks[0]],np.nanmin((v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:])),np.nanmax((v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:])),linestyles='--',colors='tab:orange')
plt.vlines(xvals_test[peaks[1]],np.nanmin((v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:])),np.nanmax((v_cm[frame_frac-i0,yind,:]-v_cm[frame_frac-1-i0,yind,:])),linestyles='--',colors='tab:orange')
Lkappa_test = xvals_test[peaks[1]] - xvals_test[peaks[0]]
plt.title('$L_{\kappa}$ = '+str(Lkappa_test)+' m')
plt.xlabel('x (m)')
plt.ylabel('Vz($t_{frac}$) - Vz($t_{frac}$-$\Delta$t)')
plt.show()

# %% affichage de differents profils avec correction angulaire vs y 
# et variation echelle horizontale vs y

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

kappa_c_velocity_vals = []

yindices = np.array([6,7,8])

ind_inf_fit = 24
ind_sup_fit = 35

plt.figure()
for yind in yindices:
    #print(xvals)
    xvals =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
    xvals_fit = xvals[ind_inf_fit:ind_sup_fit]
    #fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_angle_corrected[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit]-v_angle_corrected[frame_frac-1-i0,yind,ind_inf_fit:ind_sup_fit],2)

    plt.title('profiles at frame '+str(frame_frac))
#    plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals,v_angle_corrected[frame_frac-i0,yind,:]-v_angle_corrected[frame_frac-1-i0,yind,:],label='yind='+str(yind))
    #plt.plot(xvals_fit,fit_params[0]*xvals_fit**2+fit_params[1]*xvals_fit+fit_params[2],label='fit',color='red',alpha=0.5)
    plt.plot()
    plt.legend()
#    plt.show()


plt.show()




# %%
