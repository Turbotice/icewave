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
#Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 24
N = 0
refimg = 24

date = '0506'
name_frac_file = 'image_sequence'

#camera_SN = '22458101'

f_exc = 0.9597
freq_acq = 20

frame_frac = 698

ypix_surf = 400
alpha0_deg = 19.5 # angle caméra en degrés
#computer = 'adour'

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/PIV_results/{date}_frac_image_sequence/'


path2data = general_folder

matfile = f'{path2data}PIV_processed_3passages_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

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




# %% visualize piv data
yind = 4 # il vaut mieux utiliser les coordonnées en pixels

plt.figure()
plt.title('elevation map (px) at frame 10')
plt.imshow(v[10],vmin=-3,vmax=3)
plt.colorbar()
plt.show()


plt.figure()
plt.title('profile at yind='+str(yind)+' for the first frames')
for i in range(40):
    plt.plot(v[i,yind,:],label=str(i))
plt.legend()
plt.show()




frame_plot = frame_frac - 5

plt.figure()
plt.title('elevation map (px) at frame '+str(frame_plot))
plt.imshow(v[frame_plot-i0],vmin=-5,vmax=5)
plt.colorbar()
plt.show()


plt.figure()
plt.title('profile : '+str(yind)+' , frame frac = '+str(frame_frac))
for i in range(frame_plot-i0,frame_plot+10-i0): # fenetre temporelle où il y a la fracture
    plt.plot(v[i,yind,:],label=str(i+i0))
plt.legend()
plt.show()


#%% partie pour trouver image de reference

xind_min = 11
xind_max = 48
tind_min = 0
tind_max = 20


tind_profile_max = np.argmax(np.sum(v[tind_min:tind_max,yind,xind_min:xind_max],axis=1))
tind_profile_min = np.argmin(np.sum(v[tind_min:tind_max,yind,xind_min:xind_max],axis=1))

# on en déduit, comme c'est périodique :
tind_profile_centre = np.round((tind_profile_min+tind_profile_max)/2).astype(int)
u_field_centre = u[tind_profile_centre,:,:]
v_field_centre = v[tind_profile_centre,:,:]

plt.figure()

plt.plot(v[tind_profile_min,yind,:],label=str(tind_profile_min))
plt.plot(v[tind_profile_max,yind,:],label=str(tind_profile_max))
plt.plot(v[tind_profile_centre,yind,:],label=str(tind_profile_centre))

plt.legend()
plt.show()


u_shifted = u - np.tile(u_field_centre,(u.shape[0],1,1))
v_shifted = v - np.tile(v_field_centre,(v.shape[0],1,1))

plt.figure()
plt.title('profile (shifted) at yind='+str(yind)+' for the first frames')
for i in range(100,140):
    plt.plot(v_shifted[i,yind,:],label=str(i))
plt.legend()
plt.show()

# à faire dans un second temps : utiliser v_shifted à la place de v et moyenner les les v_field_centre pour ameliorer correction

# %% courbure de la plaque juste avant la fracture ?

xvals = np.arange(v.shape[2]) * 0.075 * W/2 * 1e-2
v_converted_meters = v * 0.075 * 1e-2


# 0.075 est la valeur approximative de dcm/dpx, mais une utilisation de la variation de dcm/dpx 
# est requise pour une meilleure estimation de kappa_c


ind_inf_fit = 20
ind_sup_fit = 42
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
"""
tab_dcm = np.array([23,30,30])
tab_dpx = np.array([259,384,497])

delta_y = 1100-200 # car photos regles pas memes dimensions que photos fracture...
tab_y = np.array([1100,1185,1350]) - delta_y
plt.figure()
plt.plot(tab_y,tab_dcm/tab_dpx,'o')
popt,pcov = curve_fit(lambda x,a,b:a*x + b,tab_y,tab_dcm/tab_dpx)
plt.plot(tab_y,popt[0]*tab_y+popt[1])
plt.show()

def compute_aspect_ratio(y,popt=popt):
    return popt[0]*y+popt[1]
"""
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
d['ypix_surf'] = data_ech_frac[:,4]
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

plt.figure()
plt.plot(xvals,v[frame_frac-i0,yind,:]*0.075,label='frame '+str(frame_frac))
plt.plot(xvals,v_cm[frame_frac-i0,yind,:],label='frame '+str(frame_frac))

plt.show()


xvals_px = np.arange(v.shape[2]) * W/2

xvals_test = np.arange(v.shape[2]) * 0.075 * W/2 * 1e-2

#%% hidden cell
'''
for yind in [3,4,5,6,7]:
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
    print(kappa_c)
'''
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



ind_inf_fit = 25
ind_sup_fit = 37


kappa_c_vals = []

yindices = np.array([2,3,4,5,6])

plt.figure()
for yind in yindices:
    #print(xvals)
    xvals =  xvals_px * DCM_SUR_DPX_2[yind,:] * 1e-2
    xvals_fit = xvals[ind_inf_fit:ind_sup_fit]
    fit_params = np.polyfit(xvals_fit,v_angle_corrected[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
    #print(fit_params)
#    fit_params_2 = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
#    plt.figure()
    plt.title('profiles at frame '+str(frame_frac))
#    plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals,v_angle_corrected[frame_frac-i0,yind,:],label='yind='+str(yind))
    plt.plot(xvals_fit,fit_params[0]*xvals_fit**2+fit_params[1]*xvals_fit+fit_params[2],label='fit',color='red')
    plt.plot()
    plt.legend()
#    plt.show()

    kappa_c = 2*fit_params[0]
    print('kappac',kappa_c)
    kappa_c_vals.append(kappa_c)

plt.show()

print(np.mean(np.array(kappa_c_vals)))
print(np.std(np.array(kappa_c_vals)))


# %%
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

# %%
