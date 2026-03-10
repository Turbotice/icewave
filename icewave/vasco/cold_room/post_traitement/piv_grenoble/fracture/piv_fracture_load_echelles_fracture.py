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
W = 32
Dt = 1
i0 = 0
N = 0

date = '0506'
acq_num = 4
camera_SN = '22458101'

f_exc = 0.94
freq_acq = 20

computer = 'adour'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'

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
ymin = 224
ymax = 470 # bords de l'eau sur images

t_plot = 408 - i0 # numero de frame par rapport à la premiere frame considerée dans la piv (i0)

Vy = mat_dict['Vy']
Vx = mat_dict['Vx']
xpix = mat_dict['xpix']
ypix = mat_dict['ypix']

y_indices = np.where((ypix>ymin)&(ypix<ymax))[0]

Vy_converted_px = Vy/freq_acq

plt.figure()
plt.imshow(Vy_converted_px[t_plot,y_indices,:],extent=[np.min(xpix),np.max(xpix),np.max(ypix[y_indices]),np.min(ypix[y_indices])],vmin=-4,vmax=4)
plt.show()

plt.figure()
plt.imshow(Vy[t_plot,y_indices,:])
plt.show()

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
if computer=='Leyre':
    file_echelles_fracture = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/echelles/echelles_fracture.txt'
else:
    file_echelles_fracture = '/media/turbots/GreDisk/Gre25/Data/0507/cameras/ref_matin/echelles.txt'
data_ech_frac = np.loadtxt(file_echelles_fracture,skiprows=1,usecols=range(5))

d = {}
d['xmoy'] = data_ech_frac[:,1]
d['ymoy'] = data_ech_frac[:,2]
d['dcm'] = data_ech_frac[:,3]
d['dpx'] = data_ech_frac[:,4]
d['delta_y'] = data_ech_frac[:,5]

d['tab_ymoy_refmanip'] = d['ymoy'] - d['delta_y']


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

def compute_aspect_ratio(x,y,d=d):
    popt,pcov = linfitinterp(x,d=d)
    dcm_sur_dpx = popt[0]*y+popt[1]
    return dcm_sur_dpx    


# %% afficher en un point l'évolution au cours du temps

plt.figure()
plt.title('vitesse verticale sur images en px/sec')
for i in [7,8,9,10,11]:
    plt.plot(Vy[300:450,i,30])
plt.show()

alpha = 0.30 # inclinaison camera en radians
Vz = Vy * 1/np.cos(alpha)

plt.figure()
plt.title('Vz en m/sec')
for i in [7,8,9,10,11]:
    plt.plot(mat_dict['t'][300:450],1e-2*Vz[300:450,i,38] * compute_aspect_ratio(xpix[38],ypix[i]))
#plt.xlim(80,100)
plt.show()

plt.figure()
plt.title('(Vz/omega) en cm')
for i in [7,8,9,10,11]:
    plt.plot(Vz[300:450,i,38] * (compute_aspect_ratio(xpix[38],ypix[i])*2*np.pi*f_exc))
plt.grid()
#plt.ylim(10,20)
plt.show()



plt.figure()
plt.title('vitesse verticale en px/sec pour differentes position en x')
for j in [30,37,50]:
    plt.plot(Vz[0:450,9,j])
plt.ylim(-70,70)
plt.show()

plt.figure()
plt.title('vitesse verticale en px/sec pour differentes position en x')
for j in [30,37,50]:
    plt.plot(Vy[0:450,9,j])
plt.ylim(-70,70)
plt.show()

plt.figure()
plt.title('norme vitesse * signe(Vy) en px/sec pour differentes position en x')
for j in [30,37,50]:
    vx = Vx[0:450,9,j]
    vy = Vy[0:450,9,j]
    v = np.sqrt(vx**2+vy**2)
    plt.plot(v*np.sign(vy))
plt.ylim(-70,70)
plt.show()




# %% obtenir une meilleur estimation de la longueur d'onde
tinit2plot = 200

if len(ypix)!=len(Vy[0,:,0]):
    xpix = xpix[1:-1]
    ypix = ypix[1:-1]

array_aspect_ratios = np.zeros_like(xpix)
for i in range(len(array_aspect_ratios)):
    array_aspect_ratios[i] = compute_aspect_ratio(xpix[i],ypix[10])


for i in range(tinit2plot,tinit2plot+30,1):
    plt.plot(np.arange(len(Vz[i,10,:]))*(W/2) * array_aspect_ratios,Vz[i,10,:] * array_aspect_ratios)
#plt.xlim(42,120)
plt.xlabel('x (cm)')
plt.ylabel('Vz (cm/sec)')
plt.show()


for i in range(tinit2plot,tinit2plot+30,1):

    vx = Vx[i,10,:]
    vy = Vy[i,10,:]
    v = np.sqrt(vx**2+vy**2)
    plt.plot(np.arange(len(v))*(W/2) * compute_aspect_ratio(xpix[10],ypix[10]),v*np.sign(vy))
#plt.xlim(42,120)
plt.xlabel('x (cm)')
plt.ylabel('V * sign(Vy) (px/sec)')
plt.title('test en affichant la norme * signe(Vy)')
plt.show()
# %%
