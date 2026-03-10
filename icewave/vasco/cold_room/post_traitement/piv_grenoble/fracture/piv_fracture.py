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
i0 = 350
N = 1400

#i0 = 800
#N = 1400

date = '20241129'
acq_num = 7
camera_SN = '22458101'

f_exc = 0.94
freq_acq = 15

computer = 'dell_vasco'

if computer=='dell_vasco':
    general_folder = f'D:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
elif computer=='Leyre':
    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'



path2data = f'{general_folder}/{f_exc}Hz_{freq_acq}Hz/matData/'
matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'

with h5py.File(matfile, 'r') as fmat:
    mat_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    mat_dict = mat_to_dict(fmat['m'],fmat['m'])
# %% visualize piv data
ymin = 224
ymax = 470 # bords de l'eau sur images

t_plot = 760 - i0 # numero de frame par rapport à la premiere frame considerée dans la piv (i0)

Vy = mat_dict['Vy']
Vx = mat_dict['Vx']
xpix = mat_dict['xpix']
ypix = mat_dict['ypix']

y_indices = np.where((ypix>ymin)&(ypix<ymax))[0]

plt.figure()
plt.imshow(Vy[t_plot,y_indices,:],extent=[np.min(xpix),np.max(xpix),np.max(ypix[y_indices]),np.min(ypix[y_indices])])
plt.show()

#%% calcul echelles avec photos regle
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
#%%load file echelles fracture pour avoir en x et y les variations de dpx/dcm
if ordi=='Leyre':
    file_echelles_fracture = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/echelles/echelles_fracture.txt'
elif ordi=='dell_vasco':
    file_echelles_fracture = 'D:/Gre24/Data/20241129/echelles/echelles_fracture.txt'
data_ech_frac = np.loadtxt(file_echelles_fracture,skiprows=1,usecols=range(6))



#%% input params and input values for postprocessing

input_params = csv2dict(csv_file_path)
# on rajoute quelquues params qui dependent des parameters d'entrée :
#input_params['x0'] = input_params['sensor_width']/2 - input_params['delta_x']
#input_params['y0'] = input_params['sensor_height']/2 - input_params['delta_y']
#input_params['focale'] = input_params['focale_mm']*1e-3/input_params['pixel_size'] # distance focale convertie en nombre de pixels du capteur
input_params['alpha_0'] = input_params['alpha_0_deg']*np.pi/180

###############################""
input_values_postprocessing = {}

# pour le 29/11 acq7:

input_values_postprocessing['idx_profile'] = 16

input_values_postprocessing['wavelength_measurement_statio'] = {}
input_values_postprocessing['wavelength_measurement_statio']['framestart'] = 0
input_values_postprocessing['wavelength_measurement_statio']['framestop'] = 300
#input_values_postprocessing['wavelength_measurement_statio']['idx_profile'] = 16
input_values_postprocessing['wavelength_measurement_statio']['xinf'] = 0.3
input_values_postprocessing['wavelength_measurement_statio']['xsup'] = 1.3

input_values_postprocessing['Vz_measurement_pre_fracture'] = {}
#input_values_postprocessing['Vz_measurement_pre_fracture']['idx_profile'] = 16
input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture'] = 760
input_values_postprocessing['Vz_measurement_pre_fracture']['xinf'] = 0.4
input_values_postprocessing['Vz_measurement_pre_fracture']['xsup'] = 1.2



Vz = (-1) * np.transpose(Vy * 1/np.cos(input_params['alpha_0']),(2,1,0)) # unit : px/sec






# %% mesure longueur d'onde pour differents y (devraient être égales)

# autre methode : mesurer pour chaque profile puis moyenner

idx_profile = input_values_postprocessing['idx_profile']

framestart = input_values_postprocessing['wavelength_measurement_statio']['framestart']
framestop = input_values_postprocessing['wavelength_measurement_statio']['framestop']
xinf = input_values_postprocessing['wavelength_measurement_statio']['xinf']
xsup = input_values_postprocessing['wavelength_measurement_statio']['xsup']


xreal_profile = 1e-2 * np.arange(Vz.shape[0])*(W/2)*compute_aspect_ratio(ypix[idx_profile])
vz_profile_mpersec = Vz[:,idx_profile,:] * 1e-2 * compute_aspect_ratio(ypix[idx_profile])

x_indices2plot = np.where((xreal_profile>xinf)&(xreal_profile<xsup))[0]

wavelengths_measured = []
for i in range(framestop-framestart):
    smoothed_profile = gaussian_filter(np.abs(vz_profile_mpersec[x_indices2plot,i]),1.5)
    if np.max(smoothed_profile)>0.010:
        minima_indices , properties= find_peaks(-smoothed_profile, width=7)
        #local_minima = [(Xreal1_centers[idx_profile,i], smoothed_profile[i]) for i in minima_indices]
        plt.plot(xreal_profile[x_indices2plot],np.abs(vz_profile_mpersec[x_indices2plot,i]),'.')
        plt.plot(xreal_profile[x_indices2plot],smoothed_profile)
        print(minima_indices)
        if len(minima_indices)==2:
            for j in minima_indices:
                plt.vlines(xreal_profile[j+x_indices2plot[0]],0,np.max(np.abs(vz_profile_mpersec[x_indices2plot,i])),linestyle='--',color='red')
            dX = xreal_profile[minima_indices[1]+x_indices2plot[1]] - xreal_profile[minima_indices[0]+x_indices2plot[0]]
            wavelength = dX * 2
            print('delta X = ',dX)
            print('wavelength = ',wavelength)
            # faire evenetuellement une  liste qui dit quelles frames ont été utilisées...
            wavelengths_measured.append(wavelength)
    if i%10==0:
        plt.show()
wavelengths_measured = np.array(wavelengths_measured)

average_wavelength = np.mean(wavelengths_measured)
std_wavelength = np.std(wavelengths_measured)
print('average wavelength : '+str(average_wavelength)+' +- '+str(std_wavelength))



# %% Plot the curvature of Vz vs amplitude of Vz to see if if satisfies the linear model

idx_profile = input_values_postprocessing['idx_profile']

frame_fracture = input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture']
xinf = input_values_postprocessing['Vz_measurement_pre_fracture']['xinf']
xsup = input_values_postprocessing['Vz_measurement_pre_fracture']['xsup']
x_indices2plot = np.where((xreal_profile>xinf)&(xreal_profile<xsup))[0]


freq_acq = input_params['freq_acq']
f_exc = input_params['f_exc']
T_exc = 1/f_exc
T_exc_frames = T_exc * float(freq_acq)

#array_indices_maxamplitude = np.arange(1,100,16)#np.array([1,18,18+16,18+16+17,18+16+17+16,18+16+17+16+17,18+16+17+16+17+16])
first_maxamp_idx = 1
shift = 0

nb_profiles = 30
PROFILES = []
colors = cm.viridis(np.linspace(-0.5, 0.8, nb_profiles))
for i in range(first_maxamp_idx,nb_profiles,1):
    profile0 = vz_profile_mpersec[:,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) - 1]
    profile1 = vz_profile_mpersec[:,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int)]
    profile2 = vz_profile_mpersec[:,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) + 1]
    
    profiles = [profile0,profile1,profile2]
    array_sums = np.array([np.sum(profile0[x_indices2plot]),np.sum(profile1[x_indices2plot]),np.sum(profile2[x_indices2plot])])    
    index_choice = np.argmax(array_sums)
    plt.plot(xreal_profile[x_indices2plot],profiles[index_choice][x_indices2plot],color=colors[i])
    shift += index_choice - 1
    PROFILES.append(profiles[index_choice])
#plt.xlim(-0.35,0.35)
#plt.ylim(-0.01,0.046)
plt.show()

PROFILES = np.array(PROFILES)
amplitudes = np.max(PROFILES[:,x_indices2plot],axis=1)
indices_maxes = np.argmax(gaussian_filter(PROFILES[:,x_indices2plot],1.5,axes=(1)),axis=1) + np.min(x_indices2plot)
semilengthfit = 30
semilengthfit_broken = 10
POPT = []
PCOV = []
second_derivatives = []

for i in range(PROFILES.shape[0]):
    profile = PROFILES[i,:]

    second_derivatives.append(np.gradient(np.gradient(profile)))

   # print(indices_maxes[i]-semilengthfit)
   # print(indices_maxes[i]+semilengthfit)
    profile2fit = profile[indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    x2fit = xreal_profile[indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    
    popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    if np.sqrt(np.diag(pcov)[0]) > 0.05:
            profile2fit = profile[indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            x2fit = xreal_profile[indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    POPT.append(popt)
    PCOV.append(pcov)
    plt.plot(xreal_profile[x_indices2plot],PROFILES[i,x_indices2plot],'.',alpha=0.3,color=colors[i])
    plt.plot(x2fit,(popt[0]*0.5)*x2fit**2+popt[1]*x2fit+popt[2],color=colors[i])

plt.xlabel('x (m)')
plt.ylabel('Vz (m/s)')
plt.title('Maxima successifs de Vz sur un profil')
plt.show()

POPT = np.array(POPT)
PCOV = np.array(PCOV)
second_derivatives = np.array(second_derivatives)
curvatures = POPT[:,0].flatten()

plt.title('curvature of Vz   vs amplitude')
plt.plot(amplitudes,np.abs(curvatures),'o',label='data')
wavelength = average_wavelength
wavenumber = 2*np.pi/wavelength
ampvals = np.linspace(0,np.max(amplitudes))
plt.plot(ampvals , wavenumber**2 * ampvals,label='$k_{avg}^{2}\cdot$amplitude')
slope,var_slope = curve_fit(lambda x,a:a*x,amplitudes,np.abs(curvatures))
slope = slope[0]
err_slope = np.sqrt(var_slope[0][0])
plt.plot(ampvals , slope * ampvals,label='fit : k²='+str(np.round(slope,2)))
#plt.xlim(0,0.05)
#plt.ylim(0,0.8)
#plt.loglog()
plt.xlabel('amplitude (m/s)')
plt.ylabel('Vz curvature (m^-1 . s^-1)')
plt.legend()
plt.show()


# %% measure the amplitude (and the curvature) just before fracture
i0 = input_params['i0']
frame_fracture = input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture']
time_idx_fracture = frame_fracture - i0 # - t0 (-t0 si on decale le film)


nb_time_indices = 14
colors = cm.winter(np.linspace(0, 0.8, nb_time_indices))
maxtemp=0

extremum = 'max' # preciser si la derniere amplitude maximale pour la vitesse correspond à un min ou un max
for i in range(nb_time_indices):
    plt.plot(xreal_profile[x_indices2plot],vz_profile_mpersec[x_indices2plot,i+time_idx_fracture],color=colors[i])#,label=str(time_idx_fracture+i0+t0+i))
    if extremum=='max':
        max_i = np.max(vz_profile_mpersec[x_indices2plot,i+time_idx_fracture][10:-20])
    elif extremum=='min':
        max_i = np.min(vz_profile_mpersec[x_indices2plot,i+time_idx_fracture][10:-20])
#    plt.plot(Xreal1_centers[idx_profile,x_indices2plot],gaussian_filter(Vz[x_indices2plot,idx_profile,i+time_idx_fracture],1.5),'--',color='gray')

    if max_i>maxtemp:
        maxtemp = max_i
plt.title('idx y profile = '+str(idx_profile))
print(maxtemp,' m/s')  
#plt.legend()
plt.show()

# dans l'approximation linéaire des ondes, la mesure de l'amplitude suffit pour déduire la courbure

print('correspond for the linear model to : amplitude = ',maxtemp / (2*np.pi*f_exc),' meters')


