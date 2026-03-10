#%% -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import glob
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import time 

import h5py
import os

#kerdir = os.getcwd()
#os.chdir('/home/vzanchi/Bureau/Turbotice_git/')
import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp
#os.chdir(kerdir)
#%%
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

#%%
"""
from scipy.interpolate import griddata
from matplotlib.animation import FuncAnimation
import scipy
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

"""
#from mat73 import loadmat
# pour obtenir quelque chose en unité de longueur on fait : delta_z = Vz*Dt/fps


#%% tests manip frigo vasco

# noms des dossiers contenant les données

def fonction_test(input_params = dict):
    for key in input_params:
        globals()[key] = input_params[key]
    #print(locals().keys())
    #print(date)
    base = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    path2data = f'{base}matData/'

    matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'
    global mat_dict
    with h5py.File(matfile, 'r') as fmat:
        mat_dict = {}
        
        print('Top-level keys : ', list(fmat.keys()))

        mat_dict = mat_to_dict(fmat['m'],fmat['m'])

    print("KEYS : ",mat_dict.keys())
    # Convert Vy displacement field to vertical field Vz 

    Vy = mat_dict['Vy'] # time, y, x
    Vy_transposed = np.transpose(Vy,(2,1,0)) 
    
    if unit_input_velocity == 'pxpersec':
        fps = freq_acq
        t = mat_dict['t']
    elif unit_input_velocity == 'pxperframe':
        t = np.arange(len(mat_dict['t']))/fps  # temps en secondes
    
    if box_position_information==False:
        x = np.arange(len(mat_dict['x']))*(W/2) + W + 0.5 # pixel position of each PIV box
        y = np.arange(len(mat_dict['y']))*(W/2) + W + 0.5 # pixel position of each PIV box 
    elif box_position_information==True:
        x = mat_dict['xpix']
        y = mat_dict['ypix']
        if a!=0:
            x = x[a:-a]
            y = y[a:-a]




    # compute interpolation of pixel displacement vertical velocity field 
    Fy = RegularGridInterpolator((x,y,t), Vy_transposed)

    # compute vertical velocity Vz
    if unit_input_velocity=='pxpersec':
        Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed/fps,x,y,t,y0,h_camera,
                                                alpha_0,focale,fps,Dt)
    elif unit_input_velocity=='pxperframe':
        Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed,x,y,t,y0,h_camera,
                                                alpha_0,focale,fps,Dt)

    # Plot a plt.colormesh 
    #W = mat_dict['PIV_param']['w']
    x_edges = np.resize(x,x.size + 1)
    y_edges = np.resize(y,y.size + 1)
    x_edges[-1] = x_edges[-2] + W
    y_edges[-1] = y_edges[-2] + W
    x_edges = x_edges - W/2
    y_edges = y_edges - W/2


    Xedges,Yedges = np.meshgrid(x_edges,y_edges, indexing = 'ij')

    Xreal,Yreal = dp.projection_real_space(Xedges, Yedges, x0, y0, h_camera
                                        ,alpha_0,focale)

    fig, ax = plt.subplots()
    pc = ax.pcolormesh(Xreal,Yreal,Vz[:,:,20],shading = 'auto',vmin=-0.01,vmax=0.01) # pcolormesh object !
    fig.colorbar(pc)
    ax.set_ylim(-0.7,-0.2)
    ax.set_xlim(-1.0,1.0)

    X,Y = np.meshgrid(x,y)

    Xreal_centers,Yreal_centers = dp.projection_real_space(X, Y, x0, y0, h_camera, alpha_0, focale)

    return mat_dict,t,Xreal,Yreal,Vz,pc,Xreal_centers,Yreal_centers

###############################


input_params = {}
input_params['date'] = '20241129'
input_params['acq_num'] = 7
input_params['f_exc'] = 0.94
input_params['freq_acq'] = 15
input_params['cam_SN'] = '22458101'
input_params['delta_x'] = 0
input_params['delta_y'] = 900
input_params['imwidth'] = 2048
input_params['imheight'] = 710
input_params['x0'] = 1024 - input_params['delta_x']
input_params['y0'] = 1024 - input_params['delta_y']
input_params['pixel_size'] = 5.5e-6 # taille d'un pixel sur le capteur en metres
input_params['focale_mm'] = 8.31 # distance focale en mm (mesurée par calibration du systeme optique dans dans une gamme de distance similaire à l'expérience)
input_params['focale'] = input_params['focale_mm']*1e-3/input_params['pixel_size'] # distance focale convertie en nombre de pixels du capteur
input_params['h_camera'] = 0.45 # hauteur de la camera en metres par rapport à la surface de l'eau
input_params['alpha_0_deg'] = 17.1
input_params['alpha_0'] = input_params['alpha_0_deg']*np.pi/180
input_params['unit_input_velocity'] = 'pxpersec'
input_params['box_position_information'] = True
input_params['a'] = 0 # number of PIV boxes croped on each side
input_params['W'] = 32
input_params['Dt'] = 1
input_params['i0'] = 350
input_params['N'] = 1400
###############################

mat_dict,t,Xreal1,Yreal1,Vz1,pc1,Xreal1_centers,Yreal1_centers = fonction_test(input_params=input_params)
#Xreal2,Yreal2,Vz2,pc2,Xreal2_centers,Yreal2_centers = fonction_test(cam_SN = '40300722')

#delta_x_cam = 0.184 # distance en metres entre les deux cameras
#Xreal2_shift = Xreal2 + delta_x_cam

#plt.figure()
#plt.pcolormesh(Xreal1,Yreal1,Vz1[:,:,20],shading = 'auto')
#plt.pcolormesh(Xreal2_shift,Yreal2,Vz2[:,:,20],shading = 'auto',vmin=-0.003,vmax=0.003)
#plt.colorbar()
#plt.ylim(-0.15,0.1)
#plt.xlim(-0.17,0.32)
#plt.show()

#%% creer une video
create_video = False
if create_video==True:
    t_video = np.arange(Vz1.shape[2])

    def update(frame):
        plt.cla()  # Efface le contenu de l'axe à chaque trame
    #    plt.imshow(matrix[frame,:,:])
        plt.pcolormesh(Xreal1,Yreal1,Vz1[:,:,frame],shading = 'auto',vmin=-0.1,vmax=0.1)
    #    plt.pcolormesh(Xreal2_shift,Yreal2,Vz2[:,:,frame],shading = 'auto',vmin=-0.003,vmax=0.003)
        plt.xlabel('x (m)',fontsize=15)
        plt.ylabel('y (m)',fontsize=15)
    #    plt.ylim(-0.15,0.1)
    #    plt.xlim(-0.17,0.32)
    # Crée la figure et l'axe pour l'animation
    fig, ax = plt.subplots() 
    frames_to_include = range(0, len(t_video), 1)
    ani = FuncAnimation(fig, update, frames=frames_to_include, repeat=False,blit=False)

    output_file = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/manip_fracture/Acquisition_7/camera_22458101/0.94Hz_15Hz/video_PIV1.mp4'
    # Save the animation as a video (specify the writer and additional options)
    ani.save(output_file, writer='ffmpeg', fps=15,dpi=100)  # Adjust the FPS as needed
    plt.show()





# %% mesure longueur d'onde pour differents y (devraient être égales)
framestart = 0 # le 0 correspond à la frame numero i0
framestop = 300

#idx_profile = 8
idx_profile = 16

"""
for i in range(framestop-framestart):
    plt.plot(Xreal1_centers[idx_profile,:],np.abs(Vz1[:,idx_profile,i]))

plt.show()

#idxmin_x = 15
idxmin_x = 30
#idxmax_x = 55
idxmax_x = 110

tempavg = np.mean(np.abs(Vz1[idxmin_x:idxmax_x,idx_profile,:])/np.max(np.abs(Vz1[idxmin_x:idxmax_x,idx_profile,:])),axis=1)
minima_indices , properties= find_peaks(-tempavg, width=5)

# Get the x and y values of the minima
local_minima = [(Xreal1_centers[idx_profile,i], tempavg[i]) for i in minima_indices]

for i in minima_indices:
    plt.vlines(Xreal1_centers[idx_profile,i+idxmin_x],0,np.max(tempavg),linestyle='--',color='red')
plt.plot(Xreal1_centers[idx_profile,idxmin_x:idxmax_x],tempavg)
plt.show()

dX = Xreal1_centers[idx_profile,minima_indices[1]+idxmin_x] - Xreal1_centers[idx_profile,minima_indices[0]+idxmin_x]
print('delta X = ',dX)
print('lambda = ',dX*2)
wavelength = dX * 2
"""
# autre methode : mesurer pour chaque profile puis moyenner

wavelengths_measured = []
x_indices2plot = np.where((Xreal1_centers[idx_profile,:]>-0.45)&(Xreal1_centers[idx_profile,:]<0.50))[0]
for i in range(framestop-framestart):
    smoothed_profile = gaussian_filter(np.abs(Vz1[x_indices2plot,idx_profile,i]),1.5)
    if np.max(smoothed_profile)>0.010:
        minima_indices , properties= find_peaks(-smoothed_profile, width=7)
        #local_minima = [(Xreal1_centers[idx_profile,i], smoothed_profile[i]) for i in minima_indices]
        plt.plot(Xreal1_centers[idx_profile,x_indices2plot],np.abs(Vz1[x_indices2plot,idx_profile,i]),'.')
        plt.plot(Xreal1_centers[idx_profile,x_indices2plot],smoothed_profile)
        print(minima_indices)
        if len(minima_indices)==2:
            for j in minima_indices:
                plt.vlines(Xreal1_centers[idx_profile,j+x_indices2plot[0]],0,np.max(np.abs(Vz1[x_indices2plot,idx_profile,i])),linestyle='--',color='red')
            dX = Xreal1_centers[idx_profile,minima_indices[1]+idxmin_x] - Xreal1_centers[idx_profile,minima_indices[0]+idxmin_x]
            wavelength = dX * 2
            print('delta X = ',dX)
            print('wavelength = ',wavelength)
            # faire evenetuellement une  liste qui dit quelles frames ont été utilisées...
            wavelengths_measured.append(wavelength)
    if i%10==0:
        plt.show()
wavelengths_measured = np.array(wavelengths_measured)

average_wavelength = np.mean(wavelengths_measured)
print('average wavelength : ',average_wavelength)
# faire subpix pour améliorer
# %% choix de la frame de reference (à laquelle le profil est nul)

#idx_x = 30
idx_x = 60
#idx_y = 8
idx_y = 16

# l'idéal c'est qu'au début de la vidéo la surface soit immobile...

for i in range(15,21):
    plt.plot(Vz1[:,idx_y,i],label='i = '+str(i))
plt.legend()
plt.show()

for i in range(325,333):
    plt.plot(Vz1[:,idx_y,i],label='i = '+str(i))
plt.legend()
plt.show()


# choix de t0 pour l'integration : on plot Vz vs t en diferents points à partir du début du dataset
# pour le 29/11 par exemple on choisit i = 16 (correspondant à frame i0+i=816)

#t0 = 16

#t0 = 330

t0 = 0
#%%  convertir  Vz(x,y,t) en une élévation z(x,y,t)
N = input_params['N']
i0 = input_params['i0']
Dt = input_params['Dt']
if Vz1.shape[2] != N - i0 - t0 - Dt:
    Vz1 = Vz1[:,:,t0:]
    t = t[t0:]
print('new shape of Vz1 : ',Vz1.shape)
# premiere methode : intégrer
dt = t[1] - t[0]
cumsumVz1 = np.zeros_like(Vz1)

for i in range(Vz1.shape[0]):
    for j in range(Vz1.shape[1]):
        cumsumVz1[i,j,:] = np.cumsum(Vz1[i,j,:])*dt


cumtrapVz1 = np.zeros_like(Vz1)
for i in range(Vz1.shape[0]):
    for j in range(Vz1.shape[1]):
        cumtrapVz1[i,j,:] = cumulative_trapezoid(Vz1[i,j,:],t,initial=0)



plt.title('premiere methode (mieux pour fracture)')
plt.plot(t,Vz1[idx_x,idx_y,:],label='Vz')
#plt.plot(t,cumsumVz1[idx_x,idx_y,:])
plt.plot(t, cumtrapVz1[idx_x,idx_y,:],label='elevation $\zeta$')
plt.xlim(25,40)
########## tentatives pour corriger le shift ###################
'''
from scipy.signal import savgol_filter
filtered_cumsumVz1 = savgol_filter(cumsumVz1[idx_x,idx_y,:],200,4)
plt.plot(t , filtered_cumsumVz1)
plt.plot(t, cumsumVz1[idx_x,idx_y,:] - filtered_cumsumVz1)
from scipy.ndimage import gaussian_filter
filtered_cumsumVz1 = gaussian_filter(cumsumVz1[idx_x,idx_y,:],20)
plt.plot(t , filtered_cumsumVz1)
plt.plot(t, cumsumVz1[idx_x,idx_y,:] - filtered_cumsumVz1)
'''

################################################################

plt.ylim(-0.01,0.01)
#plt.xlim(0,15)
plt.grid()
plt.legend()
##########
plt.show()

"""
f_exc = input_params['f_exc']
T_exc = 1/f_exc

nbperiods2avg = 20
indicestimeavg = np.where((t>=0)&(t<nbperiods2avg*T_exc))[0]
tab_constant_shift_before = np.mean(cumsumVz1[:,:,indicestimeavg],axis=2)

tab_constant_shift = np.tile(tab_constant_shift_before,(1,1,cumsumVz1.shape[2]))
z1 = cumsumVz1 - tab_constant_shift

"""

# deuxieme methode : dephaser de +pi/2 le signal Vz  et diviser par omega
f_exc = input_params['f_exc']
T_exc = 1/f_exc

t_new = t - T_exc/4
z1 = (1/(2*np.pi*f_exc)) * Vz1

plt.title('deuxieme methode (moins bien pour fracture...)')
plt.plot(t,Vz1[idx_x,idx_y,:])
plt.plot(t_new,z1[idx_x,idx_y,:])

#plt.ylim(-0.005,0.005)
#plt.xlim(0,50)
plt.show()

"""FFTVZ1test = np.fft.fft(Vz1[20,9,:])
FFTVZ1test[50:550] = 0
iffttest = np.fft.ifft(FFTVZ1test)
plt.plot(Vz1[20,9,:])
plt.plot(iffttest)
plt.ylim(-0.02,0.02)
plt.show()

plt.plot(np.cumsum(Vz1[20,9,:]))
plt.plot(np.cumsum(iffttest))"""


# choix final : methode 1
z1 = cumtrapVz1

#%% afficher des profils 
idx_time_init = 0
for i in range(20):
    plt.plot(Xreal1_centers[idx_y,:],Vz1[:,idx_y,idx_time_init+i])
plt.xlim(-0.3,0.45)
plt.show()
for i in range(20):
    plt.plot(Xreal1_centers[idx_y,:],z1[:,idx_y,idx_time_init+i])
plt.xlim(-0.3,0.45)
plt.ylim(-0.02,0.02)
plt.show()


for i in range(20):
    plt.plot(Xreal1_centers[idx_y,:],np.gradient(np.gradient(z1[:,idx_y,idx_time_init+i])))
plt.xlim(-0.3,0.45)
plt.ylim(-0.002,0.002)
plt.show()


plt.imshow(np.gradient(np.gradient(z1[:,idx_y,:],axis=0),axis=0),aspect='auto',vmin=-0.005,vmax=0.005)
plt.colorbar()
#plt.xlim(200,350)
plt.ylim(15,45)
plt.show()
# %% détection de fracture
#frame_fracture = 1150 # frame où ca commence à se casser en gros (à l'oeil)
frame_fracture = 760
i0 = input_params['i0']

time_idx_fracture = frame_fracture - i0 - t0
#idxmin_x = 20
#idxmax_x = 50
#idxmin_y = 6
#idxmax_y = 11

idxmin_x = 40
idxmax_x = 80
idxmin_y = 12
idxmax_y = 22


array_idx_profile = np.arange(idxmin_y,idxmax_y+1,1)


#%% cellule pas utilisée
"""
pour quantifier le shear :
idx_profile = 10
i = 5
gradVz1prof = np.gradient(Vz1[:,idx_profile,i],Xreal1_centers[idx_profile,:]) 

plt.plot(Xreal1_centers[idx_profile,:],gradVz1prof,'o')
"""

"""
from scipy import interpolate
def deriv1(f,x,h,smooth = 1):
    return (f(x+(smooth * h))-f(x))/ (smooth * h)

t0 = 200
dt = 50

F = np.zeros(dt*2, dtype = object)
dF = np.zeros(dt*2, dtype = object)

smooth = 3
x_new = Xreal1_centers[10,smooth+1:-smooth-1]
dx = Xreal1_centers[10,1] - Xreal1_centers[10,0]
for t in range (t0 - dt, t0 + dt) :
    Fz = interpolate.interp1d( Xreal1_centers[10,:],Vz1[:,10,t - t0 - dt])
    F[t - t0 - dt] = Fz
    dF[t - t0 - dt] = deriv1(Fz, x_new, dx, smooth)


plt.figure()
for i in range(100):
    plt.plot(Xreal1_centers[10,:],F[i](Xreal1_centers[10,:]))
plt.xlim(Xreal1_centers[10,:][idxmin_x],Xreal1_centers[10,:][idxmax_x])
#plt.plot(x_new,dF[50])
"""

# %% mesure de la hauteur max à la fracture


# premiere methode : prendre tout simplement le maximum des courbes (tres approximatif)
nb_time_indices = 10
colors = cm.winter(np.linspace(0, 0.8, nb_time_indices))
maxtemp=0
for idx_profile in array_idx_profile:
    for i in range(nb_time_indices):
        plt.plot(Xreal1_centers[idx_profile,idxmin_x:idxmax_x],z1[idxmin_x:idxmax_x,idx_profile,i+time_idx_fracture],color=colors[i])#,label=str(time_idx_fracture+i0+t0+i))
        max_i = np.max(np.abs(z1[idxmin_x:idxmax_x,idx_profile,i+time_idx_fracture]))
        plt.plot(Xreal1_centers[idx_profile,idxmin_x:idxmax_x],gaussian_filter(z1[idxmin_x:idxmax_x,idx_profile,i+time_idx_fracture],1.5),'--',color='gray')

        if max_i>maxtemp:
            maxtemp = max_i
    plt.title('idx y profile = '+str(idx_profile))
    print(maxtemp,' meters')  
    plt.legend()
    plt.show()



# deuxieme methode : d'abord

idx_profile = idx_y
i = 10

#idxmin_x_fit =  idxmin_x + 8
#idxmax_x_fit =  idxmax_x - 12

idxmin_x_fit =  idxmin_x
idxmax_x_fit =  idxmax_x

popt,pcov = curve_fit(lambda x,a,b,c:a*x**2+b*x+c,Xreal1_centers[idx_profile,idxmin_x_fit:idxmax_x_fit],z1[idxmin_x_fit:idxmax_x_fit,idx_profile,i+time_idx_fracture])

plt.plot(Xreal1_centers[idx_profile,idxmin_x:idxmax_x],z1[idxmin_x:idxmax_x,idx_profile,i+time_idx_fracture],'o')
plt.plot(Xreal1_centers[idx_profile,idxmin_x_fit:idxmax_x_fit],popt[0]*Xreal1_centers[idx_profile,idxmin_x_fit:idxmax_x_fit]**2+popt[1]*Xreal1_centers[idx_profile,idxmin_x_fit:idxmax_x_fit]+popt[2])

# plt.plot(gaussian_filter(z1[idxmin_x:idxmax_x,idx_profile,i+time_idx_fracture],1.5),'o')

#plt.ylim(-0.001,0.015)
plt.show()

# %%
dt_fracture = 50
plt.pcolormesh(np.arange(time_idx_fracture-dt_fracture,time_idx_fracture+dt_fracture,1) + t0 + i0 , Xreal1_centers[idx_profile,idxmin_x:idxmax_x],z1[idxmin_x:idxmax_x,idx_profile,time_idx_fracture-dt_fracture:time_idx_fracture+dt_fracture])
plt.xlabel('frame')
plt.ylabel('xreal (m)')
plt.colorbar()
plt.show()

plt.pcolormesh(np.arange(time_idx_fracture-dt_fracture,time_idx_fracture+dt_fracture,1) + t0 + i0 , Xreal1_centers[idx_profile,idxmin_x:idxmax_x],np.gradient(np.gradient(z1[idxmin_x:idxmax_x,idx_profile,time_idx_fracture-dt_fracture:time_idx_fracture+dt_fracture],axis=0),axis=0))
plt.xlabel('frame')
plt.ylabel('xreal (m)')

# %%
plt.pcolormesh(np.arange(time_idx_fracture-dt_fracture,time_idx_fracture+dt_fracture,1) + t0 + i0 , Xreal1_centers[idx_profile,idxmin_x:idxmax_x],np.gradient(Vz1[idxmin_x:idxmax_x,idx_profile,time_idx_fracture-dt_fracture:time_idx_fracture+dt_fracture],axis=0),vmin=-0.02,vmax=0.02)
plt.xlabel('frame')
plt.ylabel('xreal (m)')
plt.colorbar()
# %% Plot the curvature of Vz vs amplitude of Vz to see if if satisfies the linear model
freq_acq = input_params['freq_acq']
T_exc_frames = T_exc * float(freq_acq)

array_indices_maxamplitude = np.arange(1,100,16)#np.array([1,18,18+16,18+16+17,18+16+17+16,18+16+17+16+17,18+16+17+16+17+16])
first_maxamp_idx = 1
shift = 0

nb_profiles = 40
PROFILES = []
colors = cm.viridis(np.linspace(-0.5, 0.8, nb_profiles))
for i in range(first_maxamp_idx,nb_profiles,1):
    profile0 = Vz1[:,idx_y,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) - 1]
    profile1 = Vz1[:,idx_y,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int)]
    profile2 = Vz1[:,idx_y,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) + 1]
    
    profiles = [profile0,profile1,profile2]
    array_sums = np.array([np.sum(profile0[idxmin_x:idxmax_x]),np.sum(profile1[idxmin_x:idxmax_x]),np.sum(profile2[idxmin_x:idxmax_x])])    
    index_choice = np.argmax(array_sums)
    plt.plot(Xreal1_centers[idx_y,:],profiles[index_choice],color=colors[i])
    shift += index_choice - 1
    PROFILES.append(profiles[index_choice])
plt.xlim(-0.35,0.35)
plt.ylim(-0.01,0.046)
plt.show()

PROFILES = np.array(PROFILES)
amplitudes = np.max(PROFILES[:,idxmin_x:idxmax_x],axis=1)
indices_maxes = np.argmax(PROFILES[:,idxmin_x:idxmax_x],axis=1) + idxmin_x
semilengthfit = 15
semilengthfit_broken = 5
POPT = []
PCOV = []
second_derivatives = []

for i in range(PROFILES.shape[0]):
    profile = PROFILES[i,:]

    second_derivatives.append(np.gradient(np.gradient(profile)))

    print(indices_maxes[i]-semilengthfit)
    print(indices_maxes[i]+semilengthfit)
    profile2fit = profile[indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    x2fit = Xreal1_centers[idx_y,indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    
    popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    if np.sqrt(np.diag(pcov)[0]) > 0.03:
            profile2fit = profile[indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            x2fit = Xreal1_centers[idx_y,indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    POPT.append(popt)
    PCOV.append(pcov)
    plt.plot(Xreal1_centers[idx_y,idxmin_x:idxmax_x],PROFILES[i,idxmin_x:idxmax_x],'.',alpha=0.1,color=colors[i])
    plt.plot(x2fit,(popt[0]*0.5)*x2fit**2+popt[1]*x2fit+popt[2],color=colors[i])
    
plt.show()

POPT = np.array(POPT)
PCOV = np.array(PCOV)
second_derivatives = np.array(second_derivatives)
curvatures = POPT[:,0].flatten()

plt.title('curvature of Vz   vs amplitude')
plt.plot(amplitudes,np.abs(curvatures),'o',label='data')
wavelength = average_wavelength
wavenumber = 2*np.pi/wavelength
ampvals = np.linspace(np.min(amplitudes),np.max(amplitudes))
plt.plot(ampvals , wavenumber**2 * ampvals,label='$k_{avg}^{2}\cdot$amplitude')
plt.xlim(0,0.05)
#plt.ylim(0,1.5)
plt.xlabel('amplitude (m/s)')
plt.ylabel('Vz curvature (m^-1 . s^-1)')
plt.legend()
plt.show()


# %% measure the amplitude (and the curvature) just before fracture

nb_time_indices = 14
colors = cm.winter(np.linspace(0, 0.8, nb_time_indices))
maxtemp=0
idx_profile = 16
extremum = 'max' # preciser si la derniere amplitude maximale pour la vitesse correspond à un min ou un max
for i in range(nb_time_indices):
    plt.plot(Xreal1_centers[idx_profile,x_indices2plot],Vz1[x_indices2plot,idx_profile,i+time_idx_fracture],color=colors[i])#,label=str(time_idx_fracture+i0+t0+i))
    if extremum=='max':
        max_i = np.max(Vz1[x_indices2plot,idx_profile,i+time_idx_fracture][10:-20])
    elif extremum=='min':
        max_i = np.min(Vz1[x_indices2plot,idx_profile,i+time_idx_fracture][10:-20])
#    plt.plot(Xreal1_centers[idx_profile,x_indices2plot],gaussian_filter(Vz1[x_indices2plot,idx_profile,i+time_idx_fracture],1.5),'--',color='gray')

    if max_i>maxtemp:
        maxtemp = max_i
plt.title('idx y profile = '+str(idx_profile))
print(maxtemp,' m/s')  
plt.legend()
plt.show()

# dans l'approximation linéaire des ondes, la mesure de l'amplitude suffit pour déduire la courbure

print('correspond for the linear model to : amplitude = ',maxtemp / (2*np.pi*f_exc),' meters')
# %%
# focale mesurée par calibration : 1511 => focale en mm = 5.5e-6 * 1511 * 1e3