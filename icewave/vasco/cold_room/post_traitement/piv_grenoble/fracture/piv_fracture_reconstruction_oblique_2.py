#%% -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import glob
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import time 

import h5py
import os

#kerdir = os.getcwd()
#os.chdir('/home/vzanchi/Bureau/Turbotice_git/')

#sys.path.append('C:/Users/Vasco/OneDrive - Université de Paris/Documents/git/icewave')
sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')

import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp
#os.chdir(kerdir)

ordi = 'dell_vasco'

#kerdir = os.getcwd()
if ordi=='Leyre':
    os.chdir('/run/user/1003/gvfs/smb-share:server=thiou.local,share=homes/vasco/Grenoble_Nov2024/post_traitement/piv_grenoble/fracture/python_functions/')
elif ordi=='dell_vasco':
    sys.path.append('Z:/vasco/Grenoble_Nov2024/post_traitement/piv_grenoble/fracture/python_functions/')
print(os.listdir())
from organize import *
#os.chdir(kerdir)


#%% fonction

# noms des dossiers contenant les données

def reconstruct_Vz_field_from_pivdata(input_params = dict):

    date = input_params['date'] 
    acq_num = input_params['acq_num']
    f_exc = input_params['f_exc']
    freq_acq = input_params['freq_acq']
    cam_SN = input_params['cam_SN']
    delta_x = input_params['delta_x']
    delta_y = input_params['delta_y']
    imwidth = input_params['imwidth']
    imheight = input_params['imheight']
    x0 = input_params['x0']
    y0 = input_params['y0']
    pixel_size = input_params['pixel_size']
    focale_mm = input_params['focale_mm']
    focale = input_params['focale']
    h_camera = input_params['h_camera']
    alpha_0_deg = input_params['alpha_0_deg']
    alpha_0 = input_params['alpha_0']
    unit_input_velocity = input_params['unit_input_velocity']
    box_position_information = input_params['box_position_information']
    a = input_params['a']
    W = input_params['W']
    Dt = input_params['Dt']
    i0 = input_params['i0']
    N = input_params['N']
    if ordi=='Leyre':
        base = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    elif ordi=='dell_vasco':
        base =f'L:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    path2data = f'{base}matData/'
    matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_Dt{Dt}_b1_W{W}_full_total_processed.mat'
    
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

# pour la manip du 29/11 acq7 :
#csv_file_path = "/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241129/manip_fracture/Acquisition_7/input_params.csv"
# pour la manip du 28/11 acq3 :
if ordi=='Leyre':
    csv_file_path = "/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/20241128/manip_fracture/Acquisition_3/input_params.csv"
elif ordi=='dell_vasco':
    csv_file_path = "L:/Gre24/Data/20241128/manip_fracture/Acquisition_3/input_params.csv"


# load input params
input_params = csv2dict(csv_file_path)
# on rajoute quelquues params qui dependent des parameters d'entrée :
input_params['x0'] = input_params['sensor_width']/2 - input_params['delta_x']
input_params['y0'] = input_params['sensor_height']/2 - input_params['delta_y']
input_params['focale'] = input_params['focale_mm']*1e-3/input_params['pixel_size'] # distance focale convertie en nombre de pixels du capteur
input_params['alpha_0'] = input_params['alpha_0_deg']*np.pi/180
###############################

mat_dict,t,Xreal1,Yreal1,Vz1,pc1,Xreal1_centers,Yreal1_centers = reconstruct_Vz_field_from_pivdata(input_params=input_params)

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




#%% input values for postprocessing


input_values_postprocessing = {}


# pour le 29/11 acq7:
input_values_postprocessing['wavelength_measurement_statio'] = {}
input_values_postprocessing['wavelength_measurement_statio']['framestart'] = 0
input_values_postprocessing['wavelength_measurement_statio']['framestop'] = 300
input_values_postprocessing['wavelength_measurement_statio']['idx_profile'] = 16
input_values_postprocessing['wavelength_measurement_statio']['xinf'] = -0.50
input_values_postprocessing['wavelength_measurement_statio']['xsup'] = 0.50

input_values_postprocessing['Vz_measurement_pre_fracture'] = {}
input_values_postprocessing['Vz_measurement_pre_fracture']['idx_profile'] = 16
input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture'] = 760
input_values_postprocessing['Vz_measurement_pre_fracture']['xinf'] = -0.30
input_values_postprocessing['Vz_measurement_pre_fracture']['xsup'] = 0.30

"""
# pour le 28/11 acq3:
input_values_postprocessing['wavelength_measurement_statio'] = {}
input_values_postprocessing['wavelength_measurement_statio']['framestart'] = 2000
input_values_postprocessing['wavelength_measurement_statio']['framestop'] = 2500
input_values_postprocessing['wavelength_measurement_statio']['idx_profile'] = 16
input_values_postprocessing['wavelength_measurement_statio']['xinf'] = -0.50
input_values_postprocessing['wavelength_measurement_statio']['xsup'] = 0.50

input_values_postprocessing['Vz_measurement_pre_fracture'] = {}
input_values_postprocessing['Vz_measurement_pre_fracture']['idx_profile'] = 16
input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture'] = 2000
input_values_postprocessing['Vz_measurement_pre_fracture']['xinf'] = -0.30
input_values_postprocessing['Vz_measurement_pre_fracture']['xsup'] = 0.30
"""
# %% mesure longueur d'onde pour differents y (devraient être égales)

# autre methode : mesurer pour chaque profile puis moyenner


framestart = input_values_postprocessing['wavelength_measurement_statio']['framestart']
framestop = input_values_postprocessing['wavelength_measurement_statio']['framestop']
idx_profile = input_values_postprocessing['wavelength_measurement_statio']['idx_profile']
xinf = input_values_postprocessing['wavelength_measurement_statio']['xinf']
xsup = input_values_postprocessing['wavelength_measurement_statio']['xsup']

x_indices2plot = np.where((Xreal1_centers[idx_profile,:]>xinf)&(Xreal1_centers[idx_profile,:]<xsup))[0]

wavelengths_measured = []
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
            dX = Xreal1_centers[idx_profile,minima_indices[1]+x_indices2plot[1]] - Xreal1_centers[idx_profile,minima_indices[0]+x_indices2plot[0]]
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
# faire subpix pour améliorer
# %% Plot the curvature of Vz vs amplitude of Vz to see if if satisfies the linear model

idx_profile = input_values_postprocessing['Vz_measurement_pre_fracture']['idx_profile']
frame_fracture = input_values_postprocessing['Vz_measurement_pre_fracture']['frame_fracture']
xinf = input_values_postprocessing['Vz_measurement_pre_fracture']['xinf']
xsup = input_values_postprocessing['Vz_measurement_pre_fracture']['xsup']
x_indices2plot = np.where((Xreal1_centers[idx_profile,:]>xinf)&(Xreal1_centers[idx_profile,:]<xsup))[0]


freq_acq = input_params['freq_acq']
f_exc = input_params['f_exc']
T_exc = 1/f_exc
T_exc_frames = T_exc * float(freq_acq)

#array_indices_maxamplitude = np.arange(1,100,16)#np.array([1,18,18+16,18+16+17,18+16+17+16,18+16+17+16+17,18+16+17+16+17+16])
first_maxamp_idx = 1
shift = 0

nb_profiles = 26
PROFILES = []
colors = cm.viridis(np.linspace(-0.5, 0.8, nb_profiles))
for i in range(first_maxamp_idx,nb_profiles,1):
    profile0 = Vz1[:,idx_profile,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) - 1]
    profile1 = Vz1[:,idx_profile,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int)]
    profile2 = Vz1[:,idx_profile,np.round(first_maxamp_idx + shift + i * T_exc_frames).astype(int) + 1]
    
    profiles = [profile0,profile1,profile2]
    array_sums = np.array([np.sum(profile0[x_indices2plot]),np.sum(profile1[x_indices2plot]),np.sum(profile2[x_indices2plot])])    
    index_choice = np.argmax(array_sums)
    plt.plot(Xreal1_centers[idx_profile,x_indices2plot],profiles[index_choice][x_indices2plot],color=colors[i])
    shift += index_choice - 1
    PROFILES.append(profiles[index_choice])
#plt.xlim(-0.35,0.35)
#plt.ylim(-0.01,0.046)
plt.show()

PROFILES = np.array(PROFILES)
amplitudes = np.max(PROFILES[:,x_indices2plot],axis=1)
indices_maxes = np.argmax(gaussian_filter(PROFILES[:,x_indices2plot],1.5,axes=(1)),axis=1) + np.min(x_indices2plot)
semilengthfit = 15
semilengthfit_broken = 5
POPT = []
PCOV = []
second_derivatives = []

for i in range(PROFILES.shape[0]):
    profile = PROFILES[i,:]

    second_derivatives.append(np.gradient(np.gradient(profile)))

   # print(indices_maxes[i]-semilengthfit)
   # print(indices_maxes[i]+semilengthfit)
    profile2fit = profile[indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    x2fit = Xreal1_centers[idx_profile,indices_maxes[i] - semilengthfit : indices_maxes[i] + semilengthfit]
    
    popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    if np.sqrt(np.diag(pcov)[0]) > 0.05:
            profile2fit = profile[indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            x2fit = Xreal1_centers[idx_profile,indices_maxes[i] - semilengthfit_broken : indices_maxes[i] + semilengthfit_broken]
            popt,pcov = curve_fit(lambda x,a,b,c : (a*0.5)*x**2+b*x+c , x2fit , profile2fit)
    POPT.append(popt)
    PCOV.append(pcov)
    plt.plot(Xreal1_centers[idx_profile,x_indices2plot],PROFILES[i,x_indices2plot],'.',alpha=0.3,color=colors[i])
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
#plt.legend()
plt.show()

# dans l'approximation linéaire des ondes, la mesure de l'amplitude suffit pour déduire la courbure

print('correspond for the linear model to : amplitude = ',maxtemp / (2*np.pi*f_exc),' meters')



# %% choix de la frame de reference (à laquelle le profil est nul)




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

plt.ylim(-0.01,0.01)
#plt.xlim(0,15)
plt.grid()
plt.legend()
##########
plt.show()

z1 = cumtrapVz1





# %% mesure de la hauteur max à la fracture


# premiere methode : prendre tout simplement le maximum des courbes (tres approximatif)
nb_time_indices = 10
colors = cm.winter(np.linspace(0, 0.8, nb_time_indices))
maxtemp=0
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


plt.pcolormesh(np.arange(time_idx_fracture-dt_fracture,time_idx_fracture+dt_fracture,1) + t0 + i0 , Xreal1_centers[idx_profile,idxmin_x:idxmax_x],np.gradient(Vz1[idxmin_x:idxmax_x,idx_profile,time_idx_fracture-dt_fracture:time_idx_fracture+dt_fracture],axis=0),vmin=-0.02,vmax=0.02)
plt.xlabel('frame')
plt.ylabel('xreal (m)')
plt.colorbar()
