#%% -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:29:21 2024
@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt 
import glob
from scipy.interpolate import RegularGridInterpolator
import time 

import h5py

import icewave.tools.datafolders as df
import icewave.drone.drone_projection as dp

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



    
#%% noms des dossiers contenant les données
date = '0211'
year = '2024'
drone_ID = 'bernache'
flight_ID = 'structured_data'

base = df.find_path('Hublot24')
base = 'H:/Rimouski_2024/Data/'
path2data = f'{base}{year}/{date}/Drones/{drone_ID}/matData/{flight_ID}/'
path2data = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Share_hublot/Data/0211/Drones/bernache/matData/stereo_001/structured_data/'
filelist_mat = glob.glob(f'{path2data}*.mat')

structured_mat = glob.glob(f'{path2data}*scaled.mat')[0]


#%%
# Load the .mat file and convert to a dictionary

with h5py.File(structured_mat, 'r') as fmat:
    mat_dict = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    mat_dict = mat_to_dict(fmat['m'],fmat['m'])
    mat_dict['PIV_param']['p_param'] = matcell2dict_PIV(mat_dict['PIV_param']['p_param'])
    mat_dict['PIV_param']['s_param'] = matcell2dict_PIV(mat_dict['PIV_param']['s_param'])

#%% Convert Vy displacement field to vertical field Vz 
key_pix = 'PIXEL'

Vy = mat_dict['Vy'] # time, y, x
Vy_transposed = np.transpose(Vy,(2,1,0)) 

fps = 1/mat_dict['SCALE']['ft']
Dt = mat_dict['PIV_param']['Dt'] # time step between two image that are compared using PIV
t = mat_dict['t']*fps
x = mat_dict[key_pix]['x_pix'] # pixel position of each PIV box
y = mat_dict[key_pix]['y_pix'] # pixel position of each PIV box 
# compute interpolation of pixel displacement vertical velocity field 
Fy = RegularGridInterpolator((x,y,t), Vy_transposed)


#%%


Vz = dp.vertical_velocity_from_pixel_displacement(Vy_transposed,x,y,t,mat_dict[key_pix]['y0'],mat_dict['DRONE']['h_drone'],
                                               mat_dict['DRONE']['alpha_0'],mat_dict['DRONE']['focale'],fps,Dt)
    



#%% Plot a plt.colormesh 
W = mat_dict['PIV_param']['w']
x_edges = np.resize(x,x.size + 1)
y_edges = np.resize(y,y.size + 1)
x_edges[-1] = x_edges[-2] + W
y_edges[-1] = y_edges[-2] + W
x_edges = x_edges - W/2
y_edges = y_edges - W/2


Xedges,Yedges = np.meshgrid(x_edges,y_edges, indexing = 'ij')

Xreal,Yreal = dp.projection_real_space(Xedges, Yedges, mat_dict[key_pix]['x0'], mat_dict[key_pix]['y0'], mat_dict['DRONE']['h_drone']
                                       ,mat_dict['DRONE']['alpha_0'],mat_dict['DRONE']['focale'])

fig, ax = plt.subplots()
c = ax.pcolormesh(Xreal,Yreal,Vz[:,:,1500],shading = 'auto')
fig.colorbar(c,ax = ax)






#%%
##########################################
################## fin code Sebastien ####
##########################################


##########################################
################## debut code vasco ######
##########################################

from scipy.interpolate import griddata
from matplotlib.animation import FuncAnimation
import scipy
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

#from mat73 import loadmat
# pour obtenir quelque chose en unité de longueur on fait : delta_z = Vz*Dt/fps


#%% tests manip frigo vasco

# noms des dossiers contenant les données

def fonction_test(cam_SN = '40307970',box_position_information=False):

    f_exc = 10
    freq_acq = 57.7
    #cam_SN = '40300722'

    imwidth = 1920
    imheight = 1200
    focale = 4706 # distance focale convertie en nombre de pixels du capteur (ici distance focal = 16 mm)

    h_camera = 0.13 # hauteur de la camera en metres par rapport à la surface de l'eau
    alpha_0_deg = 12.5
    alpha_0 = alpha_0_deg*np.pi/180

    base = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=data/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/20241105/manip_f_k_membrane/Acquisition_1/camera_{cam_SN}/{f_exc}Hz_{freq_acq}Hz/'
    path2data = f'{base}matData/'

    W = 32
    Dt = 20

    matfile = f'{path2data}PIV_processed_i00_N0_Dt{Dt}_b1_W{W}_full_total_processed.mat'

    with h5py.File(matfile, 'r') as fmat:
        mat_dict = {}
        
        print('Top-level keys : ', list(fmat.keys()))

        mat_dict = mat_to_dict(fmat['m'],fmat['m'])

    print("KEYS : ",mat_dict.keys())
    # Convert Vy displacement field to vertical field Vz 

    Vy = mat_dict['Vy'] # time, y, x
    Vy_transposed = np.transpose(Vy,(2,1,0)) 

    fps = freq_acq
    t = np.arange(len(mat_dict['t']))*fps # unit ?
    if box_position_information==False:
        x = np.arange(len(mat_dict['x']))*(W/2) + W + 0.5 # pixel position of each PIV box
        y = np.arange(len(mat_dict['y']))*(W/2) + W + 0.5 # pixel position of each PIV box 
    else:
        x = mat_dict['xpix']
        y = mat_dict['ypix']
    if  imwidth%2==0:
        x0 = imwidth/2 + 0.5
    else:
        x0 = (imwidth+1)/2 + 0.5

    if  imheight%2==0:
        y0 = imheight/2 + 0.5
    else:
        y0 = (imheight+1)/2 + 0.5



    # compute interpolation of pixel displacement vertical velocity field 
    Fy = RegularGridInterpolator((x,y,t), Vy_transposed)

    # compute vertical velocity Vz

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
    pc = ax.pcolormesh(Xreal,Yreal,Vz[:,:,20],shading = 'auto',vmin=-0.003,vmax=0.003) # pcolormesh object !
    fig.colorbar(c,ax = ax)
    #ax.set_ylim(-0.15,0.1)
    #ax.set_xlim(-0.17,0.17)

    X,Y = np.meshgrid(x,y)

    Xreal_centers,Yreal_centers = dp.projection_real_space(X, Y, x0, y0, h_camera, alpha_0, focale)

    return Xreal,Yreal,Vz,pc,Xreal_centers,Yreal_centers

Xreal1,Yreal1,Vz1,pc1,Xreal1_centers,Yreal1_centers = fonction_test(cam_SN = '40307970')
Xreal2,Yreal2,Vz2,pc2,Xreal2_centers,Yreal2_centers = fonction_test(cam_SN = '40300722')

delta_x_cam = 0.184 # distance en metres entre les deux cameras
Xreal2_shift = Xreal2 + delta_x_cam

plt.figure()
plt.pcolormesh(Xreal1,Yreal1,Vz1[:,:,20],shading = 'auto',vmin=-0.003,vmax=0.003)
plt.pcolormesh(Xreal2_shift,Yreal2,Vz2[:,:,20],shading = 'auto',vmin=-0.003,vmax=0.003)
plt.colorbar()
plt.ylim(-0.15,0.1)
plt.xlim(-0.17,0.32)
plt.show()

#%% creer une video


t = np.arange(Vz1.shape[2])

def update(frame):
    plt.cla()  # Efface le contenu de l'axe à chaque trame
#    plt.imshow(matrix[frame,:,:])
    plt.pcolormesh(Xreal1,Yreal1,Vz1[:,:,frame],shading = 'auto',vmin=-0.003,vmax=0.003)
    plt.pcolormesh(Xreal2_shift,Yreal2,Vz2[:,:,frame],shading = 'auto',vmin=-0.003,vmax=0.003)
    plt.xlabel('x (m)',fontsize=15)
    plt.ylabel('y (m)',fontsize=15)
    plt.ylim(-0.15,0.1)
    plt.xlim(-0.17,0.32)
# Crée la figure et l'axe pour l'animation
fig, ax = plt.subplots() 
frames_to_include = range(0, len(t), 1)
ani = FuncAnimation(fig, update, frames=frames_to_include, repeat=False,blit=False)

output_file = '/run/user/1003/gvfs/smb-share:server=adour.local,share=data/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/20241105/manip_f_k_membrane/Acquisition_1/video_PIV2.gif'
# Save the animation as a video (specify the writer and additional options)
ani.save(output_file, writer='ffmpeg', fps=10,dpi=300)  # Adjust the FPS as needed
plt.show()

#%% mettre dans un meme tableau les donnees des deux pcolormesh (en interpolant)
#data_all_pc = np.hstack((pc1.get_array().data.flatten(),pc2.get_array().data.flatten()))
def interpolate_twoPIVfields_oneframe(Vz1=Vz1,Vz2=Vz2,Xreal1_centers=Xreal1_centers,Xreal2_centers=Xreal2_centers,Yreal1_centers=Yreal1_centers,Yreal2_centers=Yreal2_centers,delta_x_cam=delta_x_cam,frame=1,nb_pts_interp=1000,plot=False):

    data_all_pc = np.hstack((Vz1[:,:,frame].flatten(),Vz2[:,:,frame].flatten()))
    Xreal2_centers_shift = Xreal2_centers + delta_x_cam

    Xreal_flat_all_pc = np.hstack((np.transpose(Xreal1_centers).flatten(),np.transpose(Xreal2_centers_shift).flatten()))
    Yreal_flat_all_pc = np.hstack((np.transpose(Yreal1_centers).flatten(),np.transpose(Yreal2_centers).flatten()))

    
    # Define a regular grid to interpolate the data onto
    xi = np.linspace(np.min(Xreal_flat_all_pc), np.max(Xreal_flat_all_pc), nb_pts_interp) # comment choisir le nb de points?
    yi = np.linspace(np.min(Yreal_flat_all_pc), np.max(Yreal_flat_all_pc), nb_pts_interp)
    xgrid, ygrid = np.meshgrid(xi, yi)

    # Interpolate
    zi = griddata((Xreal_flat_all_pc, Yreal_flat_all_pc), data_all_pc, (xgrid, ygrid), method='linear')  # or 'cubic', 'nearest'

    if plot==True:
        # Plot the result
        plt.figure(figsize=(20,10))
        plt.imshow(zi, extent=(np.min(Xreal_flat_all_pc), np.max(Xreal_flat_all_pc), np.min(Yreal_flat_all_pc), np.max(Yreal_flat_all_pc)), origin='lower', cmap='viridis',vmin=-0.003,vmax=0.003)
        #plt.scatter(x, y, c=data_all_pc, edgecolors='k')  # Original points
        plt.colorbar()
        plt.show()    
    return zi,xi,yi,xgrid,ygrid

#
def interpolate_twoPIVfields_allframes(Vz1=Vz1,Vz2=Vz2,Xreal1_centers=Xreal1_centers,Xreal2_centers=Xreal2_centers,Yreal1_centers=Yreal1_centers,Yreal2_centers=Yreal2_centers,delta_x_cam=delta_x_cam,nb_pts_interp=1000,plot=False):
    nb_frames = Vz1.shape[2]
    VZi = np.zeros((nb_pts_interp,nb_pts_interp,nb_frames))
    for i in range(nb_frames):
        zi,xi,yi,xgrid,ygrid = interpolate_twoPIVfields_oneframe(Vz1=Vz1,Vz2=Vz2,Xreal1_centers=Xreal1_centers,Xreal2_centers=Xreal2_centers,Yreal1_centers=Yreal1_centers,Yreal2_centers=Yreal2_centers,delta_x_cam=delta_x_cam,frame=i,nb_pts_interp=nb_pts_interp,plot=plot)
        VZi[:,:,i] = zi
    return VZi,xi,yi,xgrid,ygrid

# test :
interpolate_PIVfields_oneframe(plot=True)

#%%
VZi,xi,yi,xgrid,ygrid = interpolate_twoPIVfields_allframes()

#%% demoduler le champ VZi

def demodulate_VZi(VZi=VZi):
    #H = np.transpose(VZi,(2,0,1)) # t,x,y (alors que VZi est sous la forme (x,y,t))
    H = VZi
    Y = np.zeros((H.shape[0],H.shape[1]))

    Nframes = H.shape[2]
    t = (1/freq_acq) * np.arange(Nframes)
    phase = np.exp(-2*np.pi*1j*f_exc*t) # dans le code d'antonin il n'y a pas le "-"

    # function np.tile(arr,(l,m,n)) replicates the array arr l times along axis 0, m times along axis 1, n times along axis 2
    demod = np.tile(phase,(H.shape[0],H.shape[1],1)) # on met 1 dans la troisieme direction (axis 2) car phase a déjà la dimension H.shape[2]

    output = np.sum(H*demod, axis=2) 

    return output

complex_field = demodulate_VZi()

#%% maintenant utiliser le code fait pendant le stage pour calculer la vitesse de phase

general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=data/thiou/storageshared/Banquise/Vasco/Frigo_pmmh/20241105/manip_f_k_membrane/Acquisition_1/camera_{cam_SN}/'

ll = os.listdir(general_folder)
ld = []
list_f_exc = []
list_freq_acq = []
for l in ll:
    if ('Hz_' in l)&(l[-2:]=='Hz'):
        ld.append(l)
        print(l)
        list_f_exc.append(int(l[:l.find('Hz')]))
        if '.' in l[l.find('_')+1:-2]:
            list_freq_acq.append(float(l[l.find('_')+1:-2]))
        elif ('.' in l[l.find('_')+1:-2])==False:
            list_freq_acq.append(int(l[l.find('_')+1:-2]))



tab_f_exc = np.array(list_f_exc)
tab_freq_acq = np.array(list_freq_acq,dtype='object')
indsort = np.argsort(tab_f_exc)
tab_f_exc = tab_f_exc[indsort]
tab_freq_acq = tab_freq_acq[indsort]

print('tab_f_exc',tab_f_exc)
print('tab_freq_acq',tab_freq_acq)

#%% correlation pour spatial shift
def detect_shift_subpix(signal1=np.sin(np.arange(10)),signal2=np.sin(np.arange(10)+1),dx=1):
    cross_corr = scipy.signal.correlate(signal2,signal1)
    lags = scipy.signal.correlation_lags(len(signal2),len(signal1))
    plt.figure()
    plt.plot(signal1,label='signal1')
    plt.plot(signal2,label='signal2')
    plt.plot(lags,cross_corr,label='cross correlation')
    plt.legend()
    plt.show()

    idx_max = np.argmax(cross_corr)

    # detection subpixellaire du max :
    x1 = lags[idx_max-1]
    y1 = cross_corr[idx_max-1]
    x2 = lags[idx_max]
    y2 = cross_corr[idx_max]
    x3 = lags[idx_max+1]
    y3 = cross_corr[idx_max+1]

    A = np.array([[x1**2,x1,1],[x2**2,x2,1],[x3**2,x3,1]])
    B = np.array([y1,y2,y3])
    Coefficients = np.linalg.solve(A,B)

    x_max_subpix = -Coefficients[1]/(2*Coefficients[0])
#    dx = xdata[1]-xdata[0] # attention ne marche que si abscisses espacés régulièrement
    return x_max_subpix , x_max_subpix * dx

#test:

nb_shifts = 10
tab_shifts = np.zeros(nb_shifts)

for i in range(len(tab_shifts)):
    signal1 = np.real(complex_field[200,300:600]*np.exp(2*np.pi*1j*(1.54)))
    signal2 = np.real(complex_field[200,300:600]*np.exp(2*np.pi*1j*(1.54 + ((i+1)/nb_shifts)*(1/10))))
    shift = detect_shift_subpix(signal1=signal1,signal2=signal2)[0]
    tab_shifts[i] = shift

plt.plot(tab_shifts,'o')

#%% bof...
'''for i in range(len(tab_f_exc)):

    f_exc = tab_f_exc[i]
    freq_acq = tab_freq_acq[i]

    #W = 64
    #Dt = 4
    #Dt = tab_Dt[i]
  #  if not('v_phase' in dico_dataset):
   #     dico_dataset['v_phase'] = []
    #    dico_dataset['v_phase_err'] = []
    #    dico_dataset['xlim_fit_px'] = []
    #    dico_dataset['tab_f_exc'] = []
    #    dico_dataset['tab_freq_acq'] = []

    #if dataset==0:
    #    folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
    #elif dataset!=0:
    #    folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData'+str(dataset)+'/video_demod_W'+ str(W) +'_Dt' +str(Dt)

    #matfile_path = folder_video_demod + '/figdata_complex.mat'
    #print('data loading... '+str(f_exc))
    #Data_demod = load_complex_field(matfile_path)
    #print('data loaded!')
    #complex_field = np.copy(Data_demod.data)
    #complex_field = DCF[str(f_exc)] # affect the corresponding complex field to the variable complex_field
    print('imshow of the real part of complex field for f_exc = '+str(f_exc)+' Hz')
    plt.figure()
    #notnan_mask = np.where(np.isnan(complex_field)==False)
    #plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field[notnan_mask])))

    plt.show()

    nb_periods = 4
    nb_pts_per_period = 256
    index_profile_line = 200
    m = nb_periods*nb_pts_per_period
    matrix = np.zeros((m,complex_field.shape[1]))
 #   if remove_left_part==True:
 #       matrix = np.zeros((m,complex_field.shape[1]-5))
    for i in range(m):
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))#/np.max(np.abs(complex_field))
        #img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
 #       if remove_left_part==True:
 #           img = img[:,5:]
 #       if do_savgol==True:
 #           img = savgol_filter(img,21,3)
        p = img[index_profile_line,:]
        #plt.imshow(img)
        #plt.pause(0.01)
        matrix[i,:] = p

    print('profile vs time at line '+str(index_profile_line)+' computed.')


    X_MAXS = np.arange(matrix.shape[1])
    Y_MAXS = []
    for i in range(matrix.shape[1]):
        indices = find_peaks(matrix[:,i])[0]
        if len(indices)==3:
            indices=np.hstack((indices,np.array([np.inf])))
        i0 = indices[0]
        i1 = indices[1]
        i2 = indices[2]
        #i3 = indices[3]        
#        tab_max_along_column = np.array([i0,i1,i2])

        if len(Y_MAXS)==0:
            Y_MAXS.append(i1)
        else:
            
            index_closest = np.argmin(abs(Y_MAXS[-1] - indices))
            Y_MAXS.append(indices[index_closest])
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=nb_pts_per_period/3)):
    #    Y_MAXS.append(i0)
    #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>nb_pts_per_period/3)):
    #    Y_MAXS.append(i1)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<abs(Y_MAXS[-1]-i1))):
        #    Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>abs(Y_MAXS[-1]-i1))):
        #    Y_MAXS.append(i1)
        #else:
        #    Y_MAXS.append(np.nan)
    print('indices maxima (using find_peaks) found.')
    print('plot peaks on top of imshow of profile vs time : ')

    plt.figure()
    plt.imshow(matrix,aspect='auto')
    plt.plot(Y_MAXS,'r')
    plt.colorbar()
    #plt.ylim(np.min(Y_MAXS),np.max(Y_MAXS))
    plt.show()
    



    #  fitter la partie d'intérêt
    
    
    
    xlim_fit = (10,33) # à modifier en fonciton de ce qu'on a!
    #xlim_fit = (1,20)
#    W = 64
#    dcm = 11
#    dpx =  978
#    dcm = 10
#    dpx =  501

    x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
#    x_meters = x_px * (W/2) * (dcm*1e-2)/dpx
    t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
    print(t_px)
    t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


    def linear(x,a,b):
        return a*x+b
#    if uncertainty_pixels==True:
#        opt,pcov = curve_fit(linear,t_sec,x_meters,sigma=np.ones(len(x_meters))*(48/2)* W/2 * (dcm*1e-2)/dpx)
#    else:
    opt,pcov = curve_fit(linear,t_sec,x_meters)
    plt.plot(t_sec,x_meters,'o')
    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Position of the wave front (meters)',fontsize=15)
    plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(pcov[0][0],2))+' m/s',fontsize=15)

    #plt.errorbar(t_sec,x_meters,np.ones(len(x_meters))*(48/2)* W/2 * (dcm*1e-2)/dpx)
    t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
    plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
    plt.show()
    print('phase velocity = '+str(opt[0])+' +- '+str(pcov[0][0])+' m/s')
    

'''