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
Dt = 1
i0 = 0
N = 0

date = '0506'
acq_num = 4
camera_SN = '22458101'

f_exc = 0.94
freq_acq = 20

#computer = 'adour'

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/{date}/cameras/frac/image_sequence/'


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

y_indices = np.where((ypix>ymin)&(ypix<ymax))[0]

Vx_converted_px = Vx/freq_acq
Vy_converted_px = Vy/freq_acq

plt.figure()
plt.imshow(Vy_converted_px[t_plot,y_indices,:],extent=[np.min(xpix),np.max(xpix),np.max(ypix[y_indices]),np.min(ypix[y_indices])],vmin=-np.max(Vy_converted_px)/5,vmax=np.max(Vy_converted_px)/5)
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
Vz = -Vy * 1/np.cos(alpha) # signe - car l'axe y est inversé

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
    vnorm = np.sqrt(vx**2+vy**2)
    plt.plot(vnorm*np.sign(vy))
plt.ylim(-70,70)
plt.show()




# %% obtenir une meilleur estimation de la longueur d'onde
tinit2plot = 0
# indices x and y to plot
if W==32:
    xind = 50
    yind = 10
elif W==64:
    xind = 25
    yind = 5


if len(ypix)!=len(Vy[0,:,0]):
    xpix = xpix[1:-1]
    ypix = ypix[1:-1]

array_aspect_ratios = np.zeros_like(xpix)
for i in range(len(array_aspect_ratios)):
    array_aspect_ratios[i] = compute_aspect_ratio(xpix[i],ypix[yind])


for i in range(tinit2plot,tinit2plot+30,1):
    plt.plot(np.arange(len(Vz[i,yind,:]))*(W/2) * np.abs(array_aspect_ratios),Vz[i,yind,:] * np.abs(array_aspect_ratios))
#plt.xlim(42,120)
plt.xlabel('x (cm)')
plt.ylabel('Vz (cm/sec)')
plt.show()


for i in range(tinit2plot,tinit2plot+30,1):

    vx = Vx[i,yind,:]
    vy = Vy[i,yind,:]
    vnorm = np.sqrt(vx**2+vy**2)
    plt.plot(np.arange(len(vnorm))*(W/2) * compute_aspect_ratio(xpix[xind],ypix[yind]),v*np.sign(vy))
#plt.xlim(42,120)
plt.xlabel('x (cm)')
plt.ylabel('V * sign(Vy) (px/sec)')
plt.title('test en affichant la norme * signe(Vy)')
plt.show()

#%%
Vx






# %% velocity field with arrows
def plot_velocity_simple(Vx, Vy, time_index=0):
    """Version simple pour un affichage rapide"""
    vx = Vx[time_index, :, :]
    vy = Vy[time_index, :, :]
    
    ny, nx = vx.shape
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    
    plt.figure(figsize=(10, 6))
    plt.quiver(X, Y, vx, vy)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Champ de vitesses - t={time_index}')
    plt.show()

plot_velocity_simple(Vx,Vy,time_index=100)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

class VelocityFieldAnimation:
    def __init__(self, Vx, Vy, dt=1.0, stride=2, scale=1.0, figsize=(12, 8)):
        """
        Créer une animation du champ de vitesses
        
        Parameters:
        -----------
        Vx, Vy : arrays 3D [temps, y, x]
            Composantes de vitesse
        dt : float
            Pas de temps entre les frames (pour l'affichage)
        stride : int
            Sous-échantillonnage spatial (1 = toutes les flèches)
        scale : float
            Facteur d'échelle pour les flèches
        figsize : tuple
            Taille de la figure
        """
        self.Vx = Vx
        self.Vy = Vy
        self.dt = dt
        self.stride = stride
        self.scale = scale
        
        # Dimensions
        self.nt, self.ny, self.nx = Vx.shape
        
        # Grilles de coordonnées (sous-échantillonnées)
        x = np.arange(0, self.nx, stride)
        y = np.arange(0, self.ny, stride)
        self.X, self.Y = np.meshgrid(x, y)
        
        # Configuration de la figure
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.ax.set_xlim(0, self.nx)
        self.ax.set_ylim(0, self.ny)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_aspect('equal')
        
        # Initialisation des éléments graphiques
        self.quiver = None
        self.colormap = None
        self.title_text = self.ax.text(0.02, 0.98, '', transform=self.ax.transAxes, 
                                      fontsize=14, verticalalignment='top',
                                      bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
    def animate(self, frame):
        """Fonction d'animation appelée à chaque frame"""
        # Nettoyer les éléments précédents
        if self.quiver:
            self.quiver.remove()
        if self.colormap:
            self.colormap.remove()
            
        # Extraire les données pour ce frame
        vx = self.Vx[frame, ::self.stride, ::self.stride]
        vy = self.Vy[frame, ::self.stride, ::self.stride]
        
        # Calculer la magnitude pour le fond coloré
        vx_full = self.Vx[frame, :, :]
        vy_full = self.Vy[frame, :, :]
        magnitude = np.sqrt(vx_full**2 + vy_full**2)
        
        # Fond coloré (magnitude)
        self.colormap = self.ax.imshow(magnitude, extent=[0, self.nx, 0, self.ny], 
                                      origin='lower', cmap='viridis', alpha=0.4,
                                      vmin=0, vmax=np.max(np.sqrt(self.Vx**2 + self.Vy**2)))
        
        # Champ de vecteurs
        self.quiver = self.ax.quiver(self.X, self.Y, vx, vy, 
                                    angles='xy', scale_units='xy', scale=self.scale,
                                    color='red', alpha=0.8, width=0.003)
        
        # Mise à jour du titre
        time = frame * self.dt
        self.title_text.set_text(f'Temps: {time:.2f} | Frame: {frame}/{self.nt-1}')
        
        return [self.quiver, self.colormap, self.title_text]
    
    def create_animation(self, interval=100, save_path=None, fps=10):
        """
        Créer et optionnellement sauvegarder l'animation
        
        Parameters:
        -----------
        interval : int
            Délai entre frames en millisecondes
        save_path : str, optional
            Chemin pour sauvegarder (ex: 'velocity_field.mp4')
        fps : int
            Images par seconde pour la sauvegarde
        """
        # Créer l'animation
        anim = FuncAnimation(self.fig, self.animate, frames=self.nt,
                           interval=interval, blit=False, repeat=True)
        
        # Ajouter une colorbar
        if not hasattr(self, 'colorbar_added'):
            # Frame temporaire pour créer la colorbar
            temp_magnitude = np.sqrt(self.Vx[0]**2 + self.Vy[0]**2)
            im = self.ax.imshow(temp_magnitude, extent=[0, self.nx, 0, self.ny], 
                               origin='lower', cmap='viridis', alpha=0)
            cbar = plt.colorbar(im, ax=self.ax, shrink=0.8)
            cbar.set_label('Magnitude de vitesse', rotation=270, labelpad=20)
            im.remove()
            self.colorbar_added = True
        
        # Sauvegarder si demandé
        if save_path:
            print(f"Sauvegarde en cours vers {save_path}...")
            try:
                if save_path.endswith('.gif'):
                    anim.save(save_path, writer='pillow', fps=fps)
                elif save_path.endswith('.mp4'):
                    # Essayer ffmpeg d'abord
                    try:
                        anim.save(save_path, writer='ffmpeg', fps=fps, bitrate=1800)
                    except:
                        # Fallback vers HTML5 si ffmpeg non disponible
                        print("ffmpeg non disponible, sauvegarde en HTML...")
                        html_path = save_path.replace('.mp4', '.html')
                        anim.save(html_path, writer='html', fps=fps)
                        print(f"Sauvegardé en HTML: {html_path}")
                        return anim
                else:
                    # Autres formats
                    anim.save(save_path, fps=fps)
                print(f"Animation sauvegardée : {save_path}")
            except Exception as e:
                print(f"Erreur de sauvegarde: {e}")
                print("Essayez le format GIF ou installez ffmpeg pour MP4")
                # Sauvegarde de secours en GIF
                gif_path = save_path.rsplit('.', 1)[0] + '.gif'
                anim.save(gif_path, writer='pillow', fps=fps)
                print(f"Sauvegardé en GIF à la place: {gif_path}")
        
        return anim

# Fonction utilitaire simple
def animate_velocity_field(Vx, Vy, save_path=None, fps=10, interval=100, **kwargs):
    """
    Fonction simple pour créer rapidement une animation
    
    Usage:
    ------
    anim = animate_velocity_field(Vx, Vy, save_path='video.mp4', 
                                 stride=3, scale=0.05, fps=15)
    plt.show()
    """
    animator = VelocityFieldAnimation(Vx, Vy, **kwargs)
    return animator.create_animation(save_path=save_path, fps=fps, interval=interval)

# Exemple d'utilisation avec données factices
def test_animation():
    """Test avec des données factices"""
    # Création de données test
    nt, ny, nx = 50, 30, 40
    t = np.linspace(0, 4*np.pi, nt)
    x = np.linspace(0, 2*np.pi, nx)
    y = np.linspace(0, 2*np.pi, ny)
    
    X, Y = np.meshgrid(x, y)
    Vx = np.zeros((nt, ny, nx))
    Vy = np.zeros((nt, ny, nx))
    
    for i, ti in enumerate(t):
        Vx[i] = np.sin(X + ti) * np.cos(Y)
        Vy[i] = np.cos(X + ti) * np.sin(Y)
    
    # Créer l'animation
    anim = animate_velocity_field(Vx, Vy, stride=2, scale=0.1)
    plt.show()
    return anim

# Pour utiliser avec vos vraies données :
"""
anim = animate_velocity_field(Vx, Vy, 
                             save_path='velocity_evolution.mp4',  # GIF plus fiable
                             stride=3, scale=5.0, fps=8)  # scale plus grand = flèches plus petites
plt.show()
"""
# %%
plt.figure()
plt.plot(np.cumsum(Vy,axis=0)[:,yind,xind])
plt.title('np.cumsum at a point on ice surface, W='+str(W)+', piv_indices:yind='+str(yind)+',xind='+str(xind))
plt.show()
# %% try to interpolate the piv field 
# to perform a Lagrangian-like tracking of the ice surface

from scipy.interpolate import RegularGridInterpolator


factor_interp = 10
def interpolate_piv_field(piv_field,factor_interp):

    xinterp = np.linspace(np.min(xpix),np.max(xpix),len(xpix)*factor_interp)
    yinterp = np.linspace(np.min(ypix),np.max(ypix),len(ypix)*factor_interp)

    interp = RegularGridInterpolator((xpix, ypix), piv_field.T)

    Xnew, Ynew = np.meshgrid(xinterp, yinterp)
    points = np.vstack((Xnew.ravel(), Ynew.ravel())).T

    # Interpolation vectorisée
    piv_field_interp = interp(points).reshape(len(yinterp), len(xinterp))

    return xinterp,yinterp,piv_field_interp

#exemple:
piv_field = Vy[0]
xinterp, yinterp, piv_field_interp = interpolate_piv_field(piv_field=piv_field,factor_interp=factor_interp)
plt.figure()
plt.imshow(piv_field_interp)
plt.show()


# ok on arrive à l'interpoler pour 1 frame


# %%
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

"""# Grille eulérienne
nt, ny, nx = 100, 50, 50
t = np.linspace(0, 10, nt)
y = np.linspace(0, 1, ny)
x = np.linspace(0, 1, nx)
dt = t[1] - t[0]

# Grille 3D (x, y, t) pour définir les champs
X, Y, T = np.meshgrid(x, y, t, indexing='xy')

# Champs de vitesse
# Vx : circulation vers la droite
V_field_x = np.cos(2 * np.pi * (Y + 0.2 * T)) * np.sin(2 * np.pi * X)
# Vy : oscillation verticale
V_field_y = np.sin(2 * np.pi * (Y - 0.2 * T)) * np.cos(2 * np.pi * X)

# Re-order pour (t, y, x)
V_field_x = V_field_x.transpose(2, 1, 0)
V_field_y = V_field_y.transpose(2, 1, 0)
"""

#
dt = 1
t = np.arange(0,Vy.shape[0],dt)
nt = len(t)
#x = np.arange(Vy.shape[0])
#t = np.arange(Vy.shape[0])

















"""
spectrum = np.fft.fft(Vy_converted_px,axis=0)

spectrum[int(spectrum.shape[0]/10):-int(spectrum.shape[0]/10) , : , :] = 0

Vy_filt = np.fft.ifft(spectrum,axis=0)

spectrum = np.fft.fft(Vx_converted_px,axis=0)

spectrum[int(spectrum.shape[0]/10):-int(spectrum.shape[0]/10) , :, :] = 0

Vx_filt = np.fft.ifft(spectrum,axis=0)
"""















# Interpolateurs
interp_Vx = RegularGridInterpolator((t, ypix, xpix), Vx_converted_px, bounds_error=False, fill_value=0)
interp_Vy = RegularGridInterpolator((t, ypix, xpix), Vy_converted_px, bounds_error=False, fill_value=0)



"""interp_Vx = RegularGridInterpolator((t, ypix, xpix), Vx_filt, bounds_error=False, fill_value=0)
interp_Vy = RegularGridInterpolator((t, ypix, xpix), Vy_filt, bounds_error=False, fill_value=0)
"""

# Particules initiales

"""n_particles = 10

x0 = np.linspace(0.1, 0.9, n_particles)
y0 = np.ones(n_particles) * 0.5

X_lag[0, :] = x0
Y_lag[0, :] = y0
"""

n_particles = Vy.shape[1]*Vy.shape[2]

X_lag = np.zeros((nt, n_particles))
Y_lag = np.zeros((nt, n_particles))

initial_positions = np.meshgrid(xpix,ypix)

X_lag[0, :], Y_lag[0, :] = initial_positions[0].ravel(),initial_positions[1].ravel()






from scipy.ndimage import gaussian_filter1d
parallel=False
smooth=True


#x_lag_along_profile = X_lag[:,len(xpix)*yind:len(xpix)*(yind+1)]
#y_lag_along_profile = Y_lag[:,len(xpix)*yind:len(xpix)*(yind+1)]



for n in range(1, nt):
    t_prev = t[n - 1]  # scalaire
    x_prev = X_lag[n - 1]  # (n_particles,)

    y_prev = Y_lag[n - 1]  # (n_particles,)
    # essayer en appliquant un filtre gaussien sur les profils à chaque pas de temps
    """
    x_prev_reshaped = np.reshape(x_prev,(len(ypix),len(xpix)))
    fft_xprev_reshaped = np.fft.fft(x_prev_reshaped,axis=1)
    fft_xprev_reshaped[:,3:-3] = 0
    x_prev_reshaped_filt = np.real(np.fft.ifft(fft_xprev_reshaped,axis=1))

    #x_prev_reshaped_smoothed = gaussian_filter1d(x_prev_reshaped,sigma=3,axis=1)

    y_prev_reshaped = np.reshape(y_prev,(len(ypix),len(xpix)))
    fft_yprev_reshaped = np.fft.fft(y_prev_reshaped,axis=1)
    fft_yprev_reshaped[:,3:-3] = 0
    y_prev_reshaped_filt = np.real(np.fft.ifft(fft_yprev_reshaped,axis=1))


    #y_prev_reshaped_smoothed = gaussian_filter1d(y_prev_reshaped,sigma=3,axis=1)

    x_prev = np.ravel(x_prev_reshaped_filt)
    y_prev = np.ravel(y_prev_reshaped_filt)
    
    """
    # Interpolateurs sur tout le batch
    # ATTENTION : certains interpolateurs comme scipy's RegularGridInterpolator nécessitent des inputs (N, 3)
    points = np.stack([np.full_like(x_prev, t_prev), y_prev, x_prev], axis=-1)  # shape: (n_particles, 3)
    vx = interp_Vx(points)
    vy = interp_Vy(points)

    X_lag[n] = x_prev + dt * vx
    Y_lag[n] = y_prev + dt * vy





plot_traj=True
if plot_traj:
    # Tracer les trajectoires
    plt.figure(figsize=(6, 6))
    for i in range(n_particles):
        plt.plot(X_lag[:700, i], Y_lag[:700, i])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Trajectoires lagrangiennes (Vx et Vy)')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    plt.show()


# %%
idx_profile = 5

# profil en lagrangien (évolution de la surface de la glace si la ref glace immobile est bien faite)




x_lag_along_profile = X_lag[:,len(xpix)*idx_profile:len(xpix)*(idx_profile+1)]
y_lag_along_profile = Y_lag[:,len(xpix)*idx_profile:len(xpix)*(idx_profile+1)]


plt.figure()
plt.plot(y_lag_along_profile[9,:])
plt.show()

# fracture à frame 700
#%%
plt.figure()
for i in range(500,510,1):
    plt.plot(y_lag_along_profile[i,:])
plt.show()

# pour l'instant rien de mieux que le champ de vitesses euleriennes...
# %%

plt.close('all')




