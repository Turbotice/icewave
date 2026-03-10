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
from PIL import Image
from scipy.signal import correlate2d
#%%

%matplotlib qt

# %% load data
W = 32
#Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 400
N = 0
refimg = 400

frame_frac = 457
frame_force_application = 437 #début de l'application de la force

dcm = 8
dpx = 313 # à changer en fonction de la manip
L = 8e-2 # distance entre 2 points d'appui

date = '20250516'
time = '135503'

system_loc = 'windows_server'
disk = 'D:/Grenoble/'
general_folder = f'{disk}/Gre25/Data/PIV_results/3points_bending/{date}_{time}/'


"""
if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/PIV_results/3points_bending/{date}_{time}/'
"""

path2data = general_folder

matfile = f'{path2data}PIV_processed_displacement_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

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


xpix = mat_dict['x'][0][0][0]
ypix = mat_dict['y'][0][0][:,0]



# %%
# la vis est à peu près au pixel : xpix=645 (indice 3) et ypix=607 (indice 3)

yind = 3
xind = 3

plt.figure()
plt.plot(np.arange(v.shape[0])+i0,v[:,yind,xind])
plt.xlim(i0,frame_frac)
plt.ylim(-np.max(v[:frame_frac-i0,yind,xind])/4,np.max(v[:frame_frac-i0,yind,xind]))
plt.show()




plt.figure()
plt.title('vertical displacement (cm) vs time (frames)')
plt.plot(np.arange(v.shape[0])+i0,v[:,yind,xind] * dcm/dpx)
plt.xlim(i0,frame_frac)
plt.ylim(-(dcm/dpx)*np.max(v[:frame_frac-i0,yind,xind])/4,(dcm/dpx)*np.max(v[:frame_frac-i0,yind,xind]))
plt.show()



h_measurements = np.array([2.16,1.98,1.84,1.58,1.33,1.68,1.85,2.01,2.12,2.21])
h_avg_mm = np.mean(h_measurements)
h_std_mm = np.std(h_measurements)

h_avg = 1e-3 * h_avg_mm
h_std = 1e-3 * h_std_mm
print('h='+str(h_avg)+'+-'+str(h_std))


def compute_kappa(delta,L):
    return 8*delta/(L**2)

def compute_epsilon(kappa,h):
    return kappa*h/2
# calcul approx pour odg :
kappa_c = compute_kappa(2.5e-4,8e-2)
epsilon_c = compute_epsilon(kappa_c,h_avg)

print('courbure critique : ',kappa_c,'m^-1')
print('déformation critique : ',epsilon_c)





#%%

def importer_images_grayscale(dossier_images, extension=".tiff", taille=None):
    """
    Importe une séquence d'images, les convertit en noir et blanc, et les stocke dans un tableau NumPy.

    :param dossier_images: Chemin du dossier contenant les images
    :param extension: Extension des fichiers image (par exemple ".png" ou ".jpg")
    :param taille: Tuple (largeur, hauteur) pour redimensionner les images (facultatif)
    :return: Tableau NumPy de forme (n_images, hauteur, largeur)
    """
    fichiers = sorted([f for f in os.listdir(dossier_images) if f.endswith(extension)])
    images = []

    for fichier in fichiers:
        chemin = os.path.join(dossier_images, fichier)
        img = Image.open(chemin).convert("L")  # "L" pour grayscale (niveaux de gris)

        if taille:
            img = img.resize(taille)

        img_array = np.array(img)
        images.append(img_array)

    return np.stack(images)

# Exemple d'utilisation
chemin_dossier = "D:/Grenoble/Gre25/Data/Readmi_01/Sample_fracture/VID_20250516_135503_images"
tableau = importer_images_grayscale(chemin_dossier)

print("Forme du tableau :", tableau.shape)  # (n_frames, hauteur, largeur)


#%% cross correlation pour mesurer la force appliquée pour fléchir la plaque

plt.figure()
plt.imshow(tableau[410,:,:])
plt.show()

# on va prendre un petit carré dans l'image, à corréler avec toute l'image pour chaque step pour tracker sa position
img2corr = tableau[410,1424:1442,611:625]

tab_cross_corr = []

for i in range(tableau.shape[0]):
    cross_corr = correlate2d(tableau[i,:,:],img2corr)
    tab_cross_corr.append(cross_corr)

tab_cross_corr = np.array(tab_cross_corr)
#%%
plt.figure()
plt.imshow(tab_cross_corr[420,:,:],vmin=0)
plt.colorbar()
plt.show()

