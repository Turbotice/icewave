#%% imports modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import pickle
import os
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage
from vasco.tools.clickonfigures import get_n_points

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken
from profiles import *

#%% definition des chemins des données
#disk = 'L:'# disk is Elements on adour
disk = 'C:'
date = '0211'

#data_path = f'{disk}/Share_hublot/Data'
data_path = f'C:/Users/Vasco Zanchi/Desktop/Saguenay2024'

daily_drone_data_path = f'{data_path}/{date}/Drones'

velocity_field_path = f'{daily_drone_data_path}/exact_solution_real_field_stereo_0211_2024_rectangular_grid.h5'

# traitement stéphane : postitions des fractures (obtenue avec divergence du champ de vitesses)
fractures_positions_path = f'{daily_drone_data_path}/Results/fracture_positions.pkl'

#%% Cellule à exécuter une ceule fois (chargement champs vitesses from stereo piv)
#load dict from h5 file

dict_stereo_pivdata = rw.load_dict_from_h5(velocity_field_path)

#%% chargement des positions des fractures détectées
with open(fractures_positions_path, 'rb') as file:
    fractures_positions_data = pickle.load(file)

# %% definitions des variables utiles (champs d'élévation etc.)
dict_stereo_pivdata.keys()

vx = dict_stereo_pivdata['u'][0,:,:,:]
vy = dict_stereo_pivdata['u'][1,:,:,:]
vz = dict_stereo_pivdata['u'][2,:,:,:]

dt = dict_stereo_pivdata['t'][1] - dict_stereo_pivdata['t'][0]

ux = np.cumsum(vx, axis=2)*dt
uy = np.cumsum(vy, axis=2)*dt
uz = np.cumsum(vz, axis=2)*dt

facq_x = dict_stereo_pivdata['SCALE']['facq_x']

#%% chargement de dict_frac
%matplotlib inline

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'
if os.path.exists(path_dict2save):
    with open(path_dict2save, 'rb') as file:
        dict_frac = pickle.load(file)
image = vz[:,:,1000]
plt.imshow(image)
plt.imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')

if 'dict_frac' in globals():
    for k in dict_frac:
        if 'single_frac' in k:
            plt.plot(dict_frac[k]['idcs_single_frac'][0][0], dict_frac[k]['idcs_single_frac'][0][1], '^r', markersize=4)
plt.show()


#%% mesurer direct courbure (à moitié à la main car on doit cliquer)

"""
%matplotlib qt

def click2chooseprofile_measure_curvature(uz, index_time=1000, refpoints_optional=None):
    ''' methode 2 : faire 2 clics puis
      afficher le profil entre ces 2 clics, puis 
     faire un clic proche du max, puis fitter une parabole proche de ce max
    '''


    prof, (xstart,ystart), (xend,yend) = profile_line_on_image_2clicks(uz[:,:,index_time], refpoints_optional=refpoints_optional)

    x2plot = np.arange(len(prof))*(1/facq_x)
    coords = get_n_points(x2plot,prof)

    x_click = coords[0][0]
    idx_center = np.where(x2plot>=x_click)[0][0]
    indices2fit = np.arange(idx_center-5, idx_center+5)

    def poly2ndorder(x,a,b,c):
        return (a/2)*x**2 + b*x + c

    popt, pcov = curve_fit(poly2ndorder, x2plot[indices2fit], prof[indices2fit])

    %matplotlib qt
    plt.figure()
    plt.plot(x2plot, prof)
    plt.plot(x2plot[indices2fit], poly2ndorder(x2plot[indices2fit], popt[0], popt[1], popt[2]))
    plt.show()

    kappa_opt = popt[0]
    kappa_err = np.sqrt(np.diag(pcov)[0])

    return (kappa_opt, kappa_err), (prof, (xstart,ystart), (xend,yend), x_click, indices2fit, popt, pcov)

# maintenant appliquer cette méthode pour toutes les fractures qui sont dans dict_frac
for k in dict_frac:

    idcs_single_frac = dict_frac[k]['idcs_single_frac']
    time_frac_approx = dict_frac[k]['times_frac_sec_approx'][0]
    index_time_frac_approx = np.where(dict_stereo_pivdata['t']>=time_frac_approx)[0][0]

    print(index_time_frac_approx)

    (kappa_opt, kappa_err), other_infos = click2chooseprofile_measure_curvature(uz, index_time=index_time_frac_approx, refpoints_optional=idcs_single_frac)

    dict_frac[k]['measure_kappac_tip_approx'] = {}
    dict_frac[k]['measure_kappac_tip_approx']['kappa_opt'] = kappa_opt
    dict_frac[k]['measure_kappac_tip_approx']['kappa_err'] = kappa_err
    dict_frac[k]['measure_kappac_tip_approx']['other_infos'] = other_infos

"""

#%% simplement mesurer la longueur d'onde avec 2 clics fonctionne mieux...

def measure_wavelength_with2clicks(vz, facq_x, index_time=1000, vmin=-0.1, vmax=0.1, refpoints_optional=None):

    image = vz[:,:,index_time]
    coords = get_n_points_onimage(image, n_points=2, refpoints_optional=refpoints_optional, vmin=vmin, vmax=vmax)
    
    x1 = coords[0][0]
    y1 = coords[0][1]
    x2 = coords[1][0]
    y2 = coords[1][1]

    distance_px = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    distance = (1/facq_x) * distance_px

    return distance, coords

%matplotlib  qt

dict_frac['lambda_byhand'] = {}
dict_frac['lambda_byhand']['list_frac_id'] = []
dict_frac['lambda_byhand']['wavelength'] = []
dict_frac['lambda_byhand']['time_frac_approx'] = []
dict_frac['lambda_byhand']['coords_clicks'] = {}



for k in dict_frac:
    if 'dict_single_frac' in k:
        idcs_single_frac = dict_frac[k]['idcs_single_frac']
        time_frac_approx = dict_frac[k]['times_frac_sec_approx_ref_noncassee'][0]
        index_time_frac_approx = np.where(dict_stereo_pivdata['t']>=time_frac_approx)[0][0]

        print(index_time_frac_approx)
        wavelength, coords_clicks = measure_wavelength_with2clicks(vz, facq_x, index_time=index_time_frac_approx, refpoints_optional=idcs_single_frac)

        dict_frac['lambda_byhand']['list_frac_id'].append(k)
        dict_frac['lambda_byhand']['wavelength'].append(wavelength)
        dict_frac['lambda_byhand']['time_frac_approx'].append(time_frac_approx)
        dict_frac['lambda_byhand']['coords_clicks'][k] = coords_clicks
    else:
        pass

############## save dict #####################

path_dict2save = f'{daily_drone_data_path}/Results/traitement_vasco/dict_results_frac.pkl'

save_input = input('save ? (y/n)')
if save_input=='y':
    pickle.dump(dict_frac, open(path_dict2save, "wb"))
else:
    pass



#%%

from find_wave_front import *



for k in dict_frac:
    if 'dict_single_frac' in k:
        idcs_single_frac = dict_frac[k]['idcs_single_frac']
        time_frac_approx = dict_frac[k]['times_frac_sec_approx_ref_noncassee'][0]
        index_time_frac_approx = np.where(dict_stereo_pivdata['t']>=time_frac_approx)[0][0]

        matrice2d = uz[:,:,index_time_frac_approx]

        peaks2d, matrice2d_smoothed = find_wave_fronts_on_image(matrice2d, sigma_smooth=5, axis=0, plot=False)

        plt.figure()
        plt.imshow(matrice2d_smoothed)
        plt.imshow(peaks2d)
        plt.plot(idcs_single_frac[:,0], idcs_single_frac[:,1], '^k')
        plt.show()



#%% FFT2 in space to measure angle and approx wavelength (semble pas marcher)
"""
%matplotlib qt

def measure_wavelength_fourier(vz, index_time=1000, padding_factor=5, tolerance_px=5, plot_profiles=True):

    # mesure de l'angle de l'onde pour une image
    image = vz[:,:,index_time]

    image_croped = image[:,:60] # car sur la gauche l'angle est un peu different


    absfft2shifted = np.fft.fftshift(np.abs(np.fft.fft2(image_croped, s=(vz.shape[0]*padding_factor,vz.shape[1]*padding_factor))))

    # detection of maxima au voisinage d'un endroit où on a cliqué
    coords = get_n_points_onimage(absfft2shifted)
    #tolerance_px = 5
    xindclic = int(coords[0][0])
    yindclic = int(coords[0][1])


    maxval = np.max(absfft2shifted[yindclic-tolerance_px:yindclic+tolerance_px, xindclic-tolerance_px:xindclic+tolerance_px])
    npwh = np.where(absfft2shifted==maxval)
    kyindmax = npwh[0][0]
    kxindmax = npwh[1][0]
    print(npwh)

    kyzeroind = int(vz.shape[0]*padding_factor/2)
    kxzeroind = int(vz.shape[1]*padding_factor/2)

    diff_kxind = kxindmax - kxzeroind
    diff_kyind = kyindmax - kyzeroind 

    diff_kx = diff_kxind * (2*np.pi*facq_x/image_croped.shape[1])/padding_factor
    diff_ky = diff_kyind * (2*np.pi*facq_x/image_croped.shape[0])/padding_factor

    k_norm = np.sqrt(diff_kx**2+diff_ky**2)


    dominating_angle_wave_deg = 180/np.pi * np.arctan(diff_kyind/diff_kxind)
    print('angle wave =',dominating_angle_wave_deg,'deg')

    plt.figure()
    plt.imshow(absfft2shifted)
    plt.plot(kxindmax, kyindmax,'r+')
    plt.show()

    wavelength = 2*np.pi/k_norm
    print('wavelength = ',wavelength,' m')

    if plot_profiles:
        angle_deg = dominating_angle_wave_deg
        profiles = all_profiles(image, angle_deg)
        plot_profilelines_onimage(image, profiles)


    return dominating_angle_wave_deg, wavelength


measure_wavelength_fourier(vz,index_time=500)
measure_wavelength_fourier(vz,index_time=1000)
measure_wavelength_fourier(vz,index_time=1500)
measure_wavelength_fourier(vz,index_time=2000)
measure_wavelength_fourier(vz,index_time=2500)
measure_wavelength_fourier(vz,index_time=3000)
"""
#%matplotlib qt
"""
# faire un tableau de profil dans la direction de propagation de l'onde

"""

#%% test avec gradient : pas l'air de marcher
"""
%matplotlib qt



duzdy, duzdx = np.gradient(uz[:,:,index_time])

norm = np.sqrt(duzdx**2 + duzdy**2)

yarr = np.arange(uz.shape[0])
xarr = np.arange(uz.shape[1])

Xarr, Yarr = np.meshgrid(xarr,yarr)

plt.figure()
plt.imshow(uz[:,:,index_time])
plt.quiver(
    Xarr, Yarr, duzdx/norm, duzdy/norm,
    norm,                    # couleur selon ||∇f||
    cmap='plasma',
    scale=100,               # ajuster la taille des flèches
    width=0.001,
    pivot='mid',
    alpha=0.9
)
"""