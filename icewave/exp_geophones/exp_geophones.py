# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 12:08:48 2023

@author: Banquise
"""

#%% INTIALISATION

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt2d
from scipy.signal import find_peaks
from skimage.measure import profile_line
import scipy.fft as fft
import scipy.io as io
import os
import cv2
import h5py

import baptiste.display.display_lib as disp
import baptiste.experiments.import_params as ip
import baptiste.signal_processing.fft_tools as ft
import baptiste.image_processing.image_processing as imp
import baptiste.math.fits as fits
import baptiste.math.RDD as rdd
import baptiste.files.save as sv
import baptiste.tools.tools as tools
import baptiste.signal_processing.smooth_signal as smooth
import baptiste.math.fits as fit

import icewave.geophone.fonctions_extract_triangles_flexion as fct

import mat73 as m7

date = '230516'
nom_exp = 'UQARG'

dico, params, loc = ip.initialisation(date)

loc_h = 'D:\Banquise\Baptiste\Resultats_video\d221104\d221104_PIVA6_PIV_44sur026_facq151Hz_texp5000us_Tmot010_Vmot410_Hw12cm_tacq020s/'
loc_h = 'W:\Banquise\\Data_DD_UQAR\\'
loc_h = 'W:\\Banquise\\Baptiste\\baie_HAHA_drone\\'
loc_h = 'D:\\Banquise\\Baptiste\\Drone\\Traitements\\'

loc_h = 'D:\\Banquise\\Baptiste\\Resultats_video\\d231013\\exp_analogue_geo\\test_2_13_10_23\\image_sequence_og\\'
loc_h = 'D:\\Banquise\\Baptiste\\Resultats_video\\d231013\\exp_analogue_geo\\test_2_13_10_23\\traitement\\sorted\\'

params['special_folder'] = 'Traitement_exp_geo'

path = loc_h

params['mmparpixel'] = 0.3990 #test sanalogue banquise
params['facq'] = 300

def make_triangle(baricentre, size, inclination, display = True, bg_image = False):

    r = size * (0.865 - 0.577/2)
    theta = np.array([-inclination + 60, -inclination + 180, -inclination + 300])
    x = r * np.cos(np.radians(theta)) + baricentre[0]
    y = r * np.sin(np.radians(theta)) + baricentre[1]
    
    if display :
        plt.figure()
        if type(bg_image) != bool :
            plt.imshow(bg_image)
            
        plt.plot(x, y, 'ro')
        plt.grid(False)  # Afficher une grille
        plt.title("Triangle équilatéral avec côté de {} et angle de {} degrés".format(size, inclination))
        plt.axis('equal')

    return x, y

params['PIV_file'] = 'Height_in_micron_surf_100001to105001.mat'


#%% Anciens triangles et leur importation

# params['PIV_file'] = 'Height_in_micron_surf_100001to105001.mat'
# data = m7.loadmat(loc_h + params['PIV_file'])
# #GRAND
# pts_1 = data['H'][11,443,:]      #G4
# pts_2 = data['H'][304,149,:]    #G6
# pts_3 = data['H'][410,550,:]    #G12


# #MOYEN
# pts_1 = data['H'][217,498,:]      #G4
# pts_2 = data['H'][359,357,:]    #G6
# pts_3 = data['H'][410,550,:]    #G12

# #PETIT
# pts_1 = data['H'][313,524,:]      #G4
# pts_2 = data['H'][388,453,:]    #G6
# pts_3 = data['H'][410,550,:]    #G12


triangle = 'grand'

grand_triangle0 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_1_10001.txt')
grand_triangle1 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_10001_20001.txt')
grand_triangle2 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_20001_30001.txt')
grand_triangle3 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_30001_40001.txt')
grand_triangle4 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_40001_50001.txt')
grand_triangle5 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_50001_60001.txt')
grand_triangle6 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_60001_70001.txt')
grand_triangle7 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_70001_80001.txt')
grand_triangle8 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_80001_90001.txt')
grand_triangle9 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_90001_100001.txt')
grand_triangle10 = np.loadtxt(loc_h + 'Triangle_' + triangle + '_100001_105001.txt')

# data_tot = np.concatenate((grand_triangle1, grand_triangle2, grand_triangle3))

data_tot = np.concatenate((grand_triangle0, grand_triangle1, grand_triangle2,grand_triangle3, grand_triangle4, grand_triangle5,grand_triangle6, grand_triangle7, grand_triangle8,grand_triangle9, grand_triangle10))

params['triangle'] = triangle

#%% Créer et afficher un triangle

import matplotlib.pyplot as plt
import numpy as np
import os
import cv2

loc_img = 'D:\\Banquise\\Baptiste\\Resultats_video\\d231013\\exp_analogue_geo\\test_2_13_10_23\\image_sequence_og\\'

liste_img = os.listdir(loc_img)

u = cv2.imread(loc_img + liste_img[0])

x,y = make_triangle([275,250], 380, 15, bg_image = u)
plt.grid('off')


#%% Import several triangles

f = h5py.File(loc_h + params['PIV_file'], 'r')

grand_triangle = []
grand_triangle.append([loc_h + 'Height_in_micron_surf_1to10001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_10001to20001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_20001to30001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_30001to40001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_40001to50001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_50001to60001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_60001to70001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_70001to80001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_80001to90001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_90001to100001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_100001to105001.mat'])


liste_inc = [-30,-15,0,15,30]
liste_size = [20,50,100,200,400]

params['liste_size'] = liste_inc #taille triangle en pixels
params['liste_inclinaison'] = liste_size #angle inclinaison triangle

data_tot = np.zeros((len(liste_inc)*len(liste_size)))

i = 0
j = 0
for inclination in liste_inc :
    print("loop " + str(j + i * len(liste_size)) + " sur " + str(len(liste_inc)* len(liste_size)))
    
    for size in liste_size :
        x,y = make_triangle([275,250], size, inclination, display = False)
        
        for aaa in range (11):

            G4 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[2]),int(y[2])]    #G4
            G6 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[1]),int(y[1])]    #G6
            G12 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[0]),int(y[0])]    #G12
                  
            data_triangle = np.vstack((G4, G4))
            data_triangle = np.vstack((data_triangle, G4))
            data_triangle = np.vstack((data_triangle, G6))
            data_triangle = np.vstack((data_triangle, G6))
            data_triangle = np.vstack((data_triangle, G6))
            data_triangle = np.vstack((data_triangle, G12))
            data_triangle = np.vstack((data_triangle, G12))
            data_triangle = np.vstack((data_triangle, G12))
            
            # Met dans le bon sens
            data_triangle = np.rot90(data_triangle)
            data_triangle = np.flip(data_triangle, 0)
            
            if aaa == 0:
                data_tot0 = data_triangle
            else :
                data_tot0 = np.concatenate((data_tot0,data_triangle) )
            print('data ' + str(aaa) + ' sur 11')
        
        data_tot[j + i * len(liste_size)] = data_tot0
        
        
        j += 1
    i += 1

#%% Import one Triangle

params['size'] = 200 #taille triangle en pixels
params['inclinaison'] = 15 #angle inclinaison triangle

x,y = make_triangle([275,250], params['size'], params['inclinaison'], display = False)


grand_triangle = []
grand_triangle.append([loc_h + 'Height_in_micron_surf_1to10001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_10001to20001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_20001to30001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_30001to40001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_40001to50001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_50001to60001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_60001to70001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_70001to80001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_80001to90001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_90001to100001.mat'])
grand_triangle.append([loc_h + 'Height_in_micron_surf_100001to105001.mat'])


for aaa in range (11):

    G4 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[2]),int(y[2])]    #G4
    G6 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[1]),int(y[1])]    #G6
    G12 = h5py.File(grand_triangle[aaa][0], 'r')['H'][:,int(x[0]),int(y[0])]    #G12
    
    
    data_triangle = np.vstack((G4, G4))
    data_triangle = np.vstack((data_triangle, G4))
    data_triangle = np.vstack((data_triangle, G6))
    data_triangle = np.vstack((data_triangle, G6))
    data_triangle = np.vstack((data_triangle, G6))
    data_triangle = np.vstack((data_triangle, G12))
    data_triangle = np.vstack((data_triangle, G12))
    data_triangle = np.vstack((data_triangle, G12))
    
    # Met dans le bon sens
    data_triangle = np.rot90(data_triangle)
    data_triangle = np.flip(data_triangle, 0)
    
    if aaa == 0:
        data_tot0 = data_triangle
    else :
        data_tot0 = np.concatenate((data_tot0,data_triangle) )
    print(aaa)

data_tot = data_tot0

for i in range (data_tot.shape[1]):
    data_tot[:,i] = data_tot[:,i] - np.mean(data_tot[:,i])

#%%PARAMS


params['facq'] = 300
params['lm'] = 120
params['fmax_analyse'] = 50


params['L'] = params['size'] * params['mmparpixel'] /1000
params['theta_f'] = 0

params['n_bins'] = 24
params['nb_blocks'] = 2
params['weight'] = True

params['savefolder'] = loc_h + '\\results\\'

#%% 4-6

params['voie']='Z'
params['num_geo'] = '4-6'

params['liste_f_garde'] = [[0,100], [26.7,30], [79.5, 100]] #[[0,22.2]] #17-14 Z
params['liste_f_2pi'] = [[40,70]] #17-14 Z
params['liste_f_moins2pi'] = [[100,1100]]#[[33.5,36]] #17-14 Z
params['f_min'] = 0
params['f_max'] = 12
params['select_data'] = False
params['cut_data'] = False

params,f1,histogram1=fct.histo(params,data_tot, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f1, phase1 = fct.detect(params, f1, histogram1, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = False)

phase1 = phase1 + 2 * np.pi

#%% 6-12

params['voie']='Z'
params['num_geo'] = '6-12'

params['liste_f_garde'] = [[4.7,17]]
params['liste_f_2pi'] = [[0,4.7]]
params['liste_f_moins2pi'] = [[17,36]]
params['f_min'] = 0
params['f_max'] = 12
params['select_data'] = False
params['cut_data'] = False

params,f2,histogram2=fct.histo(params,data_tot, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f2, phase2 = fct.detect(params, f2, histogram2, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = False)

phase2 = phase2 + 2 * np.pi

#%% 12-4

params['voie']='Z'
params['num_geo'] = '12-4'

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 0
params['f_max'] = 90
params['select_data'] = False
params['cut_data'] = False

params,f3,histogram3=fct.histo(params,data_tot, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f3, phase3 = fct.detect(params, f3, histogram3, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = False)

#%% FIND ANGLE
save = False


'''Créé un tableau avec les fréquence communes des deux géophones'''

f1_new, f2_new, phase1_new, phase2_new = tools.commun_space(f1, f2, phase1, phase2, pas_f = params['facq']/params['lm'])

# f1_new, f3_new, phase1_new, phase3_new = tools.commun_space(f1, f3, phase1, phase3, pas_f = params['facq']/params['lm'])

#affiche les deux voies (verifier que phase tend vers 0 pour f tend vers 0)
params = disp.figurejolie(params, nom_fig = 'Phase4_6etphase6_12')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f1_new, phase1_new, color = 2, exp = False, legend = 'G4-G6')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f2_new, phase2_new, color = 3, exp = False, legend = 'G6-G12')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f3_new, phase3_new, color = 4, exp = False, legend = 'G12-G4')

#formule pour theta12 en fonction des phases (4-6/6-12)
theta_f = np.arctan( -1 * (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

#formule pour theta12 en fonction des phases (4-6/12-4)
# theta_f = np.arctan( (1/np.sqrt(3)) * 2 * ( phase3_new/phase1_new + 0.5))

#formule pour 23 en fonction des phases (6-12/12-4)
# theta_f = np.arctan( (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

params = disp.figurejolie(params, nom_fig = 'theta_f_rawdata')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 5, legend = r'$\theta (f)$ raw')

theta_f = smooth.savgol(theta_f, 5, 2)
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 2, exp = False, legend = r'$\theta (f)$ smoothed')
# if save :
#     sv.save_graph (params["savefolder"], nom_fig = 'theta_f_rawdata', params = params)
# #on smooth
# theta_f = smooth.savgol(theta_f, 50, 1)

# params = disp.figurejolie(params, nom_fig = 'theta_f_smooth')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 5)

if save :
    sv.save_graph (params["savefolder"], nom_fig = 'theta_f_smooth', params = params)

params['theta_f'] = theta_f

#On trouve le omega et k
omega, k = fct.omega_k(params, f1_new, phase1_new)

params = disp.figurejolie(params, nom_fig = 'omega(k)')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot(r'k ($m^{-1}$)', r'$\omega$ (Hz)', k, omega, color = 2, exp = True, legend = 'RDD')

#%% FIT FLEXION  

save = False
type_fit = 'flexion'      
if type_fit == 'flexion':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt,pcov = fit.fit(rdd.RDD_flexion, k, omega, display = True, err = False, nb_param = 1, p0 = [2e-4], bounds = [0,1e6], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(params['savefolder'], type_fit, params = params)
        









#%% Traitement RDD classique

params['PIV_file'] = 'Height_in_micron_surf_100001to105001.mat'
data_totales = m7.loadmat(loc_h + params['PIV_file'])

data = data_totales['H'][:,:500,:500]


kacqx = 2 * np.pi / params['mmparpixel']
kacqy = 2 * np.pi / params['mmparpixel']

[ny,nx,nt] = data.shape

save_demod = False

fmin = 0.1 #2/nt * params['facq'] #fréquence minimale résolue pour avoir au moins 2 périodes par fréquence
fmax = 50
nb_f = 10
padding = [10,10]    #puissance de 2 pour le 0 padding
k_xx = []
k_yy = []
kk = []
theta = []
fff = []
t = np.linspace(0,nt / params['facq'], nt)

# mparpixel = 0.07988064791133845 * ratio_PIV #mesuré avec la taille du bateau...


plotplot = True
nb_plot = 5

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

huhu = 0
for i in np.linspace(fmin, fmax, nb_f) :
    huhu += 1
    if np.mod(huhu-1,1)==0:
        print('iteration ' + str(huhu) + ' sur ' + str(nb_f))
    demod = ft.demodulation(t,data,i)
    demod_padding = ft.add_padding(demod, padding)
    Y_FFT, k_x, k_y = ft.fft_bapt(demod_padding - np.mean(demod_padding), kacqx, kacqy)
    nx_FFT = np.shape(Y_FFT)[0]
    ny_FFT = np.shape(Y_FFT)[1]
    max_fft, maxx = ft.max_fft(Y_FFT, f1 = k_x,f2 = k_y, display = False)
    k_xx.append(max_fft[0])
    k_yy.append(max_fft[1])
    kk.append(cart2pol(max_fft[0],max_fft[1])[0])
    theta.append(max_fft[1]/max_fft[0])
    fff.append(i)
    


    if save_demod :
        sv.save_mat(demod, loc_h, title = "champ_demod_f_" + str(i))
    
    if plotplot :
        if int(huhu/nb_f * nb_plot ) == (huhu/nb_f * nb_plot):
            disp.figurejolie()
            ft.plot_fft(Y_FFT, k_x, k_y, tcbar = r'Démodulation à f = ' + str(round(i, 3)) + " Hz")
            plt.plot(max_fft[0], max_fft[1], 'ro')
 
k_xx = np.asarray(k_xx)
k_yy = np.asarray(k_yy)
kk = np.asarray(kk)
theta = np.asarray(theta)
fff = np.asarray(fff)

#%% Affichage/save RDD classique

save = False

params = disp.figurejolie(params, nom_fig = 'theta_de_f') 
params[str(params['num_fig'][-1])]['data'] = disp.joliplot(r"f (Hz)",r"$\theta$", fff, theta, color = 1)
 
if save :
    sv.save_graph (path, params[str(params['num_fig'][-1])]['nom_fig'], params = params)


params = disp.figurejolie(params, nom_fig = 'k_de_f') 
params[str(params['num_fig'][-1])]['data'] = disp.joliplot(r"f (Hz)",r"K (m$^{-1}$)", fff, kk, color = 2)
 
if save :
    sv.save_graph (path, params[str(params['num_fig'][-1])]['nom_fig'], params = params)

#%% Comparaison
#%% COrrélation

lencor = 105000

plt.figure()
plt.plot(data_tot[:500,0])
plt.plot(data_tot[:500,3])

plt.figure()

from scipy.signal import correlate
corr = correlate(data_tot[:lencor,0], data_tot[:lencor,3])
x = np.linspace(-lencor/params['facq'],lencor/params['facq'], 2 * lencor - 1 )

plt.plot(x, corr)
