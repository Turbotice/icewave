# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:51:57 2023

@author: Banquise
"""

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import scipy.fft as fft
import os
import glob
from scipy.signal import correlate

import baptiste.files.save as sv
import baptiste.display.display_lib as disp
import baptiste.tools.tools as tools


def demodulation(t,s, fexc):
    c = np.nanmean(s*np.exp(1j * 2 * np.pi * t[None,None,:] * fexc),axis=2)
    return c

params = {}


#%% IMPORT FILES

date = '20230313'


folder = 'W:\Banquise/Rimouski_2023/Data/Geophones\\' +date+'\\'


# ordinateur perso
savefolder = 'W:\Banquise\Rimouski_2023\Traitements_donnees\\baptsarahantonin\\'



if not os.path.exists(savefolder):
    os.makedirs(savefolder)

filelist = glob.glob(folder+'*.txt')
print(filelist)
filename = filelist[3]

name = os.path.basename(filename).split('.')[0]

data = np.loadtxt(filename, skiprows = 22)

#%% PARAMETRES

params['facq'] = 1000
params['lm'] = 10000

params['fmax_analyse'] = 40

params['index_geophone'] = [[0, 3],[6, 0],[3, 6]]

size = data.shape[0]
params['index'] = 0



#%% MAIN


# for www in range (0,len(params['index_geophone'])) :
#     params['index'] = www

method = 1

t = np.linspace(0, params['lm'] / params['facq'], params['lm'])

params["phi"] = np.zeros((int(size/params['lm']), int((params['fmax_analyse']*params['lm']/params['facq'])+1)))

for i in range (0, int(size/params['lm'])) :
    if np.mod(i, 10) == 0 :
        print ('iteration ' + str(i) + ' sur ' + str( int(size/params['lm'])))
    
    data_cut = data[i * params['lm']:(i+1) * params['lm'], params['index_geophone'][params['index']]]
    
    Y = fft.fft(data_cut, axis = 0)
    
    for j in range (1, int((params['fmax_analyse']*params['lm']/params['facq'])+1)) :

        if method == 1 :
        
            Y1= np.zeros(Y.shape[0], dtype= 'complex64')
            Y2= np.zeros(Y.shape[0], dtype= 'complex64')
             
            Y1[j] = Y[j, 0]
            Y2[j] = Y[j, 1]
            
            Y1 = np.real(fft.ifft(Y1))
            Y2 = np.real(fft.ifft(Y2))
            
            corr = correlate(Y1,Y2)
            
            lag = np.argmax(corr)
            
            a = 5
            
            z = corr[lag-a:lag+a]
            
            x = [u for u in range (lag-a,lag+a) - lag ]
            
            p = np.polyfit(x,z,2)    #on peut utiliser fminsearch pour fitter par une fonction quelconque
            
            maxmax = -p[1]/(2*p[0])
            
            lag = lag + maxmax - a - len(corr)/2
            
            period = params['lm'] / j
                
            phi = lag / period * 2 * np.pi
            
            params['phi'][i,j] = phi
            

        if method == 2 :
            
            demod1 = demodulation(t, data_cut[:,0], j)
            demod2 = demodulation(t, data_cut[:,1], j)
            


#%% Phase shift et histogram
save = False

params = disp.figurejolie(params = params, nom_fig = 'phase_shift_f_ttmorceaux')
phase_shift = params['phi']

params['n_morc'] = np.shape(phase_shift)[0]
params['n_f'] = np.shape(phase_shift)[1]

while phase_shift[phase_shift> np.pi] != np.array([]) or phase_shift[phase_shift< -np.pi] != np.array([]) :
    phase_shift[phase_shift< -np.pi] += 2 * np.pi
    phase_shift[phase_shift> np.pi] += - 2 * np.pi

f = np.linspace(1, params['n_f'], params['n_f']) / params['lm'] * params['facq']
for u in range (0, params['n_morc']):
    plt.plot(f, phase_shift[u,:], '+')

if save :
    sv.save_graph(savefolder, 'phase_shift_f_ttmorceaux', params = params)



params = disp.figurejolie(params, nom_fig = 'phase_shift_f_morceaux')

n_phase = np.linspace(0, params['n_morc'], params['n_morc'])

params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, n_phase ,color = 5, table = np.flip(np.rot90(phase_shift), 0))

if save :
    sv.save_graph(savefolder, 'phase_f_morceau', params = params)


params['n_bins'] = 24
histogram = np.zeros(( params['n_bins'], params['n_f']))

params = disp.figurejolie(params = params, nom_fig = 'histogram')

for k in range (0,params['n_f']) :
    if np.mod(k, 100) == 0:
        print("iteration " + str(k) + " sur " + str(params['n_f']))
    
    h = plt.hist(phase_shift[:,k], params['n_bins'])
    
    h = np.asarray(h)
    
    h[1] = (h[1][1:]+h[1][:-1]) / 2
    
    histogram[:,k] = h[0]

params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé')
boites = np.linspace(0,params['n_bins'],params['n_bins'])
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, boites, table = np.flip(np.rot90(histogram), 0))
plt.clim(0,30)

if save :
    sv.save_graph(savefolder, 'phase_shift_histogram', params = params)

#%% Detect points
save = False

params['nb_blocks'] = 12

test_detect = np.flip(histogram, 0)
for i in range (0, params['nb_blocks'] - 1) :
    test_detect = np.concatenate((test_detect, np.flip(histogram, 0)), axis = 0)

params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé_fois4')

len_phase = np.linspace(-np.pi, 2 * np.pi * params['nb_blocks'] - np.pi, params['nb_blocks'] * params['n_bins'])
plt.pcolormesh(f,len_phase, test_detect)
cbar = plt.colorbar()
plt.clim(0,10)

Y = fft.fft(test_detect, axis = 0)

phase = - np.angle(Y[params['nb_blocks'],:])

phase = np.unwrap(phase) - np.pi + 2* np.pi

params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 2, exp = False)

if save :
    sv.save_graph(savefolder, 'phase_en_fonction_de_f_branches_plus_detection', params = params)

params = disp.figurejolie(params = params, nom_fig = 'phase_en_fonction_de_f')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 2, exp = False)

if save :
    sv.save_graph(savefolder, 'phase_en_fonction_de_f_branche_detectée', params = params)
    
    
#%% Selection data
select_data = True
cut_data = False

if select_data :
    params['liste_f_garde'] = [[0,22.3],[30.3,32]]
    params['liste_f_2pi'] = [[23,29.5]]
    params['liste_f_moins2pi'] = [[33.5,36]]
    
    f_new = []
    phase_new = []
    
    for i in range (0, len(params['liste_f_garde'])):
        i_f1 = int(params['liste_f_garde'][i][0] * params['lm'] / params['facq'])
        i_f2 = int(params['liste_f_garde'][i][1] * params['lm'] / params['facq'])
        f_new.extend(f[i_f1:i_f2])
        phase_new.extend(phase[i_f1:i_f2])
    
    for i in range (0, len(params['liste_f_2pi'])):
        i_f1 = int(params['liste_f_2pi'][i][0] * params['lm'] / params['facq'])
        i_f2 = int(params['liste_f_2pi'][i][1] * params['lm'] / params['facq'])
        f_new.extend(f[i_f1:i_f2])
        phase_new.extend(phase[i_f1:i_f2]+ 2 * np.pi)
        
    for i in range (0, len(params['liste_f_moins2pi'])):
        i_f1 = int(params['liste_f_moins2pi'][i][0] * params['lm'] / params['facq'])
        i_f2 = int(params['liste_f_moins2pi'][i][1] * params['lm'] / params['facq'])
        f_new.extend(f[i_f1:i_f2])
        phase_new.extend(phase[i_f1:i_f2] - 2 * np.pi)
    
    f_new,phase_new = tools.sort_listes(f_new,phase_new, reverse=False)

elif cut_data :
    params['f_max'] = 10       #frequence ou on coupe le signal en Hz
    params['f_min'] = 0       #frequence ou on coupe le signal en Hz
    
    params['f_max'] = int(params['f_max'] * params['lm'] / params['facq'])
    params['f_min'] = int(params['f_min'] * params['lm'] / params['facq'])
    
    f_new = f[params['f_min']:params['f_max']]
    phase_new = phase[params['f_min']:params['f_max']]
    
    

else :
    f_new = f
    phase_new = phase

   
#%% Unité physique
save = False

params['L'] = 20            #distance géophones en metres
params['theta'] = 11 / 180 * np.pi
if params['index'] == 0 :
    params['theta'] = params['theta'] + np.pi/6
    
if params['index'] == 1 :
    params['theta'] = - params['theta'] + np.pi/6
    
if params['index'] == 0 :
    params['theta'] = - params['theta'] + np.pi/2
    


if np.cos(params['theta']) <= 0 :
    params['theta'] = params['theta'] + np.pi

omega = f_new * 2 * np.pi
k = phase_new / params['L'] /np.cos(params['theta'])

params = disp.figurejolie(params = params, nom_fig = 'k_de_omega')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'k ($m^{-1}$)',r'$\omega$ (Hz)', k, omega, color = 4, exp = False)

if save :
    sv.save_graph(savefolder, 'k_de_omega', params = params)
    
#%% FIT
save = False
save_matlab = False

import baptiste.math.fits as fit
import baptiste.math.RDD as rdd

type_fit = 'pesante_flexion'


if type_fit == 'pesante_flexion':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion, k, omega, display = True, err = False, nb_param = 2, p0 = [0.3, 1e4], bounds = [0,[0.4,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['hxrho'] = popt[0]
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt[1]
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)

type_fit = 'power_law'
if type_fit == 'power_law':
    params = disp.figurejolie(params = params, nom_fig = 'k_de_omega')
    popt = fit.fit_powerlaw(k[:],omega[:], display = True, xlabel = '', ylabel = '', legend = '')
    params[str(params['num_fig'][-1])]['power'] = popt
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)

save = False
type_fit = 'pesante_flexion_depth_2m'
if True :#type_fit == 'pesante_flexion_depth':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion_depth, k, omega, display = True, err = False, nb_param = 2, p0 = [0.3, 1e4], bounds = [0,[0.4,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['hxrho'] = popt[0]
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt[1]
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
        
type_fit = 'flexion'      
if type_fit == 'flexion':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    # popt, pcov = fit.fit_ransac (k, omega, thresh = 0.1, display = True, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$', newfig = False)
    popt,pcov = fit.fit(rdd.RDD_flexion, k, omega, display = True, err = False, nb_param = 1, p0 = [1e4], bounds = [0,50e10], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
    
      
if save_matlab :
    sv.save_mat(k, savefolder, title = 'k_' + str(params['index']))
    sv.save_mat(k, savefolder, title = 'omega_' + str(params['index']))
    full_params = {}
    for i in params.keys() :
        if not i.isdigit() :
            full_params[i] = params[i]
    sv.save_mat(full_params, savefolder, title = 'parametres_' + str(params['index']))




#%% SAVE ALL

save = False
if save :
    sv.save_all_figs(savefolder, params)
    




























    
    

