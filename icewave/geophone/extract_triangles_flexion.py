# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:51:57 2023

@author: Banquise
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt2d
from scipy.signal import find_peaks
from skimage.measure import profile_line
import scipy.fft as fft
import scipy.io as io
import os
import glob
from scipy.signal import correlate

import baptiste.files.save as sv
import baptiste.display.display_lib as disp


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
# sur adour
#base = '/media/turbots/DATA/thiou/labshared2/'

filelist = glob.glob(folder+'*.txt')
print(filelist)
#230313_1714_048.txt'
filename = filelist[3]
#print(filename)

name = os.path.basename(filename).split('.')[0]

data = np.loadtxt(filename, skiprows = 22)

#%% PARAMETRES

params['facq'] = 1000
params['lm'] = 10000

params['fmax_analyse'] = 40

params['index_geophone'] = [[0, 3],[0, 6],[3, 6]]

size = data.shape[0]
params['index'] = 0



#%% MAIN

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
            
            peaks, u = find_peaks(corr)
            
            diff_peaks = peaks[1:] - peaks[:-1]
            
            period = np.median(diff_peaks)
            
            lag = np.argmax(corr)
            
            a = 5
            
            z = corr[lag-a:lag+a]
            
            x = [u for u in range (lag-a,lag+a) - lag ]
            
            p = np.polyfit(x,z,2)    #on peut utiliser fminsearch pour fitter par une fonction quelconque
            
            maxmax = -p[1]/(2*p[0])
            
            lag = lag + maxmax - a - len(corr)/2
                
            phi = lag / (period * 2 * np.pi)
            
            params['phi'][i,j] = phi
            

        
        if method == 2 :
            demod1 = demodulation(t, data_cut[:,0], j)
            demod2 = demodulation(t, data_cut[:,1], j)
            


#%% Phase shift et histogram

phase_shift = params['phi']
   
phase_shift[phase_shift< -np.pi] += 2 * np.pi
phase_shift[phase_shift> np.pi] += - 2 * np.pi

params = disp.figurejolie(params, nom_fig = 'phase_shift')
f = np.linspace(1, np.shape(phase_shift)[1], np.shape(phase_shift)[1]) / np.shape(data_cut)[0] * 1000
n_phase = np.linspace(0, np.shape(phase_shift)[0], np.shape(phase_shift)[0])

params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, n_phase ,color = 5, table = np.flip(np.rot90(phase_shift), 0))

# sv.save_graph(savefolder, 'phase_f_morceau', params = params)


params['n_bins'] = 24
histogram = np.zeros(( params['n_bins'], np.shape(phase_shift)[1]))

params = disp.figurejolie(params = params, nom_fig = 'histogram')

for k in range (0,np.shape(phase_shift)[1]) :
    print(k)
    
    h = plt.hist(phase_shift[:,k], params['n_bins'])
    
    h = np.asarray(h)
    
    h[1] = (h[1][1:]+h[1][:-1]) / 2
    
    histogram[:,k] = h[0]

params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé')
boites = np.linspace(0,params['n_bins'],params['n_bins'])
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, boites, table = np.flip(np.rot90(histogram), 0))

# sv.save_graph(savefolder, 'phase_shift_histogram', params = params)

#%% Detect points

test_detect = np.concatenate((np.flip(histogram, 0),np.flip(histogram, 0),np.flip(histogram, 0),np.flip(histogram, 0)), axis = 0)

params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé_fois4')

plt.pcolormesh(test_detect)
cbar = plt.colorbar()
plt.clim(0,30)

Y = fft.fft(test_detect, axis = 0)

phase = - np.angle(Y[1,:])

phase = np.unwrap(phase)+ 2 * np.pi

params = disp.figurejolie(params = params, nom_fig = 'phase_en_fonction_de_f')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 4, exp = False)
sv.save_graph(savefolder, 'phase_en_fonction_de_f', params = params)

#%% test

params = disp.figurejolie(params = params, nom_fig = 'phase_en_fonction_de_f')
Y = fft.fft(test_detect, axis = 0)

for i in range (0,Y.shape[0]):
    print(i)
    phi = - np.angle(Y[i,:])
    
    phi = np.unwrap(phi) 
    
    plt.plot(f, phi)
    
#%% Unité physique

phase[110:230] = phase[110:230] + 2 * np.pi


params['L'] = 20
params['theta'] = 0
params['f_coup'] = 159

omega = f[:params['f_coup']] * 2 * np.pi

k = phase[:params['f_coup']] / params['L'] /np.cos(params['theta'])

plt.figure()
plt.plot(k, omega)

#%% FIT

import baptiste.math.fits as fit
import baptiste.math.RDD as rdd

popt, pcov = fit.fit(rdd.RDD_flexion, k, omega, display = True, err = False, nb_param = 1, p0 = [6e6], bounds = False, zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')








































    
    

