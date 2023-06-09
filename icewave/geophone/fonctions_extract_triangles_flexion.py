# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 14:06:28 2023

@author: Banquise
"""


import os
import glob
import numpy as np

import baptiste.display.display_lib as disp
import baptiste.tools.tools as tools

import scipy.fft as fft
from scipy.signal import correlate
import matplotlib.pyplot as plt

def import_files (date,num_fichier, folder, savefolder):
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    
    filelist = glob.glob(folder+'*.txt')
    filename = filelist[num_fichier]
    data = np.loadtxt(filename,skiprows=22)
    print(filename)
    return data


def direction(params):
    
    if params['voie'] =='Z':
        params['index_geophone'] = [[0, 3],[3, 6],[6, 0]]
    elif params['voie'] == 'L':
        params['index_geophone'] = [[1, 4],[4, 7],[1, 7]]
    elif params['voie']=='T' :
        params['index_geophone'] = [[2, 5],[5, 8],[2, 8]]
    if params['num_geo'] == '4-6':
        params['index'] = 0
    elif params['num_geo'] == '6-12':
        params['index'] = 1
    elif params['num_geo'] == '12-4':
        params['index'] = 2
    return params

def demodulation(t,s, fexc):
    c = np.nanmean(s*np.exp(1j * 2 * np.pi * t[None,None,:] * fexc),axis=2)
    return c


def demod(y,f0,facq):
    n = len(y)
    dt = 1/facq
    t = np.arange(0,n*dt,dt)
    
    return np.mean(y*np.exp(1j * 2 * np.pi * t * f0),axis=0,dtype='complex64')

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def dephasage(y1,y2,Np=2**15,nf = 500):
    n =len(y1)
    N = int(np.floor(n/Np))

    R = np.zeros(nf+1)
    Theta = np.zeros(nf+1)

    freqs = np.linspace(0,1,nf+1)
    
    for j,f0 in enumerate(freqs):
        c1 = np.zeros(N,dtype='complex64')
        c2 = np.zeros(N,dtype='complex64')

        Z = np.zeros(N,dtype='complex64')
        for i in range(N):
            c1[i] = demod(y1[i*Np:(i+1)*Np],f0)
            c2[i] = demod(y2[i*Np:(i+1)*Np],f0)

            z = c1[i]*np.conj(c2[i])#/np.abs(c1[i]*np.conj(c2[i]))
            Z[i] = z

        Zm = np.mean(Z)

        [r,theta] = cart2pol(np.real(Zm),np.imag(Zm))

        R[j]=r
        Theta[j]=theta
    return freqs,R,Theta

def main(params,data, method=1):

    params = direction(params)
    size = np.shape(data)[0]
    params["phi"] = np.zeros((int(size/params['lm']), int((params['fmax_analyse']*params['lm']/params['facq']))))
    params["hist_weight"] = np.zeros((int(size/params['lm']), int((params['fmax_analyse']*params['lm']/params['facq']))))
    
    for i in range (0, int(size/params['lm'])) :
        if np.mod(i, 10) == 0 :
            print ('iteration ' + str(i) + ' sur ' + str( int(size/params['lm'])))
        
        data_cut = data[i * params['lm']:(i+1) * params['lm'], params['index_geophone'][params['index']]]
        
        Y = fft.fft(data_cut, axis = 0)
        
        for j in range (1, int((params['fmax_analyse']*params['lm']/params['facq']))+1) :
    
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
                p = np.polyfit(x,z,2)   
                maxmax = -p[1]/(2*p[0])
                lag = lag + maxmax - a - len(corr)/2
                period = params['lm'] / j    
                phi = lag / period * 2 * np.pi 
                params['phi'][i,j-1] = phi
                params['hist_weight'][i,j-1] = np.max(corr)
        
        
    return params

def histo(params, data, add_nom_fig = ''):
    params = main(params, data)
     
    params = disp.figurejolie(params = params, nom_fig = 'phase_shift_f_ttmorceaux_' + add_nom_fig)
    phase_shift = params['phi']
    
    params['n_morc'] = np.shape(phase_shift)[0]
    params['n_f'] = np.shape(phase_shift)[1]
    
    n_phase = np.linspace(0, params['n_morc'], params['n_morc'])
    f = np.linspace(1, params['n_f'], params['n_f']) / params['lm'] * params['facq']
    histogram = np.zeros(( params['n_bins'], params['n_f']))

    
    while True in (phase_shift> np.pi) or True in (phase_shift < -np.pi) :
        phase_shift[phase_shift< -np.pi] += 2 * np.pi
        phase_shift[phase_shift> np.pi] += - 2 * np.pi
    
    for u in range (0, params['n_morc']):
        plt.plot(f, phase_shift[u,:], '+')
        
    
    params = disp.figurejolie(params, nom_fig = 'phase_shift_f_morceaux_'+ add_nom_fig)
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, n_phase ,color = 5, table = np.flip(np.rot90(phase_shift), 0))
    
    
    params = disp.figurejolie(params = params, nom_fig = 'histogram_'+ add_nom_fig)
    for k in range (0,params['n_f']) :
        h = plt.hist(phase_shift[:,k], params['n_bins'], weights = params['hist_weight'][:,k], density = True)
        h = np.asarray(h)
        h[1] = (h[1][1:]+h[1][:-1]) / 2
        histogram[:,k] = h[0]
    
    
    params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé_'+ add_nom_fig)
    boites = np.linspace(0,params['n_bins'],params['n_bins'])
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, boites, table = np.flip(np.rot90(histogram), 0))
    plt.clim(0,np.quantile(histogram,0.9))
    params['phase_shift'] = phase_shift
    return params, f, histogram


def detect(params,f,histogram, add_nom_fig = ''):
    
    test_detect = np.flip(histogram, 0)
    for i in range (0, params['nb_blocks'] - 1) :
        test_detect = np.concatenate((test_detect, np.flip(histogram, 0)), axis = 0)
    
    params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé_fois4_'+ add_nom_fig)
    
    len_phase = np.linspace(-np.pi, 2 * np.pi * params['nb_blocks'] - np.pi, params['nb_blocks'] * params['n_bins'])
    plt.pcolormesh(f,len_phase, test_detect)
    plt.colorbar()
    plt.clim(0,np.quantile(test_detect, 0.9))
    
    Y = fft.fft(test_detect, axis = 0)
    
    phase = - np.angle(Y[params['nb_blocks'],:])
    phase = np.unwrap(phase) - np.pi + 2* np.pi
    
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 2, exp = False)
    
    
    params = disp.figurejolie(params = params, nom_fig = 'phase_en_fonction_de_f_'+ add_nom_fig)
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 2, exp = False)
    
    if params['select_data'] or params['cut_data'] :
        
        params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammé_fois4_select_data_'+ add_nom_fig)
        plt.pcolormesh(f,len_phase, test_detect)
        plt.colorbar()
        plt.clim(0,np.quantile(test_detect, 0.9))
        
        f_new, phase_new = select_data(params, f, phase)
        
        params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f_new, phase_new, color = 2, exp = False)
        
        return f_new, phase_new
    
    else :
        return f, phase

def select_data(params,f, phase) :

    f_new = []
    phase_new = []
    
    if params['select_data'] :
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
    
   

    elif params['cut_data'] :
        
        f_max = int(params['f_max'] * params['lm'] / params['facq'])
        f_min = int(params['f_min'] * params['lm'] / params['facq'])
        
        f_new = f[f_min:f_max]
        phase_new = phase[f_min:f_max]
        
   
    
    else :
        f_new = f
        phase_new = phase
    
    return f_new, phase_new

def omega_k(params, f, phase):
    omega = f * 2 * np.pi
    k = phase / params['L'] /np.cos(params['theta_f'])
    return omega, k
