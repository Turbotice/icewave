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
import baptiste.files.save as sv

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

def main(params,data):

    params = direction(params)
    size = np.shape(data)[0]
    
    
    params['lm'] = np.array([50000])
    params['liste_f'] = np.array([0.06])
    
    while params['liste_f'][-1] < params['fmax_analyse']:
        if params['lm'][-1] > 3000 :
            params['lm'] = np.append(params['lm'], params['lm'][-1] / 1.1)
            params['liste_f'] = np.append(params['liste_f'], 3 * params['facq'] / params['lm'][-1])
        else :
            params['lm'] = np.append(params['lm'], params['lm'][-1])
            params['liste_f'] = np.append(params['liste_f'], params['liste_f'][-1] + 2 * params['facq'] / params['lm'][-1])
        
    params["phi"] = np.zeros(len(params['liste_f']), dtype = object)
    params["hist_weight"] = np.zeros(len(params['liste_f']), dtype = object)
    
    
    for j in range (len(params['liste_f'])) :
        
        if np.mod(j, 10) == 0 :
            print ('iteration ' + str(j) + ' sur ' + str(len(params['liste_f'])) )
        
        # params['corrs'] = np.zeros( (int(params['fmax_analyse']*params['lm']/params['facq']),2 * Y.shape[0] - 1) )
        
        phi = np.array([])
        corrs = np.array([])
        
        for i in range (1, int(size/params['lm'][j])) :
            
            
            
            data_cut = data[i * int(params['lm'][j]): (i+1) * int(params['lm'][j]), params['index_geophone'][params['index']]]
                                                                      
            Y = fft.fft(data_cut, axis = 0)
            
            Y1= np.zeros(Y.shape[0], dtype= 'complex64')
            Y2= np.zeros(Y.shape[0], dtype= 'complex64')
             
            Y1[int(params['liste_f'][j] * params['facq'] * params['lm'][j] / size)] = Y[int(params['liste_f'][j] * params['facq'] * params['lm'][j] / size), 0]
            Y2[int(params['liste_f'][j] * params['facq'] * params['lm'][j] / size)] = Y[int(params['liste_f'][j] * params['facq'] * params['lm'][j] / size), 1]
            
            Y1 = np.real(fft.ifft(Y1))
            Y2 = np.real(fft.ifft(Y2))
            
            corr = correlate(Y1,Y2)
            # params['corrs'][j-1, :] = corr
            lag = np.argmax(corr)
            a = 5
            z = corr[lag-a:lag+a]
            x = [u for u in range (lag-a,lag+a) - lag ]
            p = np.polyfit(x,z,2)   
            maxmax = -p[1]/(2*p[0])
            lag = lag + maxmax - a - len(corr)/2
            period = params['lm'][j] / i   
            phi = np.append(phi, lag / period * 2 * np.pi)
            corrs = np.append(corrs, np.max(corr))
            
        params['phi'][j] = phi
        params['hist_weight'][j] = corrs
        
    return params




def histo(params, data, add_nom_fig = '', full_display = True):
    
    params = main(params, data)
    size = np.shape(data)[0]
    
    if full_display :
        params = disp.figurejolie(params = params, nom_fig = 'phase_shift_f_ttmorceaux_' + add_nom_fig)
     
    #On interpolle phi pour pouvoir l'afficher
    params['og_phi'] = params['phi'].copy()
    
    params['phi'] = []
    

    for i in range (np.shape(params['og_phi'])[0]):
        x = np.linspace(0,int(size/params['lm'][i])-2, int(size/params['lm'][-1]) - 1)
        x_old = np.linspace(0,int(size/params['lm'][i] - 2) , int(size/params['lm'][i]) - 1)
        params['phi'].append(np.interp(x, x_old, params['og_phi'][i]))
                                 
                                 
    params['phi'] = np.asarray(params['phi'])
    
    params['n_morc'] = np.shape(params['phi'])[1]
    params['n_f'] = np.shape(params['liste_f'])[0]
    
    n_phase = np.linspace(0, params['n_morc'], params['n_morc'])
    f = params['liste_f']
    histogram = np.zeros(( params['n_bins'], params['n_f']))

    while True in (params['phi'] > 2*np.pi) or True in (params['phi'] < 0) :
        params['phi'][params['phi'] < 0] += 2 * np.pi
        params['phi'][params['phi'] > 2*np.pi] += - 2 * np.pi
    
    if full_display :
        for u in range (0, params['n_morc']):
            plt.plot(f, params['phi'][u,:], '+')
        
    if full_display :
        params = disp.figurejolie(params, nom_fig = 'phase_shift_f_morceaux_'+ add_nom_fig)
        params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','n_phase', f, n_phase ,color = 5, table = np.flip(np.rot90(params['phi']), 0))
    
    
    params = disp.figurejolie(params = params, nom_fig = 'histogram_'+ add_nom_fig)
    for k in range (0,params['n_f']) :
        if params['weight'] :
            h = plt.hist(params['phi'][:,k], params['n_bins'], weights = params['hist_weight'][:,k], density = True)
        else :
            h = plt.hist(params['phi'][:,k], params['n_bins'])
        h = np.asarray(h)
        h[1] = (h[1][1:]+h[1][:-1]) / 2
        histogram[:,k] = h[0]
    
    x = np.linspace(0, 2 * np.pi, params['n_bins'] - 1)
    disp.set_axe_pi(3, x)
    
    if full_display :
        params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammed_'+ add_nom_fig)
        boites = np.linspace(0,params['n_bins'],params['n_bins'])
        params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f(Hz)',r'$\phi$', f, boites, table = np.flip(np.rot90(histogram), 0))
        
        plt.clim(0,np.quantile(histogram,0.9))
        x = np.linspace(0, params['n_bins'], params['n_bins'] - 1)
        disp.set_axe_pi(3, x, axxe = 'y')

    return params, f, histogram


def detect(params,f,histogram, add_nom_fig = '', save = False):
    
    test_detect = np.flip(histogram, 0)
    for i in range (0, params['nb_blocks'] - 1) :
            test_detect = np.concatenate((test_detect, np.flip(histogram, 0)), axis = 0)
    
    params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammed_fois_'+ add_nom_fig)
    
    len_phase = np.linspace(0, 2 * np.pi * params['nb_blocks'], params['nb_blocks'] * params['n_bins'])
    plt.pcolormesh(f,len_phase, test_detect)
    plt.colorbar()
    plt.clim(0,np.quantile(test_detect, 0.9))
    
    x = len_phase
    disp.set_axe_pi(params['nb_blocks'] * 2 + 1, x, axxe = 'y')
    
    Y = fft.fft(test_detect, axis = 0)
    
    phase = - np.angle(Y[params['nb_blocks'],:])
    phase = np.unwrap(phase)
    
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f(Hz)', r'$\phi$', f, phase, color = 2, exp = False)
    
    if params['select_data'] == False and params['cut_data'] == False and save :

        sv.save_graph (params["savefolder"], nom_fig = 'phase_shift_histogrammed_fois_'+ add_nom_fig, params = params)
        
    params = disp.figurejolie(params = params, nom_fig = 'phase_en_fonction_de_f_'+ add_nom_fig)
    params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f, phase, color = 2, exp = False)
    
    if params['select_data'] or params['cut_data'] :
        
        params = disp.figurejolie(params = params, nom_fig = 'phase_shift_histogrammed_fois_select_data_'+ add_nom_fig)
        plt.pcolormesh(f,len_phase, test_detect)
        plt.colorbar()
        plt.clim(0,np.quantile(test_detect, 0.9))
        
        f_new, phase_new = select_data(params, f, phase)
        
        params[str(params['num_fig'][-1])]['data'] = disp.joliplot( 'f','phase', f_new, phase_new, color = 2, exp = False)
        
        if save :
            sv.save_graph (params["savefolder"], nom_fig = 'phase_shift_histogrammed_fois_select_data_'+ add_nom_fig, params = params)
        
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

# def select_data(params,f, phase) :

#     f_new = []
#     phase_new = []
    
#     if params['select_data'] :
#         for i in range (0, len(params['liste_f_garde'])):
#             i_f1 = np.where(f == params['liste_f_garde'][i][0] )
#             i_f2 = np.where(f == params['liste_f_garde'][i][1] )
#             f_new.extend(f[i_f1:i_f2])
#             phase_new.extend(phase[i_f1:i_f2])
        
#         for i in range (0, len(params['liste_f_2pi'])):
#             i_f1 = np.where(f == params['liste_f_2pi'][i][0] )
#             i_f2 = np.where(f == params['liste_f_2pi'][i][1] )
#             f_new.extend(f[i_f1:i_f2])
#             phase_new.extend(phase[i_f1:i_f2]+ 2 * np.pi)
            
#         for i in range (0, len(params['liste_f_moins2pi'])):
#             i_f1 = np.where(f == params['liste_f_moins2pi'][i][0] )
#             i_f2 = np.where(f == params['liste_f_moins2pi'][i][1] )
#             f_new.extend(f[i_f1:i_f2])
#             phase_new.extend(phase[i_f1:i_f2] - 2 * np.pi)
#         f_new,phase_new = tools.sort_listes(f_new,phase_new, reverse=False)
    
   

#     elif params['cut_data'] :
        
#         f_max = int(params['f_max'] * params['lm'] / params['facq'])
#         f_min = int(params['f_min'] * params['lm'] / params['facq'])
        
#         f_new = f[f_min:f_max]
#         phase_new = phase[f_min:f_max]
   
    
#     else :
#         f_new = f
#         phase_new = phase
    
#     return f_new, phase_new

def omega_k(params, f, phase):
    omega = f * 2 * np.pi
    k = phase / params['L'] /np.cos(params['theta_f'])
    return omega, k
