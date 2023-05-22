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
        print ('iteration ' + i)
    
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
#%% test

u = np.argmax(corr)

y0 = u

a = 5

z = corr[y0-a:y0+a]

x = [u for u in range (y0-a,y0+a) - y0 ]

p = np.polyfit(x,z,2)    #on peut utiliser fminsearch pour fitter par une fonction quelconque

maxmax = -p[1]/(2*p[0])

laglag = u + maxmax - a


plt.figure()

plt.plot(corr)
plt.plot(u, np.max(corr), 'ko')





#%% Phase shift

phase_shift = params['phi']

# if phase_shift < -np.pi :
#     phase_shift += 2 *np.pi
    
phase_shift[phase_shift< -np.pi] += 2 * np.pi
phase_shift[phase_shift> np.pi] += - 2 * np.pi

# figure(1)
# plot([1:size(phase_shift,2)]/size(cut,1)*1000,phase_shift,'+')

# for j=1:size(phase_shift,2)
#     j
#     figure(2)
#     h=histogram(phase_shift(:,j),24);
#     value(:,j)=h.Values;
# end
# close

plt.figure()
plt.pcolormesh(phase_shift - params['phi'])





#%% test 

i = 0
method = 1
print (i)

data_cut = data[i * params['lm']:(i+1) * params['lm'], params['index_geophone'][params['index']]]

Y = fft.fft(data_cut, axis = 0)

f= np.linspace(0, params['facq'], params['lm'])

plt.figure()
plt.plot(f, np.abs(Y[:,0]))
plt.plot(f, np.abs(Y[:,1]))

for j in range (1, int((params['fmax_analyse']*params['lm']/params['facq'])+1)) :
    print(j)
    if method == 1 :
    
        Y1= np.zeros(Y.shape[0], dtype= 'complex64')
        Y2= np.zeros(Y.shape[0], dtype= 'complex64')
         
        Y1[j] = Y[j, 0]
        Y2[j] = Y[j, 1]
        
        Y1 = np.real(fft.ifft(Y1))
        Y2 = np.real(fft.ifft(Y2))
        
        corr = correlate(Y1,Y2)
        
        peaks = find_peaks(corr)
        
        diff_peaks = peaks[:-1] - peaks[1:]
        
        period = np.median(diff_peaks)
        
        lag = np.argmax(corr)
        
        y0 = lag
        
        a = 5
        
        z = corr[y0-a:y0+a]
        
        x = [u for u in range (y0-a,y0+a) - y0 ]
        
        p = np.polyfit(x,z,2)    #on peut utiliser fminsearch pour fitter par une fonction quelconque
        
        lag = -p[1]/(2*p[0]) - len(corr)/2
            
        phi = lag / (period * 2 * np.pi)
        
        params['phi'][i,j] = phi
        

    
    if method == 2 :
        demod1 = demodulation(t, data_cut[:,0], j)
        demod2 = demodulation(t, data_cut[:,1], j)
    
    
    
    # sig1 = fft.ifft(np.append(demod1, np.zeros(1000)))
    # sig2 = fft.ifft(np.append(demod2, np.zeros(1000)))
    


    
    
    

    
    
    

    
    

