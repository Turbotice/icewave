#%%
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
import sys
from scipy.io import loadmat
sys.path.append('../../post_traitement/piv_grenoble/relation_dispersion/functions')
from disprel import load_complex_field, extract_data_mat_file

sys.path.append('../../../manips_Bs_As/tracking_force_displacement/functions')
from utils import load_data_PIV_without_coord, reshape_array, clean_array

#%%

def extract_amplitude(matdata_path,ind_yprof,ind_xinf,ind_xsup, plot=True):
    u, v, _ = load_data_PIV_without_coord(matdata_path)

    u = np.float64(u)
    v = np.float64(v)
    if plot:
        plt.figure()
        for i in range(v.shape[0]):
            plt.plot(np.arange(ind_xinf,ind_xsup),v[i,ind_yprof,ind_xinf:ind_xsup]) # maximum correspond à 2A, où A est l'amplitude de l'élévation
        plt.show()

    # on fait un histogramme des maximas
    Maximas = np.max(np.abs(v[:,ind_yprof,ind_xinf:ind_xsup]), axis=0)
    hist, bin_edges = np.histogram(Maximas, bins=5, density=False)
    idxmax = np.argmax(hist)
    valmax = (bin_edges[idxmax+1]+bin_edges[idxmax])/2
    #print(valmax)
    valmax_err = (bin_edges[1]-bin_edges[0])/np.sqrt(12)
    #print(valmax_err)
    if plot:
        plt.figure()
        plt.hist(Maximas, bins=5, density=False)
        plt.errorbar([valmax],[np.max(hist)],xerr=[valmax_err],ecolor='r',marker='o',color='r')
        plt.title('Histogram of peak to peak amplitude of the elevation, in pixels')
        plt.show()

    # on en déduit l'amplitude correspondante

    amplitude = valmax/2
    amplitude_err = valmax_err/2
    #print('Amplitude=', amplitude, '+-', amplitude_err, 'px')

    return amplitude,amplitude_err

# fichier d'exemple
#matdata_path = r"R:\Gre25\Data\0527\cameras\manip_relation_dispersion\Acquisition_3\camera_40307970\80Hz_157.001Hz\matData\PIV_processed_i00_N0_Dt50_b1_W64_full.mat"

# test
#extract_amplitude(matdata_path,2)

# %%
