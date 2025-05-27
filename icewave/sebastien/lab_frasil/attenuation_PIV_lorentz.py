# -*- coding: utf-8 -*-
"""
Created on Mon May 26 21:34:49 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable

import h5py
import pickle
import os
import glob

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

def get_data(file2load):
    """ Load PIV data from .mat file
    Input : file2load, str, name of .mat file where data from PIV are saved
    Output : data, python dictionnary obtained using mat2py module """
    
    with h5py.File(file2load,'r') as fmat:
        data = mat2py.mat_to_dict(fmat['m'],fmat['m'])

    data = mat2py.transpose_PIVmat_fields(data)
    return data

def show_velocity_field(V,x,y,colormap,figname):

    field = V.T
    
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots(figsize = (12,9))
    c = ax.imshow(field, cmap = parula_map , aspect = 'equal', norm = 'linear', origin = 'lower',interpolation = 'gaussian',
                  extent = (x.min(),x.max(),y.min(),y.max()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$V_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)
    ax.set_xlabel(r'$x \; \mathrm{(m)}$', labelpad = 5)
    ax.set_ylabel(r'$y \; \mathrm{(m)}$', labelpad = 5)
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
def save_TF_spectrum(TF_spectrum,freq,figname):

    fig, ax = plt.subplots()
    set_graphs.set_matplotlib_param('single')
    ax.plot(freq,TF_spectrum)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-5,1e-2])
    
    ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 15)
    ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle _{x,y}(f) \; \mathrm{(u.a.)}$',labelpad = 5)
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    return 
    

def plot_FFT2D_and_max_peak(shift,kx,ky,peak_coords,limits,figname):
    """ Plot FFT 2D and maximum peak 
    Inputs : - shift, numpy.array 2D, contain 2D FFT [nkx,nky]
             - kx, numpy.array, values of wavevector x-coordinate
             - ky, numpy.array, values of wavevector y-coordinate 
             - peak_coords, tuple or array (2,), coordinates of peak (wavevector units)
             - limits, array, limits used by ax.set_xlim and ax.set_ylim : [xmin, xmax, ymin, ymax]
             - figname, str, name under which figure will be saved 
             """
    
    xmin = limits[0]
    xmax = limits[1]
    ymin = limits[2]
    ymax = limits[3]
    
    fig, ax = plt.subplots(figsize = (12,9))
    c = ax.imshow(abs(shift).T,cmap = parula_map,aspect = 'equal',norm = 'linear',origin = 'lower',
                  extent = (kx.min(),kx.max(),ky.min(),ky.max()))
    
    ax.plot(peak_coords[0],peak_coords[1],'ro')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\hat{V_x} (x,y) \; \mathrm{(u.a.)}$',labelpad = 1)
    ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    ax.set_ylabel(r'$k_y \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
def interpolate_zeros(cut,kx):
    """ Set data around kx = 0 to zeros, than replace these values by interpolated values from extracted profile 
    Inputs : - cut, np.array, 1D, profile of FFT 
             - kx, np.array, wavevector array (N,), kx = 0 is obtained at kx[N//2] 
    Outputs : - cut_full, np.array,1D, profile of FFT with interpolated values for cut[idx0 - 1 : idx0 + 2] 
              - kx_full, np.array, 1D, array of interpolated wavevector """

    idx0 = len(kx)//2 # index at which kx = 0
    new_cut = cut
    new_cut[idx0 - 1 : idx0 + 2] = 0
    
    first_part = cut[:idx0 - 1]
    second_part = cut[idx0 + 2:]
    cut_model = np.concatenate([first_part,second_part])
    
    first_part = kx[:idx0 - 1]
    second_part = kx[idx0 + 2:]
    kx_model = np.concatenate([first_part,second_part])
    
    # interpolate
    I = scipy.interpolate.interp1d(kx_model,cut_model,kind = 'linear')
    
    k2interp = kx[idx0 - 1 : idx0 + 2]
    kx_full = np.concatenate([first_part,k2interp,second_part])
    cut_full = I(kx_full)
    
    return cut_full,kx_full

def lorentzian_fit(cut,k,k0_idx,bounds_alpha):
    """ Perform Lorentzian fit along a profile. 
    Inputs : - cut, 1D numpy array, profile along which lorentzian fit must be performed
             - k, 1D numpy array or list, list of wavevectors
             - k0_idx, index of initial guess for center of lorentzian
             - bounds_alpha, boundaries between which lorentzian width must be checked """
    
    y_exp = (cut - cut.min())/(cut.max() - cut.min())

    bounds_kx0 = (k[k0_idx - 15],k[k0_idx + 15])
    bounds_curvefit = ([bounds_kx0[0],bounds_alpha[0]],[bounds_kx0[1],bounds_alpha[1]])
    # fit by a lorentzian
    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),k,y_exp,
                                         bounds = bounds_curvefit)
    err_coeff = np.sqrt(np.diag(pcov))
    
    return popt,err_coeff

def save_lorentzian_plot(cut,k,popt,figname):
    """ Save lorentzian plot """
    
    label_fit = r'$\alpha = ' + f'{popt[1]:.1f}' + '\; \mathrm{(m^{-1})}$'
    print(label_fit)
    k_fit = np.linspace(k.min(),k.max(),1000)
    yth = lorentzian(k_fit,popt[0],popt[1])*(cut.max() - cut.min()) + cut.min()
    
    fig, ax = plt.subplots()
    ax.plot(k,cut,'o-')
    ax.plot(k_fit,yth,'r',label = label_fit)
    ax.legend()
    
    ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
    ax.set_ylabel(r'$|\hat{V_x}|(k_x) \; \mathrm{(u.a.)}$',labelpad = 5)
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    plt.close(fig)
    return 


def process_folder(folder):
      
    file2load = glob.glob(f'{folder}/*scaled.mat')[0]
    
    # Define fig_folder and collect variables
    current_path = f'{folder}/'
    fig_folder = f'{current_path}Plots/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
    
    data = get_data(file2load)
    
    f_ex = data['EXP']['f_ex']
    h = data['EXP']['h']
    ID = data['ID']
    ID_txt = f'h_{h}mm_{ID}'.replace('.','p')
    print(data.keys())
    
    # Plot histogram of velocity 
    figname = f'{fig_folder}Histogram_PIV_close_wavemaker_{ID_txt}'
    FT.histogram_PIV(data['Vx']/data['SCALE']['scale_V'],data['PIV_param']['w'],figname)
    
    # Show velocity field 
    frame = 0
    V = data['Vx'][:,:,frame]
    figname = f'{fig_folder}Vx_frame{frame}_{ID_txt}'
    show_velocity_field(V,data['x'],data['y'],parula_map,figname)
    
    # Perform time FFT 
    TF_spectrum,freq,FFT_t = FT.temporal_FFT(data['Vx'],data['SCALE']['facq_t'],padding_bool = 1,add_pow2 = 1,output_FFT = True)
    figname = f'{fig_folder}TF_spectrum_{ID_txt}'
    save_TF_spectrum(TF_spectrum,freq,figname)
    
    # Show demodulated field 
    # find peak
    p = np.argmax(TF_spectrum)
    f_max = freq[p]
    print(f'Detected f_ex = {f_max:.2f}')
    
    demod_field = np.mean(data['Vx']*np.exp(1j*2*np.pi*f_max*data['t']),axis = -1)
    real_field = np.real(demod_field)
    
    figname = f'{fig_folder}Demodulated_Vx_{ID_txt}'
    show_velocity_field(real_field, data['x'], data['y'], parula_map, figname)
    
    # Perform FFT 2D of demodulated field
    add_pow2 = [1,1]
    facq = (1/data['SCALE']['fx'],1/data['SCALE']['fx']) # acquisition frequency for each direction x and y 
    shift,kx,ky = FT.fft_2D(demod_field,facq,add_pow2)
    
    # find maximum of the 2D FFT
    idx_max = np.argmax(abs(shift).flatten())
    unravel_coords = np.unravel_index(idx_max,shift.shape)
    peak_coords = [kx[unravel_coords[0]],ky[unravel_coords[1]]]
    limits = [-100,100,-100,100]
    figname = f'{fig_folder}Spatial_FFT_{ID_txt}'
    plot_FFT2D_and_max_peak(shift,kx,ky,peak_coords,limits,figname)
    
    # extract data along kx  
    cut = abs(shift)[:,unravel_coords[1]]
    cut,kx = interpolate_zeros(cut,kx)
    
    # fit by a lorentzian
    bounds_alpha = (0,1e2)
    popt,err_coeff = lorentzian_fit(cut, kx, unravel_coords[1], bounds_alpha)
    # save Lorentzian fit plot
    figname = f'{fig_folder}Lorentzian_fit_kx_{ID_txt}'
    save_lorentzian_plot(cut, kx, popt, figname)
    
    # save results 
    dict_res = {}
    dict_res['f_demod'] = f_max
    dict_res['k0'] = abs(popt[0])
    dict_res['err_k0'] = abs(err_coeff[0])
    dict_res['alpha'] = popt[1]
    dict_res['err_alpha'] = abs(err_coeff[1])
    dict_res['lorentz'] = {}
    dict_res['lorentz']['kx'] = kx
    dict_res['lorentz']['y'] = cut
    dict_res['PIV_param'] = data['PIV_param']
    dict_res['h'] = float(data['EXP']['h'])
    dict_res['f_ex'] = float(data['EXP']['f_ex'])
    dict_res['amplitude'] = float(data['EXP']['amplitude'])
    dict_res['ID'] = ID
    
    file2save = f'{current_path}PIV_attenuation_results_{ID_txt}.pkl'
    with open(file2save,'wb') as pf:
        pickle.dump(dict_res,pf)
        
    return 


def main(h = 7.5,date = '2024_07_11'):
    path2data = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/matData/'
    folderlist = glob.glob(f'{path2data}*')
    
    with ProcessPoolExecutor(max_workers=2) as executor:
        list(tqdm(executor.map(process_folder, folderlist), total=len(folderlist)))
        
    # for folder in folderlist :
    #     process_folder(folder)
        
    return 
    
if __name__ == '__main__':
    main()



