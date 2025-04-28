# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 13:36:48 2025

@author: sebas
"""

import numpy as np
import os
import glob  
import h5py
import pickle
import scipy
import math
import cv2 as cv
# from skimage.filters import meijering, sato, frangi, hessian

import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
from matplotlib.path import Path
import matplotlib.animation as animation

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT
import icewave.tools.interactive_scatter_filter as scatter_filter
import icewave.tools.datafolders as df

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = mcolors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION 

def affine(x,a,b):
    y = a*x + b
    return y

def powerlaw_fit(x,y,err_y):
    """ Fit data using a power law, taking into account standard deviation of y """
    log_x = np.log(x)
    log_y = np.log(y)
    err_log_y = err_y/y
    
    popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(log_x,a,b),log_x,log_y,sigma = err_log_y,absolute_sigma = True)
    err_affine = np.sqrt(np.diag(pcov))
    beta = popt[0]
    err_beta = err_affine[0]
    B = np.exp(popt[1])
    err_B = B*err_affine[1]
    
    coeffs = (beta,B)
    err_coeffs = (err_beta,err_B)
    return coeffs,err_coeffs

def bound_harmonicN(k,N,h_w):
    """ Compute waves omegaN associated to bound wave of order N"""
    
    omegaN = np.sqrt(N*9.81*k*np.tanh(h_w*k/N))
    return omegaN

def unnormed_gaussian(x,x0,alpha):
    y = np.exp(-0.5*((x - x0)/alpha)**2)
    return y

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

def phase_velocity_shallow(k,hw,g = 9.81):
    """ Compute phase velocity based on shallow water equation"""
    c_phi = np.sqrt(g*np.tanh(k*hw)/k)
    return c_phi

def group_velocity_shallow(k,hw,g = 9.81):
    """ Compute group velocity based on shallow water equation """
    c_phi = phase_velocity_shallow(k,hw)
    c_g =0.5*c_phi*(1 + 2*k*hw/np.sinh(2*k*hw))
    
    return c_g


def time2space_attenuation(time_att,k,water_depth,error_bar = 1):
    """ Compute space attenuation from time attenuation using shallow water dispersion equation
    Inputs : - time_att, numpy.ndarray : N x 2 array, column 0 : time attenuation coefficient,
    column 1 : error on this coefficient
             - k, wavevector array
             - water_depth, tuple or array of size (2,), first value : water depth, second value : error on water depth
             - error_bar, boolean : default is 1, computes error_bar or not
    Outputs : - alpha, numpy.ndarray, spatial attenuation coefficients
              - err_alpha, numpy.ndarray, error on spatial attenuation coefficient
    """
    
    lambda_time = time_att[:,0]
    hw = water_depth[0]
    c_g = group_velocity_shallow(k,hw)
    alpha = lambda_time/c_g
    
    # compute errorbars
    if error_bar :
        print('Computing error bars, be sure that errorbars are second column of time_att and water_depth')
        err_lambda = time_att[:,1]
        err_hw = water_depth[1]
        err_khw = k*err_hw
        c_phi = phase_velocity_shallow(k, hw)
        error_c_phi = c_phi*err_khw/np.sinh(2*k*hw)
        error_c_g = np.sqrt((0.5*np.sinh(2*k*hw) + k*hw)**2 * error_c_phi**2 + 
                            c_phi**2 * (1 - 2*k*hw/np.tanh(2*k*hw))**2 * err_khw**2)/np.sinh(2*k*hw)
        err_alpha = np.sqrt(err_lambda**2 + lambda_time**2 * (error_c_g/c_g)**2)/c_g   
    
    return alpha,err_alpha


def animation_demodulated_field(selected_FFT,freq,colormap,cmax,time_interval = 1e3):
    
    """ Create animation of real demodulated fields for successive frequencies 
    Inputs : - selected_FFT, numpy array, [nx,ny,nf] Time Fourier transform of a 2D field
             - freq, numpy array, 1D, array of frequencies
             - colormap, cmap object, colormap used to plot real fields
             - cmax, float, value used to scale colorbar from -cmax to +cmax
             
    Outputs : - ani, matplotlib animation which can be saved in a .mp4 format
        """
    nb_frames = len(freq) # total number of frames of animation
    fig, ax = plt.subplots(figsize = (12,9))
    
    selected_freq = freq[0]
    print(f'frequency = {selected_freq}')
    
    field = selected_FFT[0]
    real_field = np.real(field)
    # show initial matrix 
    c = ax.imshow(real_field.T,cmap = colormap,aspect = 'equal', origin = 'lower', interpolation = 'gaussian',
              vmin = -cmax,vmax = cmax,extent = (data['x'].min(),data['x'].max(),data['y'].min(),data['y'].max()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$|\hat{V}_x| (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)
    
    ax.set_xlabel(r'$x \; \mathrm{(m)}$', labelpad = 5)
    ax.set_ylabel(r'$y \; \mathrm{(m)}$', labelpad = 5)
    
    def update(frame) : 
        selected_freq = freq[frame]
        print(f'frequency = {selected_freq}')
        
        field = selected_FFT[:,:,frame]
        real_field = np.real(field)
        c.set_data(real_field.T)
        ax.set_title('$f = '+ f'{selected_freq:.3f}'+ '\; \mathrm{(Hz)}$')
        
        return c
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update, frames=nb_frames, interval=time_interval)
    # plt.show()
    print('Animation computed')
    return ani




def indices2fit(y,x,p,rel_height):
    """ Compute indices of points to be fitted, based on index of a peak p and the y-values to fit
    Inputs : y, numpy.ndarray, values to fit
             x, numpy.ndarray, abscisse coordinates
             p, index of the detected peak, y[p] is therefore a local maximum
             rel_height, float between 0 and 1, relative prominence at which peak width should be computed if left or right limits
             are too large
    Outputs : points2fit, indices of point to be fitted"""

    # compute left/right limits
    decreasing_curve = 1
    left_limit = p
    right_limit = p
    while decreasing_curve :
        next_left = left_limit-1
        next_right = right_limit+1
        if np.logical_or(y[next_left] >= y[left_limit],
                         y[next_right] >= y[right_limit]) : 
            decreasing_curve = 0
        else : 
            left_limit = next_left
            right_limit = next_right

        if np.logical_or(left_limit == 0,right_limit == len(x)-1):
            decreasing_curve = 0
            widths = scipy.signal.peak_widths(y,[p],rel_height = rel_height)
            left_limit = math.floor(widths[2][0])-1
            right_limit = math.ceil(widths[3][0])+1
            
    # compute left/right distance to peak
    left_ips = x[p] - x[left_limit]
    right_ips = x[right_limit] - x[p] 
    # min_ips = min(left_ips,right_ips)
    points2fit = np.where(np.logical_and(x > x[p] - left_ips,x < x[p] + right_ips))[0]
    # points2fit = np.where(np.logical_and(Afk['k'] > Afk['k'][p] - min_ips,Afk['k'] < Afk['k'][p] + min_ips))[0]
    
    return points2fit


def matrix_peak_correlation_detection(M,x,y,x_range,y_range,model_signal,detec_param):
    """ Perform peak detection using gaussian correlation 
    Inputs : - M, numpy.ndarray [ny,nx]: matrix over which we perform peak detection
             - x, array type (nx,) : detection is performed at fixed x
             - y, array type (ny,): correlation is performed with respect to y coordinate
             - x_range,y_range, array (2,) : range of x and y over which we look for maxima
             - model_signal, array type (ny,): reference array used to perform correlation
             - detec_param, dict : parameters dictionnary to be used for find_peaks """
    
    delta_y = np.diff(y)[0]
    
    mask_y = np.where(np.logical_and(y < y_range[1],y > y_range[0]))[0]
    y = y[mask_y]
    M = M[mask_y,:]

    idx_min = np.argmin(np.abs(x_range[0] - x))
    idx_max = np.argmin(np.abs(x_range[1] - x))
    indices = np.arange(idx_min,idx_max,step = 1)

    # parameters for find peaks
    min_prominence = detec_param['prominence']
    rel_height = detec_param['rel_height']
    
    S = {'idx':[],'peaks':[],'A':[],'y':[],'x':[]}
    # create dictionnary in which results are saved 
    
    fig, ax = plt.subplots()
    for i in range(len(indices)):
        ax.clear()
        
        idx = indices[i]
        
        detected_x = x[idx]
        print(detected_x)
        section = M[:,idx]
        label_x = f'x = {detected_x:.2f}'
        
        norm = section/section.max()
        
        ax.plot(y,norm,'-o',label = label_x)
        x_model = np.arange(0,len(model_signal),step = 1)*delta_y
        ax.plot(x_model,model_signal,'-k')
        
        # Apply gaussian filter 
        correlated = scipy.signal.fftconvolve(norm,model_signal,mode = 'same')
        norm_correlated = correlated/correlated.max()
        ax.plot(y,norm_correlated)
        # find peaks of correlation
        peaks,properties = scipy.signal.find_peaks(norm_correlated,prominence = min_prominence)
        prominences = scipy.signal.peak_prominences(norm_correlated,peaks)
        widths = scipy.signal.peak_widths(norm_correlated,peaks,rel_height = rel_height)
        ax.vlines(x=y[peaks], ymin=norm_correlated[peaks] - prominences[0],
                ymax = norm_correlated[peaks],color = 'g')
        ax.hlines(y=widths[1], xmin=widths[2]*delta_y,
                xmax=widths[3]*delta_y,color = 'g')
        
        ax.plot(y[peaks],norm_correlated[peaks],'d')
        
        ax.set_xlabel(r'$y \; \mathrm{(Hz)}$')
        ax.set_ylabel(r'$A_{norm}$')
        ax.legend()
        
        plt.pause(1.0)
        # save parameters in dictionnary
        for i,p in enumerate(peaks) :
            S['idx'].append(idx)
            S['peaks'].append(p)
            S['A'].append(norm[p]*section.max())
            S['x'].append(detected_x)
            S['y'].append(y[p])
        
        detected_y = y[p]
        print(f'y = {detected_y}')
        # f_txt = f'{detected_f:.2f}'.replace('.','p')
        # figname = f'{fig_folder}Section_f{f_txt}_corr_gaussian_sigma{gaussian_width}'
        # plt.savefig(figname + '.pdf', bbox_inches='tight')
        # plt.savefig(figname + '.png', bbox_inches='tight')
        
    for key in S.keys():
        S[key] = np.array(S[key])

    return S

#%% Load data for a given date and experiment
date = '0226'
drone_ID = 'mesange'
exp_ID = '14-waves_007'

main_path = df.find_path(disk = 'Elements',year = '2024')
path2data = f'{main_path}{date}/Drones/{drone_ID}/matData/{exp_ID}/'

filelist = glob.glob(f'{path2data}*scaled.mat')
print(f'Files available : {filelist}')
file2load = filelist[0]
with h5py.File(file2load,'r') as fmat:
    print('Top-level keys : ', list(fmat.keys()))
    data = mat2py.mat_to_dict(fmat['m'],fmat['m'])

data = mat2py.transpose_PIVmat_fields(data)

    
#%% Reduce drone noise and compute FFT spectrum 
Vs = FT.supress_quadratic_noise(np.transpose(data['Vx'],(1,0,2)),data['x'],data['y'])
Vs = np.transpose(Vs,(1,0,2))

#%% Load image associated to this movie
img_list = glob.glob(f'{main_path}{date}/Drones/{drone_ID}/{exp_ID}/*exemple.tiff')
file2load = img_list[0]

img = cv.imread(file2load)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

fig, ax = plt.subplots()
ax.imshow(img)

xmin = 0
xmax = img.shape[1] - 1

crop = img[:,xmin:xmax,:]
fig, ax = plt.subplots()
ax.imshow(crop)

frasil_boxes = np.where(np.logical_and(data['PIXEL']['x_pix'] < xmax,data['PIXEL']['x_pix'] > xmin))[0]
idx_start = frasil_boxes[0]
idx_end = frasil_boxes[-1]  
  
#%% Define fig_folder and save folder
Dt = int(data['PIV_param']['Dt'])
fig_folder = f'{path2data}Plots/attenuation_from_disp_relation_Dt{Dt}/split_study_{xmin}to{xmax}/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

save_folder =f'{path2data}Results/pix_{xmin}to{xmax}/'
if not os.path.isdir(save_folder):
    os.mkdir(save_folder)

#%% Compute histogram
figname = f'{fig_folder}Histogram_{date}_{drone_ID}_{exp_ID}'
FT.histogram_PIV(Vs[idx_start:idx_end,:,:],data['PIV_param']['w'],figname)

#%% Compute FFT spectrum 
TF_spectrum,freq,FFT_t = FT.temporal_FFT(Vs,data['SCALE']['facq_t'],padding_bool = 1,add_pow2 = 0,output_FFT = True)

# save data in a structure
FFT_spectrum = {}
FFT_spectrum['TF_spectrum'] = TF_spectrum
FFT_spectrum['f'] = freq
#%%
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle _{x,y}(f) \; \mathrm{(u.a.)}$',labelpad = 5)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_ylim([1e-5,1e0])

figname = f'{fig_folder}TF_spectrum_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute space time spectrum 

N = Vs.shape[2]
Efk = FT.space_time_spectrum(Vs[idx_start:idx_end,:,:],1/data['SCALE']['fx'],data['SCALE']['facq_t'],add_pow2 = [1,1,0])

#%% Plot Afk 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
Amin = 5e-5 # change to adjust colormap
Amax = 1e-2 # change to adjust colormap
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
kbounds = [0.05,3.3] # bounds for k axis on Efk plot
fbounds = [0.05,1.3] # bounds for f axis on Efk plot
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}Efk_raw_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Look at demodulated fields

aspect_ratio = data

frequency_range = [0.1,0.8]
indices_freq = np.where(np.logical_and(Efk['f'] > frequency_range[0],Efk['f'] < frequency_range[1]))[0]

# create array of frequencies and amplitude 
freq = Efk['f'][indices_freq]
selected_FFT = FFT_t[:,:,indices_freq]

ani = animation_demodulated_field(selected_FFT, freq, parula_map, 4e-1)

# Save the created animation
file2save = f'animation_real_field_fmin_{frequency_range[0]}_fmax_{frequency_range[1]}'
file2save = file2save.replace('.','p')
file2save = f'{fig_folder}{file2save}.mp4'
ani.save(file2save)
print('Animation saved')


#%% Method for peaks detection
# step 1 - find peaks using gaussian convolution
# step 2 - filter and keep only relevant peaks
# step 3 - for each relevant peak, compute points on which gaussian will be fitted
# step 4 - save attenuation coefficient for each detected peak 

#%% Select a profile for a given wavevector and perform peak detection with gaussian fit over frequency 
Afk = Efk.copy()
# Build Gaussian to be convoluted (can change gaussian width to detect more or less peaks)
gaussian_width = 6
N = 8
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.6} # parameters for find_peaks
wavevector_range = [0.01,2.5] # range of wavevector spanned
frequency_range = [0,1.5] # range of frequency over which we look for peaks

S = matrix_peak_correlation_detection(Afk['E'], Afk['k'], Afk['f'], wavevector_range, frequency_range, gaussian, detec_param)
# change key name
S['k'] = S['x']
S['f'] = S['y']
del S['x']
del S['y']

#%% Filter detected peaks superposed with space-time spectrum
# Filter detected peaks
set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(S['k'], S['f'],S)

fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_respectto_f_fixed_k.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(filtered_properties,pf)

print('DONE.')
#%% Select a peak and cascade down to compute width over which we perform fitting
Afk = Efk.copy()

# load filtered properties 
file2load = f'{fig_folder}Filtered_peaks_respectto_f_fixed_k.pkl'
with open(file2load,'rb') as pf:
   filtered_properties = pickle.load(pf)
    
rel_height = detec_param['rel_height']

# build a dictionnary to save attenuation data
m = {'A':[],'f':[],'k':[],'lambda':[],'err_f':[],'err_lambda':[],'d':[]}

mask_f = np.where(np.logical_and(Efk['f'] < frequency_range[1],Efk['f'] > frequency_range[0]))[0]
Afk['f'] = Efk['f'][mask_f]
Afk['E'] = Efk['E'][mask_f,:]

indices = np.where(np.logical_and(Efk['k'] > wavevector_range[0],Efk['k'] < wavevector_range[1]))[0]
fig, ax = plt.subplots()
# fig,ax1 = plt.subplots()

sub_folder = f'{fig_folder}Lorentzian_fit/'
if not os.path.isdir(sub_folder):
    os.mkdir(sub_folder)

for idx in indices:

    ax.clear()
    
    detected_k = Afk['k'][idx]
    print(f'k  = {detected_k:.2f}')
    label_k = f'{detected_k:.2f}'
    
    # select peaks detected at this frequency
    mask = np.where(filtered_properties['idx'] == idx)[0]
    peaks = filtered_properties['peaks'][mask]
    
    section = Afk['E'][:,idx]
    norm = section/section.max()
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    
    ax.plot(Afk['f'],norm,'-o',label = label_k)
    ax.plot(Afk['f'],norm_correlated)
    
    for p in peaks:

        # ax1.clear()
        
        points2fit = indices2fit(norm_correlated,Afk['f'],p,rel_height)
        # plot points
        ax.plot(Afk['f'][p],norm_correlated[p],'s')
        ax.plot(Afk['f'][points2fit],norm[points2fit],'d')  
        if section[p] > section[points2fit].min():
            
            x = Afk['f'][points2fit]*2*np.pi # pulsation
            y_exp = (section[points2fit] - section[points2fit].min())/(section[points2fit].max() - section[points2fit].min())
            
            # gaussian fit 
            popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                                 bounds = ([x[0],1e-4],[x[-1],1]))
            err_coeff = np.sqrt(np.diag(pcov))
            print(f'f0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')
            
            xfit = np.linspace(x[0],x[-1],len(x)) #pulsation
            yth = lorentzian(xfit, popt[0], popt[1])
            
            # ax1.plot(x,y_exp,'o')
            # ax1.plot(xfit,yth)
            
            yplot = yth*(section[points2fit].max() - section[points2fit].min()) + section[points2fit].min()
            yplot = yplot/section.max()
            ax.plot(xfit/2/np.pi,yplot,'r-')  
            
            # compute distance to fit
            dist2fit = np.sum((y_exp - yth)**2)/np.sum(y_exp**2)
            m['d'].append(dist2fit)
            m['A'].append(Afk['E'][p,idx])
            m['f'].append(popt[0]/2/np.pi)
            m['k'].append(detected_k)
            m['lambda'].append(popt[1])
            m['err_f'].append(err_coeff[0]/2/np.pi)
            m['err_lambda'].append(err_coeff[1])
        
    # fancy plot and save it 
    ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 5)
    ax.set_ylabel(r'$A_{norm}$',labelpad = 5)
    ax.legend()
    
    k_txt = f'{detected_k:.2f}'.replace('.','p')
    figname = f'{sub_folder}Lorentzian_fit_respectto_f_fixed_k_{k_txt}'
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
        
    plt.pause(0.5)
        
for key in m.keys():
    m[key] = np.array(m[key])    
    
#%% Compute water depth hw

popt,pcov = scipy.optimize.curve_fit(lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi,m['k'],m['f'],
                                     sigma = m['err_f'],absolute_sigma = True)

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.errorbar(m['k'],m['f'],yerr = m['err_f'],fmt = '.',color = 'w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)
ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

hw = popt[0]
err_hw = np.sqrt(np.diag(pcov)[0])
title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'

k_fit = np.linspace(kbounds[0],kbounds[1],100)
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi

ax.plot(k_fit,y_exp,'r',label = title)
# y_bathy = bound_harmonicN(k_fit, 1, 3.63)/2/np.pi
# ax.plot(k_fit,y_bathy,'k--',label = 'bathy')
ax.legend()

print(f'H = {hw:.2f} ± {err_hw:.2f}')

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}Efk_hw_{hw_txt}_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# save water depth measurement
m['hw'] = hw
m['err_hw'] = err_hw

#%% Compute spatial attenuation coefficient 

time_att = np.column_stack((m['lambda'],m['err_lambda'])) # time attenuation coefficient and standard deviation
water_depth = [hw , err_hw] # water depth and standard deviation
alpha,err_alpha = time2space_attenuation(time_att, m['k'], water_depth) # compute spatial attenuation

xbounds = [0.1,1]
ybounds = [1e-3,1e0]
x_fit = np.linspace(xbounds[0],xbounds[1],100)

# fit by a powerlaw
coeffs,err_coeffs = powerlaw_fit(m['f'], alpha, err_alpha)
yth = coeffs[1]*x_fit**coeffs[0]

label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'

fig, ax = plt.subplots()
ax.errorbar(m['f'],alpha,yerr = err_alpha,xerr = m['err_f'],fmt = 'o')
ax.plot(x_fit,yth,linewidth = 2,label = label_th)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_from_temporal_attenuation_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  

# Save data 
m['alpha'] = alpha
m['err_alpha'] = err_alpha
m['power_law'] = {}
m['power_law']['B'] = coeffs[1]
m['power_law']['err_B'] = err_coeffs[1]
m['power_law']['beta'] = coeffs[0]
m['power_law']['err_beta'] = err_coeffs[0]

file2save = f'{save_folder}attenuation_from_Efk_fixed_k_respect2f_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(m,pf)
    
    
#%% Create main_results structure and save it 

main_results = {}
main_results['date'] = date
main_results['drone_ID'] = drone_ID
main_results['exp_ID'] = exp_ID
main_results['DRONE'] = data['DRONE']
main_results['SCALE'] = data['SCALE']
main_results['GPS'] = data['GPS']
main_results['t0_UTC'] = data['t0_UTC']

main_results['FFT_spectrum'] = FFT_spectrum
main_results['attenuation'] = m

file2save = f'{save_folder}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf :
    pickle.dump(main_results,pf)

print('Main results file saved !')







###############################################################################################################
################### Compute attenuation coefficient using fixed frequency #########################################################################
###############################################################################################################

#%% Detect harmonics for a given frequency and perform peak detection
# Build Gaussian to be convoluted
gaussian_width = 2
N = 6
gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
detec_param = {'prominence':1e-2,'rel_height':0.4} # parameters for find peaks function
frequency_range = [0.1,1.1] 
wavevector_range = [0.05,3.0]


S = matrix_peak_correlation_detection(Efk['E'].T, Efk['f'], Efk['k'], frequency_range, wavevector_range, gaussian, detec_param)
# change key name
S['f'] = S['x']
S['k'] = S['y']
del S['x']
del S['y']


#%% Show detected peaks superposed with space-time spectrum

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(S['k'], S['f'],S)

# save filtered_properties
file2save = f'{fig_folder}Filtered_peaks_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(filtered_properties,pf)
    
#%%
fig, ax = plt.subplots()
ax.plot(filtered_properties['k'],filtered_properties['f'],'.w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

figname = f'{fig_folder}Detected_fk_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Select a profile for a given frequency and perform peak detection with gaussian fit
Afk = Efk.copy()

# load filtered properties 
file2load = f'{fig_folder}Filtered_peaks_fixed_f_respect2k_{date}_{drone_ID}_{exp_ID}.pkl'
with open(file2load,'rb') as pf:
   filtered_properties = pickle.load(pf)
    
rel_height = detec_param['rel_height']
# build a dictionnary for peaks detected at fixed frequency
dict_f = {'A':[],'f':[],'k':[],'alpha':[],'err_k':[],'err_alpha':[],'d':[]}

mask_k = np.where(np.logical_and(Afk['k'] < wavevector_range[1],Afk['k'] > wavevector_range[0]))[0]
Afk['k'] = Efk['k'][mask_k]
Afk['E'] = Efk['E'][:,mask_k]

indices = np.where(np.logical_and(Afk['f'] > frequency_range[0],Afk['f'] < frequency_range[1]))[0]
fig, ax = plt.subplots()
# fig,ax1 = plt.subplots()

sub_folder = f'{fig_folder}Lorentzian_fit_respectto_k/'
if not os.path.isdir(sub_folder):
    os.mkdir(sub_folder)

for idx in indices:

    ax.clear()
    
    detected_f = Afk['f'][idx]
    print(f'k  = {detected_f:.2f}')
    label_k = f'{detected_f:.2f}'
    
    # select peaks detected at this frequency
    mask = np.where(filtered_properties['idx'] == idx)[0]
    peaks = filtered_properties['peaks'][mask]
    
    section = Afk['E'][idx,:]
    norm = section/section.max()
    correlated = scipy.signal.fftconvolve(norm,gaussian,mode = 'same')
    norm_correlated = correlated/correlated.max()
    
    ax.plot(Afk['k'],norm,'-o',label = label_k)
    ax.plot(Afk['k'],norm_correlated)
    
    for p in peaks:

        # ax1.clear()
        
        points2fit = indices2fit(norm_correlated,Afk['k'],p,rel_height)
        # plot points
        ax.plot(Afk['k'][p],norm_correlated[p],'s')
        ax.plot(Afk['k'][points2fit],norm[points2fit],'d')  
        if points2fit.size > 0 :
            if section[p] > section[points2fit].min():
                
                x = Afk['k'][points2fit]
                y_exp = (section[points2fit] - section[points2fit].min())/(section[points2fit].max() - section[points2fit].min())
                
                try : 
                    # gaussian fit 
                    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                                         bounds = ([x[0],1e-4],[x[-1],1]))
                    err_coeff = np.sqrt(np.diag(pcov))
                    print(f'f0 = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')
                    
                    xfit = np.linspace(x[0],x[-1],len(x))
                    yth = lorentzian(xfit, popt[0], popt[1])
                    
                    # ax1.plot(x,y_exp,'o')
                    # ax1.plot(xfit,yth)
                    
                    yplot = yth*(section[points2fit].max() - section[points2fit].min()) + section[points2fit].min()
                    yplot = yplot/section.max()
                    ax.plot(xfit,yplot,'r-')  
                    
                    # compute distance to fit
                    dist2fit = np.sum((y_exp - yth)**2)/np.sum(y_exp**2)
                    dict_f['d'].append(dist2fit)
                    dict_f['A'].append(Afk['E'][idx,p])
                    dict_f['k'].append(popt[0])
                    dict_f['f'].append(detected_f)
                    dict_f['alpha'].append(popt[1])
                    dict_f['err_k'].append(err_coeff[0])
                    dict_f['err_alpha'].append(err_coeff[1])
                    
                except RuntimeError :
                    print(f'Fit did not work for frequency f = {detected_f} Hz')
        
    # fancy plot and save it 
    ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$',labelpad = 5)
    ax.set_ylabel(r'$A_{norm}$',labelpad = 5)
    ax.legend()
    
    f_txt = f'{detected_f:.3f}'.replace('.','p')
    figname = f'{sub_folder}Lorentzian_fit_respectto_k_fixed_f_{f_txt}'
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
        
    plt.pause(0.5)
        
for key in dict_f.keys():
    dict_f[key] = np.array(dict_f[key])    
    
#%% 
fig, ax = plt.subplots()
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

ax.errorbar(dict_f['k'],dict_f['f'],xerr = dict_f['err_k'],fmt = 'o')
ax.plot(filtered_properties['k'],filtered_properties['f'],'o')    
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)
#%% DISCRIMINATE DIFFERENT HARMONICS
N = np.arange(1,5)
hw = 4
nb_points = len(dict_f['k'])
k_list = np.linspace(dict_f['k'].min(),dict_f['k'].max(),nb_points)
harmonics = np.zeros((k_list.shape[0],N.shape[0]))
for i in range(len(N)):
    harmonics[:,i] = bound_harmonicN(k_list, N[i], hw)

fig, ax = plt.subplots()
ax.plot(dict_f['k'],dict_f['f'],'o')
for i in range(len(N)):
    ax.plot(k_list,harmonics[:,i]/2/np.pi)

# build a distance to each harmonic
D = np.zeros((nb_points,len(N)))
for i in range(len(N)):
    D[:,i] = np.abs(dict_f['f'] - bound_harmonicN(dict_f['k'], N[i], hw)/2/np.pi)

closest_harmonic = np.argmin(D,axis = 1) + 1

# correct closestharmonic
closest_harmonic[np.where(dict_f['k']<0.5)] = 1

mask = np.where(np.logical_and(closest_harmonic == 3,dict_f['k'] < 0.7))
closest_harmonic[mask] = 2

mask = np.where(np.logical_and(closest_harmonic == 4,dict_f['k'] < 1.03))
closest_harmonic[mask] = 3

dict_f['closest_harmonic'] = closest_harmonic

fig,ax = plt.subplots()
for n in N:
    mask = dict_f['closest_harmonic'] == n
    ax.plot(dict_f['k'][mask],dict_f['f'][mask],'o')

#%% create structure specific to harmonic 1
dict_harm_1 = {}
mask = dict_f['closest_harmonic'] == 1
for key in dict_f.keys():
    dict_harm_1[key] = dict_f[key][mask]

# get rid off wrong points 
filtered_x,filtered_y,dict_harm_1 = scatter_filter.interactive_scatter_filter(dict_harm_1['k'], dict_harm_1['f'],dict_harm_1)


#%% Find water depth using detection at fixed frequency
# fit hw over harmonic 1
ktofit = dict_harm_1['k']
ftofit = dict_harm_1['f']

popt,pcov = scipy.optimize.curve_fit(lambda k,hw : bound_harmonicN(k, 1, hw)/2/np.pi,ktofit,ftofit)
k_fit = np.linspace(kbounds[0],kbounds[1],100)

set_graphs.set_matplotlib_param('single')

fig, ax = plt.subplots(figsize = (12,9))
ax.plot(ktofit,ftofit,'.w')
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower',interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(kbounds)
ax.set_ylim(fbounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

hw = popt[0]
err_hw = np.sqrt(np.diag(pcov)[0])
title = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
y_exp = bound_harmonicN(k_fit, 1, hw)/2/np.pi
ax.plot(k_fit,y_exp,'r',label = title)
print(f'H = {hw:.2f} ± {err_hw:.2f}')
ax.legend()

hw_txt = f'{hw:.2f}'.replace('.','p')
figname = f'{fig_folder}Efk_hw_{hw_txt}_fixed_f_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Keep only relevant points

set_graphs.set_matplotlib_param('single')
filtered_x,filtered_y,filtered_properties = scatter_filter.interactive_scatter_filter(dict_harm_1['f'], dict_harm_1['alpha'],
                                                                                      dict_harm_1)

for key in filtered_properties.keys():
    dict_harm_1[key] = filtered_properties[key]

#%% Compute spatial attenuation for harmonic 1

x = dict_harm_1['f']
y = dict_harm_1['alpha']
err_y = dict_harm_1['err_alpha']

# bounds for plot
xbounds = [0.1,1]
ybounds = [1e-3,1e0]
x_fit = np.linspace(xbounds[0],xbounds[1],100)

# power law fit taking into account uncertainties over y 
coeffs,err_coeffs = powerlaw_fit(x, y, err_y)
yth = coeffs[1]*x_fit**coeffs[0]
label_th = r'$y = ' + f'{coeffs[1]:.2f}'+ 'f^{' + f'{coeffs[0]:.3f}' + '}$'

fig, ax = plt.subplots()
# plot points detected at fixed f
ax.errorbar(x,y,yerr = err_y,fmt = 'o')
ax.plot(x_fit,yth,linewidth = 2,label = label_th)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xbounds)
ax.set_ylim(ybounds)
ax.legend()

figname = f'{fig_folder}Spatial_attenuation_from_fixed_f_{date}_{drone_ID}_{exp_ID}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')  

#%%
# store hw in dict_f
dict_harm_1['hw'] = hw
dict_harm_1['err_hw'] = err_hw

# store power law parameters 
dict_harm_1['power_law'] = {}
dict_harm_1['power_law']['B'] = coeffs[1]
dict_harm_1['power_law']['err_B'] = err_coeffs[1]
dict_harm_1['power_law']['beta'] = coeffs[0]
dict_harm_1['power_law']['err_beta'] = err_coeffs[0]




#%% Save data in main_results 

# load data from main_results 
file_main_results = f'{save_folder}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
if os.path.isfile(file_main_results):
    print('Main results already exists')
    with open(file_main_results,'rb') as pf :
        main_results = pickle.load(pf)
else :
    print('Main results file does not exist !')
    
main_results['attenuation_fixed_f'] = {}
main_results['attenuation_fixed_f']['all_harms'] = dict_f
main_results['attenuation_fixed_f']['harm_1'] = dict_harm_1

print(main_results.keys())

# save main_results 
with open(file_main_results,'wb') as pf :
    main_results = pickle.dump(main_results,pf)
print('Data detected at fixed f also saved in main_results !')


###############################################################
#%% Compare attenuation and dispersion relation 
###############################################################

file2load = f'{save_folder}main_results_{date}_{drone_ID}_{exp_ID}.pkl'
if os.path.isfile(file_main_results):
    print('Main results already exists')
    with open(file_main_results,'rb') as pf :
        main_results = pickle.load(pf)
    print('Main results loaded')
else :
    print('Main results file does not exist !')


#%% Compare dispersion relation

set_graphs.set_matplotlib_param('single')
x_bounds = [0.05,3.3]
y_bounds = [0.05,1.3]
x_fit = np.linspace(x_bounds[0],x_bounds[1],100)

fig, ax = plt.subplots(figsize = (12,9))
# plot space-time spectrum
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', vmin = Amin,vmax = Amax,
              origin = 'lower',interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))
cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# plot points detected at fixed k
ax.scatter(main_results['attenuation']['k'],main_results['attenuation']['f'],color = 'w',s = 15)
hw = main_results['attenuation']['hw']
err_hw = main_results['attenuation']['err_hw']
yth = bound_harmonicN(x_fit, 1, hw)/2/np.pi
label_th = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
ax.plot(x_fit,yth,'r--',label = label_th)

# plot points detected at fixed f 
ax.scatter(main_results['attenuation_fixed_f']['all_harms']['k'],main_results['attenuation_fixed_f']['all_harms']['f'], s = 15,
           color = 'tab:orange',marker = 'd')
hw = main_results['attenuation_fixed_f']['harm_1']['hw']
err_hw = main_results['attenuation_fixed_f']['harm_1']['err_hw']
yth = bound_harmonicN(x_fit, 1, hw)/2/np.pi
label_th = r'$H = ' + f'{hw:.2f}' + '\pm' + f' {err_hw:.2f}' +'$'
ax.plot(x_fit,yth,'--',color = 'k',label = label_th)

ax.set_xlim(x_bounds)
ax.set_ylim(y_bounds)

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
ax.legend()

figname = f'{fig_folder}Comparison_dispersionrelation_time_detection_VS_space_detection'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')



#%% Compare attenuation coefficients

set_graphs.set_matplotlib_param('single')
x_bounds = [1e-1,1e0]
y_bounds = [1e-3,1e0]
x_fit = np.linspace(x_bounds[0],x_bounds[1],100)

fig, ax = plt.subplots(figsize = (12,9))
# plot coefficients from space detection (fixed k)
ax.errorbar(main_results['attenuation']['f'],main_results['attenuation']['alpha'],
            yerr = main_results['attenuation']['err_alpha'],fmt = 'o',color = 'tab:blue')
B = main_results['attenuation']['power_law']['B']
beta = main_results['attenuation']['power_law']['beta']
yth = B*x_fit**beta
label_th = r'$\alpha(f) = ' + f'{B:.2f}'+ 'f^{' + f'{beta:.3f}' + '}$'
ax.plot(x_fit,yth,'-',color = 'tab:blue',label = label_th)

# plot coefficients from time detection (fixed f)
ax.errorbar(main_results['attenuation_fixed_f']['harm_1']['f'],main_results['attenuation_fixed_f']['harm_1']['alpha'],
            yerr = main_results['attenuation_fixed_f']['harm_1']['err_alpha'],fmt = 'o',color = 'tab:orange')
B = main_results['attenuation_fixed_f']['harm_1']['power_law']['B']
beta = main_results['attenuation_fixed_f']['harm_1']['power_law']['beta']
yth = B*x_fit**beta
label_th = r'$\alpha(f) = ' + f'{B:.2f}'+ 'f^{' + f'{beta:.3f}' + '}$'
ax.plot(x_fit,yth,'-',color = 'tab:orange',label = label_th)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(x_bounds)
ax.set_ylim(y_bounds)
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\alpha \; \mathrm{(m^{-1})}$')

ax.legend()

figname = f'{fig_folder}Comparison_spatial_attenuation_space_detection_VS_time_detection'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
















