# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:58:17 2025

@author: sebas
"""

import numpy as np
import math
import matplotlib.pyplot as plt 
import matplotlib as mpl
import matplotlib.colors as colors 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py 
import glob
import os 
import pickle
import scipy.signal
from scipy.signal import find_peaks

import cv2 as cv

# import modules 
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs

parula_map = matcmaps.parula()
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/Referee_answer/'

#%%
font_size_medium = 30
font_size_small = round(0.75*font_size_medium)
plt.rc('font', size=font_size_medium)          # controls default text sizes
plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

fig_size = (8,6)
img_quality = 100 # dpi to save images 

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% FUNCTION SECTION 

def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y

#------------------------------------------------------------------------------------------------------------------------

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

def freq_hydrostatic(h,rho):
    """ Compute resonance frequency according to hydrostatic equilibrium (without added mass)
    Inputs: - h, height of cylindrical floater
            - rho, density of floater
    Output: - f0, resonance frequency """
    
    g = 9.81
    rho_w = 1e3
    f0 = np.sqrt(rho_w*g/(rho*h))/2/np.pi
    return f0


#%%  Load data 

base = 'F:/Waves_reconstruction_wilson/8_hammering_2024_01_02/'
exp_ID = 'D6cm_h15mm'
test_ID = 'test2'

path2data = f'{base}{exp_ID}/center/{test_ID}/'
filelist = glob.glob(f'{path2data}*Height_in_cm_surf.mat')
file2load = filelist[0]

data = {}
with h5py.File(file2load, 'r') as fmat:
    
    list_keys = list(fmat.keys())
    print('Top-level keys : ', list_keys)
    for i in range(1,len(list_keys) - 1):
        key = list_keys[i]
        data[key] = mat2py.mat_to_dict(fmat[key],fmat[key])
        
    data['table'] = mat2py.mat_to_dict(fmat['table'],fmat['#refs#'])
    
    data['H_ww'] = np.flip(data['H_ww'],axis = 2)

#%% Show a frame 

frame = 250
field = data['H_ww'][frame,:,:]
extents = np.array([data['X'][0] , data['X'][-1], data['Y'][0], data['Y'][-1]])/data['fx']*1e2

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
imsh = ax.imshow(field.T,cmap = parula_map, aspect = 'equal', 
                 origin = 'lower', extent = extents, vmin = -0.05, vmax = 0.05) 

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$')


ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)


#%% Compute FFT in time averaged over space 

fps = 200 # frame rate 
idx_start = 140
idx_end = 350
TF_spectrum,freq,FFT_t = FT.temporal_FFT(np.transpose(data['H_ww'],(1,2,0))[:,:,idx_start:idx_end],fps,padding_bool = 1,
                                         add_pow2 = 1,output_FFT = True)

print('Temporal FFT computed !')
#%% Plot spectrum log-log

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([3e-5,1e-2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

#%% Plot spectrum lin-lin

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'-')

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_ylim([3e-5,1e-2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

#%% Find peak and gets its width at half maximum

rel_height = 0.5
delta_f = freq[1] - freq[0]
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
prominences = scipy.signal.peak_prominences(TF_spectrum,peaks)
widths = scipy.signal.peak_widths(TF_spectrum,peaks,rel_height = rel_height)


fig, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')
ax.plot(freq[peaks],TF_spectrum[peaks],'rd')

ax.vlines(x=freq[peaks], ymin=TF_spectrum[peaks] - prominences[0],
        ymax = TF_spectrum[peaks],color = 'g')
ax.hlines(y=widths[1], xmin=widths[2]*delta_f,
        xmax=widths[3]*delta_f,color = 'g')

# ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_ylim([3e-5,1e-2])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

alpha = widths[0]*2*np.pi*delta_f
print(alpha)

#%% Find peak and fit it by a lorentzian 

# find peaks 
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
first_peak = peaks[0] # get first peak 

# normalize spectrum
norm_spectrum =(TF_spectrum - TF_spectrum.min())/(TF_spectrum.max() - TF_spectrum.min())

# find points to fit 
points2fit = indices2fit(norm_spectrum,freq,first_peak,0.85)

# fit by a lorentzian
x = freq[points2fit]*2*np.pi # pulsation
y_exp = norm_spectrum[points2fit]
popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                     bounds = ([x[0],1e-4],[x[-1],10]))
err_coeff = np.sqrt(np.diag(pcov))
print(f'Omega = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')

xfit = np.linspace(x[0],x[-1]) 
yth = lorentzian(xfit, popt[0], popt[1])
scaled_theory = yth*(TF_spectrum.max() - TF_spectrum.min()) + TF_spectrum.min()
label_th = r'$\frac{1}{\sqrt{1 + (\omega - \omega_0)^2/\lambda^2}}$'

set_graphs.set_matplotlib_param('single')
fix, ax = plt.subplots()
ax.plot(freq,TF_spectrum,'o')
ax.plot(freq[points2fit],TF_spectrum[points2fit],'ro')
ax.plot(xfit/2/np.pi,scaled_theory,'r',label = label_th)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([3e-5,1e-2])
ax.set_xlim([3e-1,1.5e2])
# ax.set_ylim([-1e-3,1e-2])
# ax.set_xlim([-1,40])

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(cm)}$')

ax.legend(loc = 'lower left')

figname = f'TF_spectrum_impulse_lorentzian_loglog_{exp_ID}_{test_ID}_Omega_{popt[0]:.2f}_lambda_{popt[1]:.2f}'
figname = figname.replace('.','p')
figname = f'{fig_folder}{figname}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#########################################################################################
#%% -------------------- Referee answer - New Figure 2c ---------------------------------
#########################################################################################

fig, ax = plt.subplots(figsize = (8,7))
ax.plot(freq,TF_spectrum*1e1,'-')
ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(mm)}$')

# create insert
axins = inset_axes(ax, width="50%", height="50%")
axins.plot(freq,TF_spectrum*1e1,'.')
axins.set_xscale('log')
axins.set_yscale('log')
axins.set_ylim([3e-4,1e-1])

# add fit in insert
# find peaks 
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
first_peak = peaks[0] # get first peak 

# normalize spectrum
norm_spectrum =(TF_spectrum - TF_spectrum.min())/(TF_spectrum.max() - TF_spectrum.min())

# find points to fit 
points2fit = indices2fit(norm_spectrum,freq,first_peak,0.85)

# fit by a lorentzian
x = freq[points2fit]*2*np.pi # pulsation
y_exp = norm_spectrum[points2fit]
popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                     bounds = ([x[0],1e-4],[x[-1],10]))
err_coeff = np.sqrt(np.diag(pcov))
print(f'Omega = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')

xfit = np.linspace(x[0],x[-1]) 
yth = lorentzian(xfit, popt[0], popt[1])
scaled_theory = yth*(TF_spectrum.max() - TF_spectrum.min()) + TF_spectrum.min()
label_th = r'$\frac{1}{\sqrt{1 + (\omega - \Omega)^2/\lambda^2}}$'

axins.plot(freq[points2fit],TF_spectrum[points2fit]*1e1,'r.')
axins.plot(xfit/2/np.pi,scaled_theory*1e1,'r',label = label_th)

figname = f'TF_spectrum_insert_lorentzian_linear_{exp_ID}_{test_ID}'
figname = f'{fig_folder}{figname}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')














##############################################################################
#%% ----------------------- New Figure 2 -----------------------
##############################################################################

# load matfile
path2folder_tiff = 'Y:/Banquise/Sebastien/Article_Wave_floats/8_hammering_2024_01_02/D6cm_h15mm/center/'
path2matfile = path2folder_tiff + '*.mat'

matfile = glob.glob(path2matfile)
with h5py.File(matfile[0], 'r') as fmat:
    H_ww = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    H_ww = mat2py.mat_to_dict(fmat['H_ww'],fmat['H_ww'])


#%% Show elevation field 

# change dimensions order 
H = np.transpose(H_ww, [2,1,0])
H = np.flip(H,axis = 0)

#%%
fx = 2.4167*1e3 # scale in pixels/meter
x = np.arange(np.size(H,1))/fx*100 # x-coordinates in cm
y = np.arange(np.size(H,0))/fx*100

idx_frame = 279
figure, ax = plt.subplots(figsize = fig_size)
wavefield = ax.imshow(H[:,:,idx_frame],cmap = parula_map ,vmin = -0.05, vmax = 0.05,
                      aspect = 'equal',origin = 'lower', extent = [x[0],x[-1],y[0],y[-1]])
ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)
# ax.set_title(f'$t = {t[idx_frame]:.2f} \;' + r' \mathrm{(s)}$')
# # Use make_axes_locatable to create a new axis for the colorbar
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)

# # Add the colorbar to the new axis
# cb = fig.colorbar(wavefield, cax=cax)

cbar = plt.colorbar(wavefield, ax = ax,shrink = 0.78)
cbar.set_label(r'$\zeta \; \mathrm{(cm)}$',labelpad = 5)


#%% Load image
filename = path2folder_tiff + 'test1/' + 'Basler_a2A1920-160ucBAS__40232065__20240102_145223416_0175.tiff'

img = cv.imread(filename)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
img = np.flip(img,axis = 0)

 
#%% Load data for centered impulse 
path2wavefield_data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Data/Resonance_freq/'
file2load = path2wavefield_data + 'Centered_impulse_resonance_freq.pkl'
with open(file2load,'rb') as pfile:
    S_center = pickle.load(pfile)


#%% Create new figure 2 

fig, ax = plt.subplots(2,2,figsize = (12,9),layout = 'constrained')

# TF spectrum
spectrum_window = (1,0)
current_ax = ax[spectrum_window]
current_ax.plot(freq,TF_spectrum*1e1,'-')
current_ax.set_xlabel(r'$f \; \mathrm{(Hz)}$')
current_ax.set_ylabel(r'$\langle \hat{\zeta} \rangle _{x,y} \; \mathrm{(mm)}$')

# create insert
axins = inset_axes(current_ax, width="50%", height="50%")
axins.plot(freq,TF_spectrum*1e1,'.')
axins.set_xscale('log')
axins.set_yscale('log')
axins.set_ylim([3e-4,1e-1])

# add fit in insert
# find peaks 
peaks,properties = find_peaks(TF_spectrum, prominence = 1e-3,width = 3)
first_peak = peaks[0] # get first peak 

# normalize spectrum
norm_spectrum =(TF_spectrum - TF_spectrum.min())/(TF_spectrum.max() - TF_spectrum.min())

# find points to fit 
points2fit = indices2fit(norm_spectrum,freq,first_peak,0.85)

# fit by a lorentzian
x = freq[points2fit]*2*np.pi # pulsation
y_exp = norm_spectrum[points2fit]
popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),x,y_exp,
                                     bounds = ([x[0],1e-4],[x[-1],10]))
err_coeff = np.sqrt(np.diag(pcov))
print(f'Omega = {popt[0]:.2f} ± {err_coeff[0]:.2f} and alpha = {popt[1]:.2e} ± {err_coeff[1]:.2e}')

xfit = np.linspace(x[0],x[-1]) 
yth = lorentzian(xfit, popt[0], popt[1])
scaled_theory = yth*(TF_spectrum.max() - TF_spectrum.min()) + TF_spectrum.min()
label_th = r'$\frac{1}{\sqrt{1 + (\omega - \Omega)^2/\lambda^2}}$'

axins.plot(freq[points2fit],TF_spectrum[points2fit]*1e1,'r.')
axins.plot(xfit/2/np.pi,scaled_theory*1e1,'r',label = label_th)

# Image 
img_window =(0,0)
current_ax = ax[img_window]
current_ax.imshow(img,aspect = 'equal',origin = 'lower', extent = [x[0],x[-1],y[0],y[-1]])
current_ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
current_ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)

# Wave field 
wavefield_window = (0,1)
current_ax = ax[wavefield_window]

wavefield = current_ax.imshow(H[:,:,idx_frame]*10,cmap = parula_map ,vmin = -0.5, vmax = 0.5,
                      aspect = 'equal',origin = 'lower', extent = [x[0],x[-1],y[0],y[-1]])
current_ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
current_ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)

cbar = plt.colorbar(wavefield, ax = current_ax,shrink = 0.74)
cbar.set_label(r'$\zeta \; \mathrm{(mm)}$',labelpad = 5)
# # cbar.set_ticks([-1,0,1])
# # cbar.set_ticklabels([f'{l:.1e}' for l in [-5e-2,0,5e-2]])

# divider = make_axes_locatable(current_ax)
# cax = divider.append_axes("right", size="2%", pad=0.1)
# cbar = plt.colorbar(wavefield,cax = cax)
# cbar.set_label(r'$\zeta \; \mathrm{(cm)}$',labelpad = 5)



# Resonance frequencies 
resonance_window = (1,1)
current_ax = ax[resonance_window]

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# theory
rho = 953 # density of black polypropylene
h_th = np.linspace(2e-3,30e-3,100) # array of height values
yth = freq_hydrostatic(h_th, rho)

scatt = scatter = current_ax.scatter(S_center['H'],S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm,edgecolors = 'k', zorder = 1)
current_ax.plot(h_th*1e3,yth,'k--',zorder= 0 )

current_ax.set_xlim(3,22)
current_ax.set_ylim(2.5,6)

cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

current_ax.set_xlabel(r'$h \; \mathrm{(mm)}$',labelpad = 5)
current_ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5)  

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])

figname = 'Figure_2_subplots_referee_corrections_with_lorentzian_fit'
figname = f'{fig_folder}{figname}'

plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')



















