# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 17:13:37 2024

@author: sebas
"""
import numpy as np
import cmath
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import cv2 as cv
import glob
import os
import pickle 
import h5py

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT

#%% Parameters for plots
# picture on whole page width : fig_size = (8,6), font_size_medium = 22, font_size_small = 0.75*font_size_medium

# picture on half page width : fig_size = (8,6), font_size_medium = 30, font_size_small  =  0.75*font_size_medium

# quadrant of 4 figures : fig_size = (6.4,4.8), font_size_medium = 24 (or 20), font_size_small = 0.75*font_size_medium 

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

fig_folder = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Figures_article/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
    
# Create a new blues colormap
# create new Blues colormap
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% Function section 

def alpha_ellipsoid(epsilon):
    """ Compute correction of added mass when considering an ellipsoid shape 
    - epsilon : aspect ratio of the ellipsoid : h/R for a cylinder"""
    
    alpha = np.zeros(np.shape(epsilon))
    for i,epp in enumerate(epsilon):
        if epp == 1:
            alpha[i] = 1
        elif epp > 1:
            acos_matlab = -cmath.acos(epp)
            a = (cmath.sqrt(1/epp**2 - 1) - acos_matlab)/(acos_matlab - 
                                                                   epp*cmath.sqrt(1 - epp**2))
            alpha[i] = a.real
        else :
            alpha[i] = (np.sqrt(1/epp**2 - 1) - np.arccos(epp))/(np.arccos(epp) - 
                                                                   epp*np.sqrt(1 - epp**2))
        
    return alpha    
    

def freq_Kyle_ellipsoid(R,h,rho):
    """ Compute theoretical resonance frequency using Kyle's model. We here compute the added mass of
    an oblate ellipsoidal object moving in a fluid : 
        - R : object radius in meter
        - h : object height in meter
        - rho : object volumic mass in kg/m3 """

    rho_w = 1e3
    g = 9.81
    r = rho/rho_w
    
    epsilon =r*h/R 
    
    alpha = alpha_ellipsoid(epsilon)

    omega = np.sqrt(g/(r*h + 2*R*epsilon*alpha/3))
    f_th = omega/2/np.pi

    return f_th  

def freq_Kyle_thin_disk(R,h,rho):
    """ Compute theoretical resonance frequency using Kyle's model. We here compute the added mass of 
    a thin disk (nearly 2D) moving in a fluid : 
        - R : object radius in meter
        - h : object height in meter
        - rho : object volumic mass in kg/m3
    """
    rho_w = 1e3
    g = 9.81
    r = rho/rho_w
    
    omega = np.sqrt(g/(r*h + 4*R/3/np.pi))
    f_th = omega/2/np.pi
    
    return f_th

def freq_hydrostatic(h,rho):
    """ Compute resonance frequency according to hydrostatic equilibrium (without added mass)
    Inputs: - h, height of cylindrical floater
            - rho, density of floater
    Output: - f0, resonance frequency """
    
    g = 9.81
    rho_w = 1e3
    f0 = np.sqrt(rho_w*g/(rho*h))/2/np.pi
    return f0
        
########################################################################################
#%% LOAD DATA OBTAINED FROM WAVE FIELDS - NON AVERAGED 
######################################################################################## 

path2wavefield_data = 'W:/Banquise/Sebastien/Test_Disks_Wilson_Seb/'
file2load = path2wavefield_data + 'Structure_non_averaged_resonance_freq_waves.mat'

with h5py.File(file2load, 'r') as fmat:
    data_wave = {}
       
    print('Top-level keys : ', list(fmat.keys()))

    data_wave = mat2py.mat_to_dict(fmat['S_navg'],fmat['S_navg'])
    
#%% Keep only data associated to centered impulse

S_center = {}
mask = data_wave['center_test'] == 1
relevant_keys = ['R','H','f']
for key in relevant_keys:
    S_center[key] = data_wave[key][mask]
    
#%% Save structure associated to centered impulse
file2save = path2wavefield_data + 'Centered_impulse_resonance_freq.pkl'

if os.path.isfile(file2save):
    print('File already exists !')
else : 
    with open(file2save,'wb') as pfile:
        pickle.dump(S_center,pfile)
    
#%% Load data for centered impulse 
path2wavefield_data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Data/Resonance_freq/'
file2load = path2wavefield_data + 'Centered_impulse_resonance_freq.pkl'
with open(file2load,'rb') as pfile:
    S_center = pickle.load(pfile)


#%% Figure 2 : f0 as a function of H, colormap R

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig, ax = plt.subplots(figsize = fig_size)
scatter = ax.scatter(S_center['H'],S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm,edgecolors = 'k', zorder = 2)
# ax.grid(zorder = 1)

ax.set_xlim(3,22)
ax.set_ylim(2.5,6)
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

ax.set_xlabel(r'$h \; \mathrm{(mm)}$',labelpad = 5)
ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5)  

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])

figname = fig_folder + 'f_center_wavefield_VS_H_cmap_R'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
    

#%% Figure 4 : Experimental data VS Kyle elliptic model

rho = 953
f_th = freq_Kyle_ellipsoid(S_center['R']*1e-3,S_center['H']*1e-3,rho)
xmin = 2.5
xmax = 6
xrange = np.linspace(2,6,100)

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig, ax = plt.subplots(figsize = fig_size)
scatter = ax.scatter(f_th,S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm, edgecolors = 'k', zorder = 3)
ax.plot(xrange,xrange,'k--',zorder = 2)
# ax.grid(zorder = 1)

ax.set_xlim(xmin,xmax)
ax.set_ylim(2.5,6)
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

ax.set_xlabel(r'$\frac{1}{2 \pi} \sqrt{\frac{g}{\frac{\rho}{\rho_w}h + \frac{2}{3}R \epsilon \alpha(\epsilon)}}$'
              ,labelpad = 5)
ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5) 

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
    
figname = fig_folder + 'f_center_wavefield_VS_Kyle_ellipsoid'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
    

#%% Figure 8 - Supplementary - Knocking resonance freq VS Kyle spheric model


rho = 953
f_th = freq_Kyle_thin_disk(S_center['R']*1e-3,S_center['H']*1e-3,rho)
xmin = 2.5
xmax = 6
xrange = np.linspace(xmin,xmax,100)

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

fig, ax = plt.subplots(figsize = fig_size)
scatter = ax.scatter(f_th,S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm, edgecolors = 'k', zorder = 3)
ax.plot(xrange,xrange,'k--',zorder = 2)
# ax.grid(zorder = 1)

ax.set_xlim(xmin,xmax)
ax.set_ylim(2.5,6)
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

ax.set_xlabel(r'$\frac{1}{2 \pi} \sqrt{\frac{g}{\frac{\rho}{\rho_w}h + \frac{4}{3 \pi}R }}$'
              ,labelpad = 5)
ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5) 

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
    
figname = fig_folder + 'f_center_wavefield_VS_Kyle_thin_disk'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Referee answer - f0 VS h - Kyle elliptical model 

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

rho = 953 # density of black polypropylene
radius_array = np.array([7.5,10,15,20,30]) # radius in millimeter
h_th = np.linspace(2e-3,30e-3,100) # array of height values

# compute theory curves
fth_Kyle = np.zeros((len(radius_array),len(h_th)))
for i,radius in enumerate(radius_array):
    current_fth = freq_Kyle_ellipsoid(radius*1e-3, h_th , rho)
    fth_Kyle[i,:] = current_fth

# plot experimental data 
fig, ax = plt.subplots(figsize = fig_size)
scatter = ax.scatter(S_center['H'],S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm,edgecolors = 'k', zorder = len(radius_array) + 1)
# ax.grid(zorder = 1)

# plot theory curves
for i in range(fth_Kyle.shape[0]):
    current_color = new_blues(norm(radius_array[i]))
    ax.plot(h_th*1e3,fth_Kyle[i,:],'--',color = current_color,
            zorder = i,lw = 2)

ax.set_xlim(3,22)
ax.set_ylim(2.5,6)
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

ax.set_xlabel(r'$h \; \mathrm{(mm)}$',labelpad = 5)
ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5)  

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])

# figname = fig_folder + 'f_center_wavefield_VS_H_cmap_R_with_Kyle_ellipsoid_theory'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Referee answer - f0 VS hydrostatic model







#####################################################################################################
################################### WAVE EMISSION ###################################################

#%% Figure 6 - Waves emission, transverse direction 

path2transdata = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Data/'
file2load = path2transdata + 'Structure_TransverseData.pkl'

with open(file2load,'rb') as pfile:
    data_trans = pickle.load(pfile)
    

#%% Save data_trans in a pickle file 
path2transdata = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Waves_float/Data/'
file2save = path2transdata + 'Structure_TransverseData.pkl'

if os.path.isfile(file2save):
    print('File already exists !')
else : 
    with open(file2save,'wb') as pfile:
        pickle.dump(data_trans,pfile)


#%%

mask = data_trans['H'] == 10

keys_R = np.unique(data_trans['R'][mask])[:-1]
Rmin = 5
Rmax = 20

# norm = (keys_R - Rmin)/(Rmax - Rmin)
bounds = np.array([7,8,12,18,22])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

fig, ax = plt.subplots(figsize = fig_size)

for i in range(len(keys_R)):
    current_R = keys_R[i]
    current_mask = np.logical_and(mask,data_trans['R'] == current_R)
    
    y = data_trans['Amax'][current_mask]/data_trans['Aw'][current_mask]
    x = data_trans['f_ex'][current_mask]/data_trans['f_center'][0][current_mask]
    
    current_color = new_blues(norm(current_R))
    ax.plot(x,y,'o-',c = current_color, lw = 2, mec = 'k', ms = 8)
        

# ax.grid(zorder = 1)

ax.set_xlim(0.55,1.7)
ax.set_ylim(0,0.8)
ax.set_xlabel(r'$f/f_0$',labelpad = 5)
ax.set_ylabel(r'$A_o/A_w$',labelpad = 5) 


# figname = fig_folder + 'ratio_Ao_Aw_VS_ratio_fex_f0'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')


#%% Compute frequency at which Ao/Aw is minimal 

keys_R = np.unique(data_trans['R'])
keys_H = np.unique(data_trans['H'])

S = {'R':[],'H':[],'f_min':[],'f0':[],'f_min_error':[],'f0_error':[]}


for i,r in enumerate(keys_R):
    for j,h in enumerate(keys_H):
        
        mask = np.logical_and(data_trans['R'] == r, data_trans['H'] == h)
        if np.any(mask):
            ratio_AoAw = data_trans['Amax'][mask]/data_trans['Aw'][mask]
            f_ex = data_trans['f_ex'][mask]
            
            idx_min = np.argmin(ratio_AoAw)
            # if idx_min < len(f_ex) - 1:
            #     idx_subpix,max_subpix = FT.subpix_precision(ratio_AoAw, idx_min)
            #     f_min = f_ex[idx_min] + (idx_subpix - idx_min)*(f_ex[idx_min] - f_ex[idx_min-1])
            # else :
            f_min = f_ex[idx_min]
            if idx_min < len(f_ex) - 1:
                error_f_min = 0.5*(abs(f_min - f_ex[idx_min - 1]) + abs(f_min - f_ex[idx_min + 1]))
                error_f_min = 0.5*error_f_min
            else : 
                error_f_min = 0.5*abs(f_min - f_ex[idx_min - 1])
            
            S['R'].append(r)
            S['H'].append(h)
            S['f_min'].append(f_min)
            S['f_min_error'].append(error_f_min)
            S['f0'].append(data_trans['f_center'][0,mask][0])
            S['f0_error'].append(data_trans['f_center'][1,mask][0]) 
    

for key in S.keys():
    S[key] = np.array(S[key])

#%% Plot frequency of minimal emission VS resonance frequency

xmin = 2.5
xmax = 6
xrange = np.linspace(xmin,xmax,100)

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)


fig, ax = plt.subplots(figsize = fig_size)
scatter = ax.scatter(S['f0'],S['f_min'], c = S['R'], s = 70,cmap = new_blues,
           norm = norm, edgecolors = 'k', zorder = 4)
err_bar = ax.errorbar(S['f0'],S['f_min'], fmt = 'none',yerr = S['f_min_error'], xerr = 2*S['f0_error'], 
                      capsize = 4,ecolor = 'k',zorder = 3)
ax.plot(xrange,xrange,'k--',zorder = 2)
# ax.grid(zorder = 1)

ax.set_xlim(xmin,xmax)
ax.set_ylim(2.5,6)

ax.set_yticks([3,4,5,6])
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

ax.set_xlabel(r'$f_0 \; \mathrm{(Hz)}$'
              ,labelpad = 5)
ax.set_ylabel(r'$f_{min}(A_o/A_w) \; \mathrm{(Hz)}$',labelpad = 5) 

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])

# figname = fig_folder + 'fmin_VS_f0'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')



#%% Supplementary Figure 10 

mask = data_trans['R'] == 15

keys_H = np.unique(data_trans['H'][mask])
Rmin = 5
Rmax = 20

# norm = (keys_R - Rmin)/(Rmax - Rmin)
bounds = np.array([2.5,7.5,12.5,17.5,22.5])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
full_oranges = mpl.colormaps['Oranges'].resampled(256)
new_oranges = colors.ListedColormap(full_oranges(np.linspace(0.2,1,256)))

fig, ax = plt.subplots(figsize = fig_size)

for i in range(len(keys_H)):
    current_H = keys_H[i]
    current_mask = np.logical_and(mask,data_trans['H'] == current_H)
    
    y = data_trans['Amax'][current_mask]/data_trans['Aw'][current_mask]
    x = data_trans['f_ex'][current_mask]/data_trans['f_center'][0][current_mask]
    
    current_color = new_oranges(norm(current_H))
    ax.plot(x,y,'o-',c = current_color, lw = 2, mec = 'k', ms = 8)
        
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=new_oranges), ax=ax)
midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])
cbar.set_label(r'$H \; \mathrm{(mm)}$',labelpad = 5)
ax.grid(zorder = 1)


ax.set_xlim(0.4,2.4)
ax.set_ylim(0,0.7)
ax.set_xlabel(r'$\frac{f}{f_0}$',labelpad = 5)
ax.set_ylabel(r'$\frac{A_o}{A_w}$',labelpad = 5) 

figname = fig_folder + 'ratio_Ao_Aw_VS_ratio_fex_f0_H_dependancy'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Figure 2 : Wave field, picture and FFT spectrum 

# path2folder_tiff = 'Y:/Banquise/Sebastien/Article_Wave_floats/8_hammering_2024_01_02/D6cm_h15mm/center/'
path2folder_tiff = 'F:/Waves_reconstruction_wilson/Article_wave_floats/'
path2matfile = path2folder_tiff + '*.mat'

matfile = glob.glob(path2matfile)

with h5py.File(matfile[0], 'r') as fmat:
    H_ww = {}
    
    print('Top-level keys : ', list(fmat.keys()))

    H_ww = mat2py.mat_to_dict(fmat['H_ww'],fmat['H_ww'])

#%% Show H_ww

# change dimensions order 
H = np.transpose(H_ww, [2,1,0])
H = np.flip(H,axis = 0)
# PARULA COLORMAP 
parula_map = matcmaps.parula()
#%%
fx = 2.4167*1e3 # scale in pixels/meter
x = np.arange(np.size(H,1))/fx*100 # x-coordinates in cm
y = np.arange(np.size(H,0))/fx*100

fps = 100 # frame per second
t = np.arange(np.size(H,2))/fps

idx_impact = 162
t_post_impact = t - t[idx_impact]
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

# figname = fig_folder + 'wavefield_R30mm_h15mm_center_test1_frame_279'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Plot associated frame 

# filename = path2folder_tiff + 'test1/' + 'Basler_a2A1920-160ucBAS__40232065__20240102_145223416_0175.tiff'
filename = path2folder_tiff + 'Basler_a2A1920-160ucBAS__40232065__20240102_145223416_0175.tiff'
img = cv.imread(filename)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
img = np.flip(img,axis = 0)
#%%

img = np.flip(img,axis = 0)
fig, ax = plt.subplots(figsize = fig_size)
ax.imshow(img,aspect = 'equal',origin = 'lower', extent = [x[0],x[-1],y[0],y[-1]])
ax.set_xlabel(r'$x \; \mathrm{(cm)}$',labelpad = 5)
ax.set_ylabel(r'$y \; \mathrm{(cm)}$',labelpad = 5)

figname = fig_folder + 'Image_R30mm_h15mm_center_test1_frame_175'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')

#%% Compute and plot FFT spectrum


cut_H = H[:,:,140:380]
fps = 200
TF_spectrum,freq = FT.temporal_FFT(cut_H,fps,1,2,output_FFT = False)


#%%
imax = np.argmax(TF_spectrum)
idx_subpix,max_subpix = FT.subpix_precision(TF_spectrum, imax)
f_res = freq[imax] + (idx_subpix - imax)*(freq[imax] - freq[imax-1])
fig, ax = plt.subplots(figsize = fig_size)
ax.loglog(freq,abs(TF_spectrum),'o',zorder = 3)
ax.loglog(f_res,np.max(TF_spectrum),'d',mfc = 'r',mec = 'k',zorder = 4)
ax.vlines(f_res,3e-5,1e-2,colors = 'r',linestyle = '--',zorder = 2)
ax.set_xlim(1e-1,1e2)
ax.set_ylim(3e-5,1e-2)
ax.grid(which = 'major',linestyle = '-',zorder = 1)

ax.set_xlabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
ax.set_ylabel(r'$\langle \hat{\zeta} \rangle_{x,y} (f)$',labelpad = 5)

#%% Merging Figure 2 

# TF spectrum
fig, ax = plt.subplots(2,2,figsize = (12,9),layout = 'constrained')
spectrum_window = (1,0)
current_ax = ax[spectrum_window]
current_ax.loglog(freq,abs(TF_spectrum),'o',zorder = 3)
current_ax.loglog(f_res,np.max(TF_spectrum),'d',ms = 8, mfc = 'r',mec = 'k',zorder = 4)
current_ax.vlines(f_res,3e-5,1e-2,colors = 'r',linestyle = '--',lw = 2, zorder = 2)
current_ax.set_xlim(1e-1,1e2)
current_ax.set_ylim(3e-5,1e-2)
current_ax.grid(which = 'major',linestyle = '-',zorder = 1)

current_ax.set_xlabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)
current_ax.set_ylabel(r'$\langle \hat{\zeta} \rangle_{x,y} (f)$',labelpad = 5)

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

cbar = plt.colorbar(wavefield, ax = current_ax,shrink = 0.79)
cbar.set_label(r'$\zeta \; \mathrm{(mm)}$',labelpad = 5)
# cbar.set_ticks([-1,0,1])
# cbar.set_ticklabels([f'{l:.1e}' for l in [-5e-2,0,5e-2]])

# Resonance frequencies 
resonance_window = (1,1)
current_ax = ax[resonance_window]

bounds = np.array([7,8,12,18,22,38])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)


scatter = current_ax.scatter(S_center['H'],S_center['f'], c = S_center['R'], s = 70,cmap = new_blues,
           norm = norm,edgecolors = 'k', zorder = 2)
# ax.grid(zorder = 1)

current_ax.set_xlim(3,22)
current_ax.set_ylim(2.5,6)
cbar = plt.colorbar(scatter)
cbar.set_label(r'$R \; \mathrm{(mm)}$',labelpad = 5)

current_ax.set_xlabel(r'$h \; \mathrm{(mm)}$',labelpad = 5)
current_ax.set_ylabel(r'$f_0 \; \mathrm{(Hz)}$',labelpad = 5)  

midpoints = [(bounds[i] + bounds[i+1])/2 for i in range(len(bounds) - 1)]
cbar.set_ticks(midpoints)
cbar.set_ticklabels([f'{mid:.1f}' for mid in midpoints])


figname = fig_folder + 'Figure2_python_subplots'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')


