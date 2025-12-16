# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 18:53:07 2025

@author: sebas
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 

import h5py

import scipy 
import pywt


import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.tools.Fourier_tools as FT
import icewave.sebastien.set_graphs as set_graphs
import icewave.das.DAS_package as DS
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather

# PARULA COLORMAP 
parula_map = matcmaps.parula()

full_yarg = mpl.colormaps['gist_yarg'].resampled(256)
new_yarg = colors.ListedColormap(full_yarg(np.linspace(0.1,1,256)))

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))


global g
g = 9.81
global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/spatio/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% FUNCTION SECTION 

def plot_spatio_temp(spatio,fiber_length,extents,cmap):
    """ Plot spatio-temporal using specific format
    Inputs: - spatio, numpy 2D array [nt,nx],
            - fiber_length, float, length of fiber, as set in Febus software
    Outputs: - fig, matplotlib figure
             - ax, matplotlib axis object
             - imsh, matplotlib imshow object
             - cbar, matplotlib colorbar object """
    
    
    normalization = 'linear'
    fig,ax = plt.subplots()
    imsh = ax.imshow(spatio.T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = cmap,
              interpolation = 'gaussian', extent = extents)
    ax.set_ylim([0,fiber_length])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)

    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
    
    return fig,ax,imsh,cbar

def annotate_axes(ax, text, fontsize=18):
    ax.text(0.5, 0.5, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize, color="darkgrey")

#%% Load DAS parameters and data

# Load parameters for DAS
main_path = 'U:/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

date = '0211'
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 1 #5 for 0211
file2load = filelist[idx_file]
print(file2load)

Nb_minutes = 10 # duration of each stack
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
format_date = '%Y-%m-%dT%H-%M-%S'
label_UTC0 = UTC_stack[0,0].strftime(format_date)

# shift curvilinear axis
offset_fiber = 37.5
s = s - offset_fiber
# decimation for 0210 and 0212
# if date in date_downsampling:
#     fs = fs/down_sampling_factor # new value of sampling frequency
#     stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
#                                                                         stack_time,UTC_stack,down_sampling_factor)
#     print(f'New value of sampling frequency, fs = {fs:.2f}')

#%% Show spatio-temporal 

chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,6))
imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax.set_ylim([0,s[-1]])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'UTC')

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

figname = f'{fig_folder}spatio_temporal_{label_UTC0}_chunk_{chunk}_general'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Show zoomed on ZEN spatio-temporal

chunk = 0
set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,6))
imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax.set_ylim([0,s[-1]])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'UTC')
xmin = datetime(2025,2,11,18,43,48)
xmax= datetime(2025,2,11,18,44,40)
ax.set_xlim([xmin,xmax])
ax.set_ylim([0,160])

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

figname = f'{fig_folder}spatio_temporal_{label_UTC0}_chunk_{chunk}_ZEN'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Show zoomed on flexural wave spatio-temporal

chunk = 0
xmin = datetime(2025,2,11,18,43,54,110000)
xmax= datetime(2025,2,11,18,43,55,500000)
t0 = xmin.timestamp()
t1 = xmax.timestamp()

# fs = 84
# t = np.linspace(0,t1-t0,round((t1-t0)*fs))

# timestamp_array = np.array([UTC.timestamp() for UTC in UTC_stack[chunk,:]])
# mask = np.logical_and(timestamp_array > t0, timestamp_array < t1)

set_graphs.set_matplotlib_param('single')
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
# fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
normalization = 'linear'
fig,ax = plt.subplots(figsize = (12,6))
imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
          interpolation = 'gaussian', extent = extents)
ax.set_ylim([0,fiber_length])

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(imsh,cax = cax)

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'$t \; \mathrm{(s)}$')

ax.set_xlim([xmin,xmax])
ax.set_ylim([0,160])

# change tick labels
xticks = ax.get_xticks()
t = np.linspace(0,t1-t0,len(xticks))
xlabels = [r'$' f'{element:.2f}' '$' for element in t]
print(xlabels)
ax.set_xticks(xticks,labels = xlabels)

cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

offset_text = cbar.ax.yaxis.get_offset_text()
offset_text.set_x(1)

figname = f'{fig_folder}spatio_temporal_{label_UTC0}_chunk_{chunk}_Z_zoomed'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Create a subplot of three graphs

set_graphs.set_matplotlib_param('single')
fig = plt.figure(figsize = (12,8),constrained_layout = True)
gs_kw = dict(width_ratios = [1,1],height_ratios = [1.4,1])
axd = fig.subplot_mosaic([['top','top'],['bottom_left', 'bottom_right']],gridspec_kw = gs_kw)

# for k, ax in axd.items():
#     annotate_axes(ax, f'axd[{k!r}]', fontsize=14)
    
chunk = 0
extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
normalization = 'linear'

for ax_key in axd.keys():
    if ax_key == 'top':
        ax = axd[ax_key]
        imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
                  interpolation = 'gaussian', extent = extents)
        ax.set_ylim([0,s[-1]])
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.1)
        cbar = plt.colorbar(imsh,cax = cax)
        
        ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
        ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
        
        imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
        ax.set_xlabel(r'UTC')
        
        cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(a.u.)}$')
        cbar.formatter.set_powerlimits((3, 3))
        cbar.update_ticks()
        
        # offset_text = cbar.ax.yaxis.get_offset_text()
        # offset_text.set_x(1)

    elif ax_key == 'bottom_left':
        
        ax = axd[ax_key]
        imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
                  interpolation = 'gaussian', extent = extents)
        ax.set_ylim([0,fiber_length])

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="2%", pad=0.1)
        # cbar = plt.colorbar(imsh,cax = cax)
        # cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
        # cbar.formatter.set_powerlimits((3, 3))
        # cbar.update_ticks()

        ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
        ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

        imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
        ax.set_xlabel(r'UTC')
        xmin = datetime(2025,2,11,18,43,48)
        xmax= datetime(2025,2,11,18,44,40)
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([0,160])

        # offset_text = cbar.ax.yaxis.get_offset_text()
        # offset_text.set_x(1)
        
    else :
        
        ax = axd[ax_key]
        xmin = datetime(2025,2,11,18,43,54,110000)
        xmax= datetime(2025,2,11,18,43,55,500000)
        t0 = xmin.timestamp()
        t1 = xmax.timestamp()

        # fs = 84
        # t = np.linspace(0,t1-t0,round((t1-t0)*fs))

        # timestamp_array = np.array([UTC.timestamp() for UTC in UTC_stack[chunk,:]])
        # mask = np.logical_and(timestamp_array > t0, timestamp_array < t1)
        imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
                  interpolation = 'gaussian', extent = extents)
        ax.set_ylim([0,s[-1]])

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="2%", pad=0.1)
        # cbar = plt.colorbar(imsh,cax = cax)
        # cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
        # cbar.formatter.set_powerlimits((3, 3))
        # cbar.update_ticks()

        ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
        ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

        imsh.set_clim([-1e4,1e4]) # [-1e4,1e4] for 0211
        ax.set_xlabel(r'$t \; \mathrm{(s)}$')

        ax.set_xlim([xmin,xmax])
        ax.set_ylim([0,160])

        # change tick labels
        xticks = ax.get_xticks()
        t = np.linspace(0,t1-t0,len(xticks))
        xlabels = [r'$' f'{element:.2f}' '$' for element in t]
        print(xlabels)
        ax.set_xticks(xticks,labels = xlabels)

        # offset_text = cbar.ax.yaxis.get_offset_text()
        # offset_text.set_x(1)


figname = f'{fig_folder}montage_spatio_temporal_{label_UTC0}_chunk_{chunk}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')