# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 15:49:58 2025

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
import imageio as iio
import cv2 as cv
import h5py
import csv
import scipy

# module for map plotting
import cmocean.cm as cmo
import alphashape
from shapely.geometry import MultiPoint
from shapely.geometry import Point, Polygon
from scipy.interpolate import griddata
import cartopy.io.shapereader as shpreader

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from shapely.ops import unary_union
import cartopy.io.shapereader as shpreader

# icewave modules
import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps
import icewave.das.DAS_package as DS


import icewave.analysis.bathy as bathy
import stephane.display.graphes as graphes
import icewave.gps.gps as gps
import icewave.geometry.tables as tables

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

global date_downsampling 
date_downsampling = ['0210','0212']
global down_sampling_factor
down_sampling_factor = 10

#%% Set fig_folder path

fig_folder = 'U:/Data/Summary/DAS/Sources_location/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% FUNCTION SECTION

def get_UTC_from_filename(filename):

    array_str = filename.split('\\')
    date_str = array_str[-1][1:-7]
    print(date_str)

    # create datetime object
    format_date = '%Y-%m-%d_%H-%M-%S'
    UTC_0 = datetime.strptime(date_str,format_date)
    
    return UTC_0

#-------------------------------------------------------------------------

def compute_source_UTC(sources_xt,UTC_array):
    """ Compute UTC time of each source, based on initial UTC time of each acquisition file
    and the time tmin at which tyhe source is observed """
    
    for i in range(len(sources_xt)):
        source = sources_xt[i]
        dt = timedelta(seconds = source['tmin'])
        UTC_source = UTC_array[source['acq_num']] + dt
        sources_xt[i]['UTC'] = UTC_source.replace(tzinfo = pytz.utc)
        
    return sources_xt

#------------------------------------------------------------------------

def linear(x,a):
    y = a*x
    return y

#--------------------------------------------------------------------------

def affine(x,a,b):
    y = a*x + b
    return y

#%% Load sources GPS position and time from Stephane's GPS
 
data_gps = {}
dates = ['0211','0212']
for date in dates:
    data_gps[date] = {}
    
    path2data = f'U:/Data/{date}/GPS/'
    gpx_filelist = glob.glob(path2data + '*.gpx')
    
    gpx_file = gpx_filelist[0]
    # get gpx
    gpx = gps.get_wpts(gpx_file)
    
    # get waypoints position
    # use Map_table to create a dictionnary
    data_gps[date] = tables.dict_from_gpx(gpx,path2data)
    
# keep only waypoints related to sources and borne
sources_gps = {}
for date in data_gps.keys():
    sources_gps[date] = {}
    for key in data_gps[date].keys():
        instrument = key[0]
        
        if instrument == 'S':
            sources_gps[date][key] = data_gps[date][key]
            
        elif key == 'borne_01':
            sources_gps[date][key] = data_gps[date][key]
            
# adjust UTC time of GPS sources
for date in sources_gps.keys():
    for key in sources_gps[date].keys():
        sources_gps[date][key]['time'] = sources_gps[date][key]['time'].astimezone(pytz.utc)
 
#%% Ludovic's sources detection 

sources_xt = {}
sources_xt['0212'] = [
    {'acq_num': 0, 'x0': 58, 'tmin': 164.30},
    {'acq_num': 0, 'x0': 85, 'tmin': 307.29},
    {'acq_num': 1, 'x0': 130, 'tmin': 53.68},
    {'acq_num': 1, 'x0': 152, 'tmin': 135.12},
    {'acq_num': 1, 'x0': 170, 'tmin': 235.50},
    {'acq_num': 1, 'x0': 190, 'tmin': 327.4},
    {'acq_num': 1, 'x0': 211, 'tmin': 404.57},
    {'acq_num': 1, 'x0': 231, 'tmin': 492.65},
    {'acq_num': 2, 'x0': 260, 'tmin': 129.66},
    {'acq_num': 2, 'x0': 290, 'tmin': 308.40},#298.80 #332.70
    {'acq_num': 2, 'x0': 327, 'tmin': 414.66},
    {'acq_num': 2, 'x0': 348, 'tmin': 495.54},
    {'acq_num': 3, 'x0': 387, 'tmin': 71.68},#89.34
    {'acq_num': 3, 'x0': 410, 'tmin': 146.79},
    {'acq_num': 3, 'x0': 426, 'tmin': 230.82},
    {'acq_num': 3, 'x0': 442, 'tmin': 310.76},
    {'acq_num': 4, 'x0': 460, 'tmin': 185.17},
    {'acq_num': 4, 'x0': 490, 'tmin': 321.03},
    {'acq_num': 4, 'x0': 505, 'tmin': 563.85},
    {'acq_num': 5, 'x0': 525, 'tmin': 46.18},
    {'acq_num': 5, 'x0': 545, 'tmin': 122.85},
    {'acq_num': 5, 'x0': 565, 'tmin': 200.61},
    {'acq_num': 5, 'x0': 585, 'tmin': 322.68},
]

sources_xt['0211'] = [
     {'acq_num': 0, 'x0': 40, 'tmin': 125.7},
     {'acq_num': 0, 'x0': 60, 'tmin': 356.1},
     {'acq_num': 0, 'x0': 90, 'tmin': 494.4},
     {'acq_num': 0, 'x0': 120, 'tmin': 560.5},
     {'acq_num': 1, 'x0': 150, 'tmin': 245.1},
     {'acq_num': 1, 'x0': 170, 'tmin': 329.6},
     {'acq_num': 2, 'x0': 190, 'tmin': 9.2},        #  NOT GOOD
     {'acq_num': 2, 'x0': 290, 'tmin': 453.4},      #  NOT GOOD
     {'acq_num': 2, 'x0': 320, 'tmin': 573.7},      #  NOT GOOD
     {'acq_num': 3, 'x0': 330, 'tmin': 16.5},       #  NOT GOOD
     {'acq_num': 3, 'x0': 365, 'tmin': 189.7}, 
     {'acq_num': 3, 'x0': 385, 'tmin': 303.7}, 
     {'acq_num': 3, 'x0': 405, 'tmin': 382.4}, 
     {'acq_num': 3, 'x0': 420, 'tmin': 468.9}, 
     {'acq_num': 3, 'x0': 445, 'tmin': 549}, 
     {'acq_num': 4, 'x0': 485, 'tmin': 106}, 
     {'acq_num': 4, 'x0': 505, 'tmin': 183},
     {'acq_num': 4, 'x0': 525, 'tmin': 265.6},]

# save sources_xt structure
file2save = f'{fig_folder}Active_sources_spatio_detection.h5'
rw.save_dict_to_h5(sources_xt, file2save)


#%% Compute sources_xt UTC time from UTC file time 

base = 'U:/Data/'

for date in ['0211','0212']:
    path2DAS_data = f'{base}{date}/DAS/*_UTC.h5'
    filelist = glob.glob(path2DAS_data)
    
    # get UTC initial time for each acquisition 
    UTC_array = np.array([get_UTC_from_filename(filename) for filename in filelist])
    
    # Compute UTC time from tmin
    sources_xt[date] = compute_source_UTC(sources_xt[date],UTC_array)

#%% Plot all sources Stephane GPS in GPS coordinates 

Lat0 = sources_gps['0212']['borne_01']['latitude']
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

color_dict = {'0211':'tab:blue','0212':'tab:orange','0210':'tab:green'}

for date in ['0211','0212']:
    current_color = color_dict[date]
    for key in sources_gps[date].keys():
        if key != 'borne_01':
            ax.plot(sources_gps[date][key]['longitude'],sources_gps[date][key]['latitude'],'o',color = current_color,
                    mec = 'k')
        else : 
            ax.plot(sources_gps[date][key]['longitude'],sources_gps[date][key]['latitude'],'s',color = current_color,
                    mec = 'k')
            
# ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x
ax.set_xlabel('Longitude (°)')
ax.set_ylabel('Latitude (°)')

# figname = f'{fig_folder}GPS_sources_map'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Plot all sources Stephane GPS in metric coordinates


# compute azimuth between first and last point of DAS waypoints 
Lat0 = sources_gps['0212']['borne_01']['latitude']
Long0 = sources_gps['0212']['borne_01']['longitude']

Lat1 = sources_gps['0211']['S_32']['latitude']
Long1 = sources_gps['0211']['S_32']['longitude']
psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(Long1 - Long0)*np.pi/180,(Lat1 - Lat0)*np.pi/180)*180/np.pi

for date in sources_gps.keys():
    for key in sources_gps[date].keys():
        lat = sources_gps[date][key]['latitude']
        long = sources_gps[date][key]['longitude']
        psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(long - Long0)*np.pi/180,(lat - Lat0)*np.pi/180)*180/np.pi
        X,Y = dp.GPS2XY(sources_gps[date][key]['latitude'],sources_gps[date][key]['longitude'], Lat0, Long0, psi)
        sources_gps[date][key]['X'] = X
        sources_gps[date][key]['Y'] = Y


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

color_dict = {'0211':'tab:blue','0212':'tab:orange','0210':'tab:green'}

for date in ['0211','0212']:
    current_color = color_dict[date]
    for key in sources_gps[date].keys():
        if key != 'borne_01':
            ax.plot(sources_gps[date][key]['X'],sources_gps[date][key]['Y'],'o',color = current_color,
                    mec = 'k')
        else : 
            ax.plot(sources_gps[date][key]['X'],sources_gps[date][key]['Y'],'s',color = current_color,
                    mec = 'k')
                

#%% Compute x distance of each source to borne_01

Lat0 = np.mean([sources_gps[date]['borne_01']['latitude'] for date in ['0212']])
Long0 = np.mean([sources_gps[date]['borne_01']['longitude'] for date in ['0212']])

for date in sources_gps.keys():
    for key in sources_gps[date].keys():
        sources_gps[date][key]['x'] = gps.distance_GPS(sources_gps[date][key]['latitude'],sources_gps[date][key]['longitude'],
                                                       Lat0,Long0)
       
    
#%% Plot x VS time for gps_sources and sources detection done by Ludovic

offset_fiber = 37.5 # amount of fiber enrolled in the DAS unit
txt_shift = 1

for date in sources_gps.keys():
    fig, ax = plt.subplots()
    
    current_color = color_dict[date]
    # plot Stephane GPS data
    x = np.array([sources_gps[date][key]['time'].timestamp() for key in sources_gps[date].keys()])
    y = np.array([sources_gps[date][key]['x'] for key in sources_gps[date].keys()])
    ax.plot(x,y,'o',color = current_color,mec = 'k',
            label = 'GPS')
    
    for i,key in enumerate(sources_gps[date].keys()):
        ax.text(x[i]*txt_shift,y[i]*txt_shift,key,fontsize = 8)
    
    # plot Ludovic's sources detection
    x = [sources_xt[date][i]['UTC'].timestamp() for i in range(len(sources_xt[date]))]
    y = [sources_xt[date][i]['x0'] - offset_fiber for i in range(len(sources_xt[date]))]
    ax.plot(x,y,'^',color = current_color,mec = 'k',
            label = 'Spatio detection')
    
    
    ax.legend()
    ax.set_xlabel(r'UTC')
    ax.set_ylabel(r'distance to IU')
    
    ax.set_title(date)
    ax.grid()

    # figname = f'{fig_folder}Comparison_GPS_sources_VS_spatio_sources_x_VS_UTC_{date}'
    # plt.savefig(figname + '.pdf', bbox_inches='tight')
    # plt.savefig(figname + '.svg', bbox_inches='tight')
    # plt.savefig(figname + '.png', bbox_inches='tight')
    
#%% Make correspondance between sources

xt2gps_table = {}
xt2gps_table['0211'] = ['borne_01','S_02','S_03','S_05','S_08','S_09','S_11','S_15','S_16','S_16','S_18','S_19','S_20',
                        'S_21','S_23','S_25','S_26','S_27']
xt2gps_table['0212'] = ['S_01','S_03','S_06','S_07','S_08','S_09','S_10','S_11',None,'S_12','S_13','S_14','S_16','S_17',
                        'S_18','S_20','S_21','S_23','S_24','S_25','S_26',None,None]

#%% Plot X vs X using correspondance 

offset_fiber = 39.5

fig, ax = plt.subplots()
for date in xt2gps_table.keys():

    relevant_keys = xt2gps_table[date]
    current_color = color_dict[date]
    y = np.array([sources_xt[date][i]['x0']- offset_fiber for i in range(len(sources_xt[date]))])
    x = []
    for key in relevant_keys:
        if key != None:
            x.append(sources_gps[date][key]['x'])
        else:
            x.append(None)
    x = np.array(x)
    ax.plot(x,y,'o',color = current_color,mec = 'k')
    
    x_th = np.linspace(-10,620,50)
    if offset_fiber is None:
        popt,pcov = scipy.optimize.curve_fit(lambda x,a,b : affine(x,a,b),x,y,nan_policy = 'omit')    
        coeff = popt
        err_coeff = np.sqrt(np.diag(pcov))
        y_th = affine(x_th,coeff[0],coeff[1])
        label_th = f'{date}: ' + '$' + f'y = {popt[0]:.3f}x + {popt[1]:.2f}' + '$'

    else :
    # fit by an affine function 
        popt,pcov = scipy.optimize.curve_fit(lambda x,a : linear(x,a),x,y,nan_policy = 'omit')
        coeff = popt[0]
        err_coeff = np.sqrt(np.diag(pcov[0]))
        y_th = linear(x_th,coeff)
        label_th = f'{date}: ' + '$' + f'y = {popt[0]:.3f}x' + '$'

    ax.plot(x_th,y_th,'--',color = current_color,label = label_th)

ax.set_xlabel(r'$x_{GPS} \; \mathrm{(m)}$')
ax.set_ylabel(r'$x_{spatio} \; \mathrm{(m)}$')

ax.set_xlim([-5,610])
ax.set_ylim([-5,610])

ax.legend()

offset_txt = f'{offset_fiber:.1f}'
offset_txt = offset_txt.replace('.','p')
figname = f'{fig_folder}Comparison_x_spatio_VS_x_gps_offset_{offset_txt}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Load a spatio 

# Load parameters for DAS
main_path = 'U:/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

colorscale = {'0211':1e4,'0212':5e3}

date = '0211'
# Load spatio-temporal 
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')
idx_file = 0
file2load = filelist[idx_file]
Nb_minutes = 10
stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)

# Show spatio-temporal 
markersize = 8

chunk = 0
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
cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
cbar.formatter.set_powerlimits((3, 3))
cbar.update_ticks()

ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)

imsh.set_clim([-colorscale[date],colorscale[date]]) # [-1e4,1e4] for 0211
ax.set_xlabel(r'UTC')
ax.legend()


#%% Superpose sources on spatio-temporal plots

# Load parameters for DAS
main_path = 'U:/'
path2DAS_param = f'{main_path}Data/parameters_Febus_2025.pkl'

colorscale = {'0211':1e4,'0212':5e3}

date = '0211'
# Load spatio-temporal 
fs,fiber_length,facq_x = DS.get_DAS_parameters(path2DAS_param,date)

# Load DAS data 
path2data = f'{main_path}Data/{date}/DAS/'
filelist = glob.glob(path2data + '*UTC.h5')

for idx_file in range(len(filelist)):
    # idx_file = 0 #5 for 0211
    file2load = filelist[idx_file]
    print(file2load)
    
    if date == '0212' and idx_file == 3:
        Nb_minutes = 9
    else:
        Nb_minutes = 10 # duration of each stack
    stack_strain,stack_time,UTC_stack,s = DS.stack_data_fromfile(file2load, fiber_length, Nb_minutes)
    format_date = '%Y-%m-%dT%H-%M-%S'
    label_UTC0 = UTC_stack[0,0].strftime(format_date)
    
    # decimation for 0210 and 0212
    if date in date_downsampling:
        fs = fs/down_sampling_factor # new value of sampling frequency
        stack_strain,stack_time,UTC_stack = DS.time_decimation_stack_strain(stack_strain,
                                                                            stack_time,UTC_stack,down_sampling_factor)
        print(f'New value of sampling frequency, fs = {fs:.2f}')
    
    # keep only sources related to this acquisition
    UTC_start = UTC_stack[0,0]
    UTC_end = UTC_stack[0,-1]

    # keep only sources within UTC range
    sub_sources_xt = [i for i in range(len(sources_xt[date])) if sources_xt[date][i]['acq_num'] == idx_file]
    sub_sources_gps = []
    for key in sources_gps[date].keys():
        test = np.logical_and(sources_gps[date][key]['time'] >= UTC_start, sources_gps[date][key]['time'] <= UTC_end)
        if test:
            sub_sources_gps.append(key)

    # Show spatio-temporal 
    markersize = 8
    
    chunk = 0
    set_graphs.set_matplotlib_param('single')
    extents = [UTC_stack[chunk,0],UTC_stack[chunk,-1],s[0],s[-1]]
    # fig, ax ,imsh, cbar = plot_spatio_temp(stack_strain[chunk,:,:], fiber_length, extents, 'seismic')
    normalization = 'linear'
    fig,ax = plt.subplots(figsize = (12,6))
    imsh = ax.imshow(stack_strain[chunk,:,:].T,origin = 'lower',aspect = 'auto',norm = normalization, cmap = 'seismic',
              interpolation = 'gaussian', extent = extents)
    
    # plot spatio sources
    current_color = 'tab:green'
    x = [sources_xt[date][idx]['UTC'] for idx in sub_sources_xt]
    y = [sources_xt[date][idx]['x0'] for idx in sub_sources_xt]
    ax.plot(x,y,'^',color = current_color,mec = 'k',ms = markersize,
            label = 'Spatio')
    
    # plot GPS sources
    x = [sources_gps[date][key]['time'] for key in sub_sources_gps]
    y = [sources_gps[date][key]['x'] + offset_fiber for key in sub_sources_gps]
    ax.plot(x,y,'o',color = current_color,mec = 'k',ms = markersize,
            label = 'GPS')
    
    ax.set_ylim([0,fiber_length])
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(imsh,cax = cax)
    cbar.set_label(r'$\dot{\epsilon} \; \mathrm{(u.a.)}$')
    cbar.formatter.set_powerlimits((3, 3))
    cbar.update_ticks()
    
    ax.set_xlabel(r'$t \; \mathrm{(s)}$',labelpad = 5)
    ax.set_ylabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
    
    imsh.set_clim([-colorscale[date],colorscale[date]]) # [-1e4,1e4] for 0211
    ax.set_xlabel(r'UTC')
    ax.legend()
    
    offset_text = cbar.ax.yaxis.get_offset_text()
    offset_text.set_x(1)
    
    figname = f'{fig_folder}Spatio_sources_comparison_{label_UTC0}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
    plt.close('all')



