# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:13:25 2026

@author: sebas
"""


import os
import re
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from datetime import datetime, time , timedelta
import pytz
import glob 
import cv2 as cv
import h5py
import scipy
import csv

import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.field.drone as field_drone
import icewave.field.gps as field_gps
import icewave.gps.gps as gps_pack

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% Set fig_folder path 
disk = 'F:/'
fig_folder = f'{disk}Rimouski_2025/DAS_article/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Load DAS results 
base = f'{disk}Rimouski_2025/Data/'

file2load = f'{base}/Summary/DAS/main_results_active_passive_V3.h5'
results_DAS = rw.load_dict_from_h5(file2load)

# folder2active = 'F:/Rimouski_2025/DAS_article/Uncertainties/'
# file2load = f'{folder2active}active_E_modulus_results_uncertainties.pkl'
# with open(file2load,'rb') as pf:
#     E_results = pickle.load(pf)
    
# file2load = f'{folder2active}active_thickness_results_uncertainties.pkl'
# with open(file2load,'rb') as pf:
#     h_results = pickle.load(pf)

# file2load = f'{folder2active}active_flex_modulus_results_uncertainties.pkl'
# with open(file2load,'rb') as pf:
#     D_results = pickle.load(pf)

#%% Load DAS GPS position 

file_water_height = f'{base}/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)
    
#%% Load geophones acquisition position

dates = ['0210','0211','0212']

GPS = {'Geophones':{},'thickness':{}}
date = '0210'
year = '2025'

path2logfile = f'{base}{date}/Geophones'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
# geophones_table_path = 'C:/Users/sebas/Github/icewave/icewave/geophone/geophones_table'
UTC_timearea = pytz.timezone('UTC')

# load geophone GPS coordinates 
new_matrix = geophone_gps.get_GPS_coordinates(path2logfile, geophones_table_path, date, year)
GPS_geo = geophone_gps.build_GPS_array(new_matrix)
# get GPS coordinates for each acquisition  
geoph_results = {'acq0001':{},'acq0002':{},'acq0003':{}}
for idx_acq,key in enumerate(list(geoph_results.keys())):
    mean_GPS = np.mean(GPS_geo[:,idx_acq,:],axis = 0)
    geoph_results[key]['lat'] = mean_GPS[1]
    geoph_results[key]['long'] = mean_GPS[0]
    
#%% Collect geophones results from MCMC inversion 

date = '0210'
path2geoph = f'{base}{date}/Geophones/'
filelist = glob.glob(f'{path2geoph}*distribution_results.pkl')
keys = list(geoph_results.keys())
for i,file2load in enumerate(filelist):
    with open(file2load,'rb') as pf:
        test = pickle.load(pf)
    key = keys[i]
    geoph_results[key].update(test)

for key in geoph_results.keys():
    for k in geoph_results[key]['E'].keys():
        geoph_results[key]['E'][k]  = geoph_results[key]['E'][k]*1e9 


#%% Load Maddie GPS results 

filename = f'{base}/Summary/DAS/Maddie_gps_measurements_DAS.pkl'
with open(filename,'rb') as pf:
    MS_gps = pickle.load(pf)
    
#%% Collect Stephane GPS waypoints 

path2records = f'{base}/Summary/records_gps_BicWin25.pkl'
with open(path2records,'rb') as pf:
    records = pickle.load(pf)
# records = {}
# for date in dates :
#     record = field_gps.get_records(date,year = '2025')
#     records[date] = record['gps']['garminSP']
  
# collect gps thickness measurements 
gps_results = {}
for date in records.keys():
    
    gps_results[date] = {'latitude':[],'longitude':[],'h':[],'time':[],'name':[],'elevation':[]}
    for key in records[date].keys():
        if 'H' in key :
            for key_property in records[date][key].keys():
                gps_results[date][key_property].append(records[date][key][key_property][0])
                
            gps_results[date]['h'].append(int(key[2:])/100)
            
    # add Maddie thickness waypoints
    if date in MS_gps:
        for key in MS_gps[date].keys():
            for elem in MS_gps[date][key]:
                gps_results[date][key].append(elem)
    
    for key in gps_results[date].keys():
        gps_results[date][key] = np.array(gps_results[date][key])


# =============================================================================
#%% Compute Flexural modulus for gps positions and geophones lines 
# =============================================================================

# define initial position of fiber 
Lat0 = DAS_water_height['lat'][0]
Long0 = DAS_water_height['long'][0]

# geophone lines
for key_acq in geoph_results.keys():
    E = geoph_results[key_acq]['E']['max']
    std_E = geoph_results[key_acq]['E']['std']
    h = geoph_results[key_acq]['h']['max']
    std_h = geoph_results[key_acq]['h']['std']
    nu = geoph_results[key_acq]['nu']['max']
    std_nu = geoph_results[key_acq]['nu']['std']   
    
    geoph_results[key_acq]['D'] = {}
    geoph_results[key_acq]['D']['max'] = E*h**3/12/(1-nu**2)
    dDdE = h**3/12/(1-nu**2)
    dDdh = 3*E*h**2/12/(1-nu**2)
    dDdnu = E*2*nu*h**3/12/((1-nu**2)**2)
    geoph_results[key_acq]['D']['std'] = np.sqrt(dDdE**2 * std_E**2 + dDdh**2 * std_h**2 + dDdnu**2 * std_nu**2)
    geoph_results[key_acq]['x'] = gps_pack.distance_GPS(geoph_results[key_acq]['lat'], geoph_results[key_acq]['long'], 
                                             Lat0, Long0)
    
# geophone lines
# model_key = '20250210'
# sub_results[model_key]['D'] = sub_results[model_key]['Young modulus (Pa)']*sub_results[model_key]['thickness (m)']**3/12/(1-
#                                                                         sub_results[model_key]['Poisson ratio']**2)

# sub_results[model_key]['x'] = gps_pack.distance_GPS(sub_results[model_key]['latitude'], sub_results[model_key]['longitude'], 
#                                          Lat0, Long0)

# gps thickness measurements 
E = np.mean(np.array([geoph_results[key]['E']['max'] for key in geoph_results.keys()]))
nu = np.mean(np.array([geoph_results[key]['nu']['max'] for key in geoph_results.keys()]))

for key_date in gps_results.keys():
    gps_results[key_date]['D'] = E*gps_results[key_date]['h']**3/12/(1-nu**2)
    gps_results[key_date]['x'] = gps_pack.distance_GPS(gps_results[key_date]['latitude'], gps_results[key_date]['longitude'], 
                                             Lat0, Long0)
    
#%% Plot ice thickness drilling measurements as mean value 

# create an array of all measured thickness
drillings = {'x':np.concatenate([gps_results[date]['x'] for date in gps_results.keys()]),
             'h': np.concatenate([gps_results[date]['h'] for date in gps_results.keys()])} 

# create three different zones 
R_bounds = np.array([[0,170],[170,310],[310,600]])
data_h = {} 
regions_name = ['R1','R2','R3']
for i in range(R_bounds.shape[0]):
    
    xmin = R_bounds[i,0]
    xmax = R_bounds[i,1]
    
    mask = np.logical_and(drillings['x'] > xmin, drillings['x'] < xmax)
    data_h[regions_name[i]] = np.array([drillings['x'][mask],drillings['h'][mask]])
    
# create an array of three points (one for each region)
avg_data = np.zeros((len(data_h),2,2))
for i,key in enumerate(data_h.keys()):
    avg = np.mean(data_h[key],axis = 1)
    std = np.std(data_h[key],axis = 1)
    avg_data[i,0,:] = avg
    avg_data[i,1,:] = std    

#%% Create subplot
offset_fiber = 37.5

selected_acq = ['acq0001','acq0003']
x_err = 45 # in meter
markersize = 8
x_boarder = np.array([170,310])

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

full_reds = mpl.colormaps['Reds'].resampled(256)
new_reds = colors.ListedColormap(full_reds(np.linspace(0.2,1,256)))

color_date = {}
color_date['active'] = {'0210':new_reds(0.2),'0211':new_reds(0.4), '0212': new_reds(0.8)}
color_date['passive'] = {'0210':new_blues(0.2),'0211':new_blues(0.6), '0212': new_blues(0.9)}

set_graphs.set_matplotlib_param('triple') 
fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize = (16,6.5),layout = 'constrained')   

# plot E VS x
# axs[0].plot(results_DAS['active']['0212']['x'],results_DAS['active']['0212']['E']*1e-9,'-o',
#             color = 'k',mec = 'k')
for date in results_DAS['active'].keys():
    if date == '0211':
        continue
    else :
        x = results_DAS['active'][date]['x']
        E = results_DAS['active'][date]['E']
        E_unc = results_DAS['active'][date]['E_unc']
        axs[0].plot(x,E*1e-9,'-o',
                    color = color_date['active'][date],mec = 'k')
        axs[0].fill_between(x,
                            (E - E_unc)*1e-9,
                            (E + E_unc)*1e-9,
                            alpha = 0.2,
                    color = color_date['active'][date])

for i,key in enumerate(selected_acq):
    if i == 0:
        axs[0].errorbar(geoph_results[key]['x'],geoph_results[key]['E']['max']*1e-9,
                yerr = geoph_results[key]['E']['std']*1e-9,xerr = x_err,fmt = '^',
            color = 'tab:green',mec = 'k',ms = markersize,elinewidth = 2,label = 'geophones')
    else :
        axs[0].errorbar(geoph_results[key]['x'],geoph_results[key]['E']['max']*1e-9,
                yerr = geoph_results[key]['E']['std']*1e-9,xerr = x_err,fmt = '^',
            color = 'tab:green',mec = 'k',ms = markersize,elinewidth = 2)
        
axs[0].set_xlabel(r'$x \; \mathrm{(m)}$')
axs[0].set_ylabel(r'$E \; \mathrm{(GPa)}$')
axs[0].set_ylim([3.5,6])


# plot h VS x
# plot drilling
axs[1].errorbar(avg_data[:,0,0],avg_data[:,0,1],xerr = avg_data[:,1,0],yerr = avg_data[:,1,1],
           fmt = 'none',ecolor = 'tab:cyan',elinewidth = 2,label = 'drilling',zorder = 3)

# for i,date in enumerate(gps_results.keys()):
#     if i == 0:
#         axs[1].plot(gps_results[date]['x'],gps_results[date]['h'],'p',
#                 color = 'tab:cyan',mec = 'k',ms = markersize,label = 'drilling',zorder = 2)
#     else:
#         axs[1].plot(gps_results[date]['x'],gps_results[date]['h'],'p',
#                 color = 'tab:cyan',mec = 'k',ms = markersize)
# plot DAS
# for date in results_DAS['active'].keys():
#     axs[1].plot(results_DAS['active'][date]['x'],results_DAS['active'][date]['h'],'-o',
#                 color = color_date['active'][date],mec = 'k',label = f'active {date}',zorder = 0)
for date in results_DAS['active'].keys():
    x = results_DAS['active'][date]['x']
    h = results_DAS['active'][date]['h']
    h_low = results_DAS['active'][date]['h_low']
    h_high = results_DAS['active'][date]['h_high']
    axs[1].plot(x,h,'-o',
                color = color_date['active'][date],mec = 'k',label = f'active {date}')
    axs[1].fill_between(x,
                        h - h_low,
                        h + h_high,
                        alpha = 0.2,
                color = color_date['active'][date])

# plot geophone
for key in selected_acq:
    axs[1].errorbar(geoph_results[key]['x'],geoph_results[key]['h']['max'],fmt = '^',
            yerr = geoph_results[key]['h']['std'],xerr = x_err, 
            color = 'tab:green',mec = 'k',ms = markersize,elinewidth = 2,zorder = 3)
# axs[1].plot(np.mean(sub_results[model_key]['x'][1:3]),np.mean(sub_results[model_key]['thickness (m)'][1:3]),'^',
#             color = 'tab:green',mec = 'k',ms = markersize)

axs[1].set_xlabel(r'$x \; \mathrm{(m)}$')
axs[1].set_ylabel(r'$h \; \mathrm{(m)}$')
axs[1].set_ylim([0.22,0.82])

# add vertical lines 


# plot D VS x
# plot DAS active
# axs[2].plot(results_DAS['active']['0212']['x'],results_DAS['active']['0212']['D'],'o-',
#             color = color_date['active']['0212'],mec = 'k')
for date in results_DAS['active'].keys():
    x = results_DAS['active'][date]['x']
    D = results_DAS['active'][date]['D']
    D_unc = results_DAS['active'][date]['D_unc']
    if date == '0211':
        D_unc[-1] = 37*1e6 # ensured no negative value of D - D_unc
    
    axs[2].plot(x,D,'-o',
                color = color_date['active'][date],mec = 'k')
    axs[2].fill_between(x,
                        D - D_unc,
                        D + D_unc,
                        alpha = 0.2,
                color = color_date['active'][date])
    
# plot DAS passive
for date in results_DAS['passive']['corrected'].keys():
    x = results_DAS['passive']['corrected'][date]['x']
    D = results_DAS['passive']['corrected'][date]['D']
    D_err = results_DAS['passive']['corrected'][date]['err_D']
    axs[2].plot(x,D,'o-',
                    color = color_date['passive'][date],mec = 'k',label = f'passive {date}')
    axs[2].fill_between(x,D - D_err,D + D_err,alpha= 0.2,color = color_date['passive'][date])
    
# plot geophones
for key in selected_acq:
    axs[2].errorbar(geoph_results[key]['x'],geoph_results[key]['D']['max'],
                    yerr = geoph_results[key]['D']['std'],xerr = x_err, fmt = '^',
            color = 'tab:green',mec = 'k',ms = markersize,elinewidth = 2)
# axs[2].plot(np.mean(sub_results[model_key]['x'][1:3]),np.mean(sub_results[model_key]['D'][1:3]),'^',
#             color = 'tab:green',mec = 'k',ms = markersize)

axs[2].set_xlabel(r'$x \; \mathrm{(m)}$')
axs[2].set_ylabel(r'$D \; \mathrm{(Pa.m^{3})}$')
axs[2].set_ylim([-1.1e6,2.6e8])
# axs[2].set_yscale('log')
# axs[2].set_ylim([3e6,4e8])

# add vertical lines 
axs[0].vlines(x_boarder,1, 10,linestyles = '--',colors = 'k',lw = 1)
axs[1].vlines(x_boarder,0, 10,linestyles = '--',colors = 'k',lw = 1)
axs[2].vlines(x_boarder,0, 1e9,linestyles = '--',colors = 'k',lw = 1)


for ax in axs:
    ax.set_xlim([0,600])

# set legend()
fig.legend(ncols = 3,
              loc='outside upper center',frameon = False)

figname = f'{fig_folder}Subplot_EhD_VS_x_with_uncertainties'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')
 

#%% compute average value of E within each region 

# R1
date = '0212'
x = results_DAS['active'][date]['x']
E = results_DAS['active'][date]['E']
mask = np.where(x < x_boarder[0])[0]
E_mean = np.mean(E[mask])
E_std = np.std(E[mask])

print(f'E(R1) = {E_mean*1e-9:.2f} ± {E_std*1e-9:.3f}')

#R2 
mask = np.where(np.logical_and(x < x_boarder[1],x > x_boarder[0]))[0]
E_mean = np.mean(E[mask])
E_std = np.std(E[mask])

print(f'E(R2) = {E_mean*1e-9:.2f} ± {E_std*1e-9:.3f}')

# R3
mask = np.where(x > x_boarder[1])[0]
E_mean = np.mean(E[mask])
E_std = np.std(E[mask])
print(f'E(R3) = {E_mean*1e-9:.2f} ± {E_std*1e-9:.3f}')

#%% Compute average value of h within each region

# R1
h_list = []
for date in results_DAS['active'].keys():
    x = results_DAS['active'][date]['x']
    h = results_DAS['active'][date]['h']
    mask = np.where(x < x_boarder[0])[0]
    for elem in h[mask]:
        h_list.append(elem)
    
h_list = np.array(h_list)
h_mean = np.mean(h_list)
h_std = np.std(h_list)
print(f'h(R1) = {h_mean:.2f} ± {h_std:.3f}')

# R2
h_list = []
for date in results_DAS['active'].keys():
    x = results_DAS['active'][date]['x']
    h = results_DAS['active'][date]['h']
    mask = np.where(np.logical_and(x < x_boarder[1],x > x_boarder[0]))[0]
    for elem in h[mask]:
        h_list.append(elem)
    
h_list = np.array(h_list)
h_mean = np.mean(h_list)
h_std = np.std(h_list)
print(f'h(R2) = {h_mean:.2f} ± {h_std:.3f}')

# R3
h_list = []
for date in results_DAS['active'].keys():
    x = results_DAS['active'][date]['x']
    h = results_DAS['active'][date]['h']
    mask = np.where(x > x_boarder[1])[0]
    for elem in h[mask]:
        h_list.append(elem)
    
h_list = np.array(h_list)
h_mean = np.mean(h_list)
h_std = np.std(h_list)
print(f'h(R3) = {h_mean:.2f} ± {h_std:.3f}')


#%% Compute average value of D within each region with active method

method = 'passive'
# R1
D_list = []
if method == 'active':
    main_dict = results_DAS[method]
else: 
    main_dict = results_DAS[method]['corrected']
for date in main_dict.keys():
    x = main_dict[date]['x']
    D = main_dict[date]['D']
    mask = np.where(x < x_boarder[0])[0]
    for elem in D[mask]:
        D_list.append(elem)
    
D_list = np.array(D_list)
D_mean = np.nanmean(D_list)
D_std = np.nanstd(D_list)
print(f'D(R1) = {D_mean*1e-6:.2f} ± {D_std*1e-6:.2f}')

# R2
D_list = []
if method == 'active':
    main_dict = results_DAS[method]
else: 
    main_dict = results_DAS[method]['corrected']
for date in main_dict.keys():
    x = main_dict[date]['x']
    D = main_dict[date]['D']
    mask = np.where(np.logical_and(x < x_boarder[1],x > x_boarder[0]))[0]
    for elem in D[mask]:
        D_list.append(elem)
    
D_list = np.array(D_list)
D_mean = np.nanmean(D_list)
D_std = np.nanstd(D_list)
print(f'D(R2) = {D_mean*1e-6:.2f} ± {D_std*1e-6:.2f}')

# R3
D_list = []
if method == 'active':
    main_dict = results_DAS[method]
else: 
    main_dict = results_DAS[method]['corrected']
for date in main_dict.keys():
    x = main_dict[date]['x']
    D = main_dict[date]['D']
    mask = np.where(x > x_boarder[1])[0]
    for elem in D[mask]:
        D_list.append(elem)
    
D_list = np.array(D_list)
D_mean = np.nanmean(D_list)
D_std = np.nanstd(D_list)
print(f'D(R3) = {D_mean*1e-6:.2f} ± {D_std*1e-6:.2f}')


#%% Estimate relative error on h when comparing active and passive methods 

E = 5.1*1e9
nu = 0.3
D = 62.6*1e6
h = (12*(1-nu**2)*D/E)**(1/3)
print(h)

# =============================================================================
#%% Archive - tests
# =============================================================================

#%% Collect geophone line results 

geoph_results = {}
path2csv = 'U:/General/Summary_geophone_lines/results_acquisitions.csv'
with open(path2csv, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        for key in row:
            if key not in geoph_results.keys():
                geoph_results[key] = [row[key]]
            else:
                geoph_results[key].append(row[key])

# convert strings to float
keys2convert = ['longitude','latitude','water level (m)','thickness (m)','Young modulus (Pa)',
                'Poisson ratio','density (kg/m^3)']
for key in keys2convert:
    current_list = []
    for value in geoph_results[key]:
        if value == '':
            current_list.append(None)
        else:
            current_list.append(float(value))
        
    geoph_results[key] = np.array(current_list)

#%% create sub dictionnary

sub_results = {}
date = '0210'
model_key = f'{year}{date}'
mask = [date_key == model_key for date_key in geoph_results['date']]
# mask = np.where(geoph_results['date'] == model_key)[0]

sub_results[model_key] = {}
for key in geoph_results.keys():
    sub_results[model_key][key] = np.array([value for i,value in enumerate(geoph_results[key]) if mask[i]])







