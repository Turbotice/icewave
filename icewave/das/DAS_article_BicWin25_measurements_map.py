# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 16:43:10 2025

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
import imageio as iio
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

plt.rcParams.update({
    "text.usetex": True}) # use latex

full_blues = mpl.colormaps['Blues'].resampled(256)
new_blues = colors.ListedColormap(full_blues(np.linspace(0.2,1,256)))

#%% FUNCTION SECTION


#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/Measurement_map/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load DAS GPS coordinates

file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

# GPS coordiantes of fiber beginning
Lat0 = DAS_water_height['lat'][0]
Long0 = DAS_water_height['long'][0]

#%% Load DAS results 

file2load = 'U:/Data/Summary/DAS/main_results_active_passive.h5'
results_DAS = rw.load_dict_from_h5(file2load)

#%% Load geophones acquisition position

dates = ['0210','0211','0212']

GPS = {'Geophones':{},'thickness':{}}
date = '0210'
year = '2025'

path2logfile = f'U:/Data/{date}/Geophones'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
UTC_timearea = pytz.timezone('UTC')

# load geophone GPS coordinates 
new_matrix = geophone_gps.get_GPS_coordinates(path2logfile, geophones_table_path, date, year)

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
    
    
#%% Plot geophones results of 0210

cmap = new_blues
norm_thick = colors.Normalize(vmin = 0.3,vmax = 1.0)

fig, ax = plt.subplots()
model_key = '20250210'

current_color = cmap(norm_thick(sub_results[model_key]['thickness (m)']))
long = sub_results[model_key]['longitude']
lat = sub_results[model_key]['latitude']

ax.scatter(long,lat,marker = 's',color = current_color)


#%% Collect Stephane GPS waypoints 

records = {}
for date in dates :
    record = field_gps.get_records(date,year = '2025')
    records[date] = record['gps']['garminSP']
  
# collect gps thickness measurements 
gps_results = {}
for date in records.keys():
    
    gps_results[date] = {'latitude':[],'longitude':[],'h':[],'time':[],'name':[],'elevation':[]}
    for key in records[date].keys():
        if 'H' in key :
            for key_property in records[date][key].keys():
                gps_results[date][key_property].append(records[date][key][key_property][0])
                
            gps_results[date]['h'].append(int(key[2:])/100)
    
    for key in gps_results[date].keys():
        gps_results[date][key] = np.array(gps_results[date][key])

#%% Plot seismology deployment 

extents = [-68.82754,-68.81270,48.345778,48.348111]

ms = 60
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

# create norm for ice thickness
norm_acq = colors.Normalize(vmin = 1,vmax = 4)
full_cmap = mpl.colormaps['Oranges'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,1,256)))

color_date = {'0210':new_blues(0.2),'0211':new_blues(0.6), '0212': new_blues(0.9)}

# plot DAS positions
ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k')

# plot geophone lines and tomo
for geo_key in new_matrix.keys():
    for acqu_key in new_matrix[geo_key].keys():
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
        
        current_color = cmap(norm_acq(acqu_key))
        ax.scatter(GPS_logs['longitude'],GPS_logs['latitude'],marker = '^',color = current_color)

# plot ice thickness
# for key_date in records.keys():
#     for key in records[key_date].keys():
#         if 'H' in key:
#             h_value = int(key[2:])/100
#             current_color = cmap(norm_thick(h_value))
#             ax.scatter(records[key_date][key]['longitude'],records[key_date][key]['latitude'], marker = '^',
#                        color = current_color)

for key_date in gps_results.keys():
    current_color = color_date[key_date]
    ax.scatter(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],
               marker = 'p',color = current_color)

ax.set_xlim([extents[0], extents[1]])
ax.set_ylim([extents[2], extents[3]])

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

Lat0 = DAS_water_height['lat'][0]
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}seismology_measurements_map'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Superpose Geophones positions, DAS position and GPS positions using thickness as colorbar 

extents = [-68.82754,-68.81270,48.345778,48.348111]

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

# create norm for ice thickness
norm_thick = colors.Normalize(vmin = 0.2,vmax = 1.0)
cmap = new_blues

# plot DAS positions
ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k')

# # plot geophone lines and tomo
# for geo_key in new_matrix.keys():
#     for acqu_key in new_matrix[geo_key].keys():
#         GPS_logs = new_matrix[geo_key][acqu_key]
        
#         if  GPS_logs['num_GPS_logs'] > 1:
#             GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
            
#         ax.scatter(GPS_logs['longitude'],GPS_logs['latitude'])

# plot ice thickness
# for key_date in records.keys():
#     for key in records[key_date].keys():
#         if 'H' in key:
#             h_value = int(key[2:])/100
#             current_color = cmap(norm_thick(h_value))
#             ax.scatter(records[key_date][key]['longitude'],records[key_date][key]['latitude'], marker = '^',
#                        color = current_color)

for key_date in gps_results.keys():
    ax.scatter(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],c = gps_results[key_date]['h'],
               marker = '^',cmap = new_blues,norm = norm_thick)
    
# plot geophone thickness measurements
model_key = '20250210'
current_color = cmap(norm_thick(sub_results[model_key]['thickness (m)']))
long = sub_results[model_key]['longitude']
lat = sub_results[model_key]['latitude']

scatter = ax.scatter(long,lat,c = sub_results[model_key]['thickness (m)'],
                     marker = 's',cmap = new_blues,norm = norm_thick)

ax.scatter(results_DAS['active']['0212']['longitude'],results_DAS['active']['0212']['latitude'],
           c = results_DAS['active']['0212']['h'], cmap = new_blues,norm = norm_thick)

# colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(scatter,cax = cax)
cbar.set_label(r'$h \; \mathrm{(m)}$')

ax.set_xlim([extents[0], extents[1]])
ax.set_ylim([extents[2], extents[3]])

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

figname = f'{fig_folder}DAS_measurements_map_thickness_colorbar'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Create a map with superposed images 

# path to drone images and records 
date = '0210'
drone_ID = 'bernache'
exp_ID = '09-ortho_004'

base = 'U:/Data/'

records = field_drone.get_records('0210',jpg = False)
 # = f'{base}{date}/Drones/{drone_ID}/{exp_ID}/'
 
rec = records['drones'][drone_ID]['Drones'][8]
print(rec.keys())

filelist = glob.glob(f'U:/PIV_images/{date}/Drones/{drone_ID}/{exp_ID}/*.tiff')

#%% Load flightrecord

file2load = f'{base}{date}/Drones/{drone_ID}/flightrecords/Flightrecord_dict.pkl'
with open(file2load,'rb') as pf:
    records_csv = pickle.load(pf)

# cut flight records using .srt file t_start and t_end 
records_csv = field_drone.cut_flightrecord(rec, records_csv)

#%% Set drone parameters 

param = {}
param['focal'] = 2700
param['alpha'] = np.pi/2

#%% Convert csv time in time since epoch 

# conversion for csvflightrecord
records_csv['datetime'] = field_drone.get_datetime_from_csvflightrecord(records_csv)
records_csv['t_epoch'] = [d.timestamp() for d in records_csv['datetime']]

# conversion of SRT records
rec['datetime'] = field_drone.get_datetime_from_srtrecord(rec)
rec['t_epoch'] = [d.timestamp() for d in rec['datetime']]

#%% interpolate GPS coordinates and azimuth along time 
I = {}
I['latitude'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.latitude'])
I['longitude'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.longitude'])
I['azimuth'] = scipy.interpolate.interp1d(records_csv['t_epoch'], records_csv['OSD.yaw [360]'])

#%% Plot several frames 

alpha_0 = param['alpha']
focale = param['focal']

fig, ax = plt.subplots(figsize = (12,6))

N = 3500 #3500
step = 200 #100
indices_list = np.arange(0,N,step = step)
extents = []
for idx in indices_list:
    img_name = filelist[idx]
    img = cv.imread(img_name)
    img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
    gray_img = cv.cvtColor(img, cv.COLOR_RGB2GRAY)

    mean_gray = np.mean(gray_img)
    norm_img = gray_img/mean_gray
    
    H = dp.get_height_from_record(rec,indice = idx//100)
    
    Lat_D = I['latitude'](rec['t_epoch'][idx//100])
    Long_D = I['longitude'](rec['t_epoch'][idx//100])
    azimuth = I['azimuth'](rec['t_epoch'][idx//100])
    
    print(idx,H,Lat_D,Long_D,azimuth)
    
    # Lat_D = rec['latitude'][idx//100]
    # Long_D = rec['longitude'][idx//100]
    # print(H,Lat_D,Long_D)
    
    GPS_D = (Lat_D,Long_D,azimuth)
    Lat,Long = dp.georeference_from_param(img,H,alpha_0,focale,GPS_D)
    
    extent = [Long.min(),Long.max(),Lat.min(),Lat.max()]
    extents.append(extent)

    # ax.pcolormesh(Long,Lat,img[:,:,0], cmap = 'gray')
    # ax.imshow(np.transpose(img,(1,0,2)),extent = extent)
    ax.pcolormesh(Long,Lat,norm_img,cmap = 'gray',vmin = 0,vmax = 1.2,rasterized = True)


# plot DAS positions
# create a norm and a line collection 
points = np.array([DAS_water_height['long'],DAS_water_height['lat']]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
distance = DAS_water_height['s']
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 2)
lc.set_array(distance) # values used for colormap
ax.add_collection(lc)
# ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k')


# superpose deployments points 
ms = 70

# create norm for ice thickness
norm_acq = colors.Normalize(vmin = 1,vmax = 4)
full_cmap = mpl.colormaps['Oranges'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,1,256)))

color_date = {'0210':new_blues(0.2),'0211':new_blues(0.6), '0212': new_blues(0.9)}

# plot geophone lines and tomo
for geo_key in new_matrix.keys():
    for acqu_key in new_matrix[geo_key].keys():
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
        
        current_color = cmap(norm_acq(acqu_key))
        ax.scatter(GPS_logs['longitude'],GPS_logs['latitude'],marker = '^',color = current_color,
                   edgecolors = 'k',zorder = 2)

# plot ice thickness
# for key_date in records.keys():
#     for key in records[key_date].keys():
#         if 'H' in key:
#             h_value = int(key[2:])/100
#             current_color = cmap(norm_thick(h_value))
#             ax.scatter(records[key_date][key]['longitude'],records[key_date][key]['latitude'], marker = '^',
#                        color = current_color)

for key_date in gps_results.keys():
    current_color = color_date[key_date]
    ax.scatter(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],
               marker = 'p',color = current_color,edgecolors = 'k',zorder = 3)


extents = [-68.821900,-68.813627,48.346588,48.348010]
ax.set_xlim([extents[0], extents[1]])
ax.set_ylim([extents[2], extents[3]])

# ax.set_xlabel(r'Longitude (°)')
# ax.set_ylabel(r'Latitude (°)')

Lat0 = DAS_water_height['lat'][0]
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

figname = f'{fig_folder}DAS_measurements_image_instrument_superposition'
plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 300)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 300)


#%% Compute Flexural modulus for gps positions and geophones lines 

# geophone lines
model_key = '20250210'
sub_results[model_key]['D'] = sub_results[model_key]['Young modulus (Pa)']*sub_results[model_key]['thickness (m)']**3/12/(1-
                                                                        sub_results[model_key]['Poisson ratio']**2)

sub_results[model_key]['x'] = gps_pack.distance_GPS(sub_results[model_key]['latitude'], sub_results[model_key]['longitude'], 
                                         Lat0, Long0)

# gps thickness measurements 
E = np.mean(sub_results[model_key]['Young modulus (Pa)'])
nu = np.mean(sub_results[model_key]['Poisson ratio'])

for key_date in gps_results.keys():
    gps_results[key_date]['D'] = E*gps_results[key_date]['h']**3/12/(1-nu**2)
    gps_results[key_date]['x'] = gps_pack.distance_GPS(gps_results[key_date]['latitude'], gps_results[key_date]['longitude'], 
                                             Lat0, Long0)
    
#%%  Plot flexural modulus using geophones and thickness measurements - gps  

scale = 'log'
ms = 7

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()

color_date = {'0210':new_blues(0.2),'0211':new_blues(0.6), '0212': new_blues(0.9)}
key_swell = 'corrected'
for date in ['0211','0212'] :
    ax.errorbar(results_DAS['passive'][key_swell][date]['x'],results_DAS['passive'][key_swell][date]['D'],
                yerr = results_DAS['passive'][key_swell][date]['err_D'],
                fmt = '-o',label = date, color = color_date[date])

ax.plot(results_DAS['active']['0212']['x'],results_DAS['active']['0212']['D'],'o-',color = 'tab:red',label = 'Active 0212')

ax.plot(sub_results[model_key]['x'],sub_results[model_key]['D'],'s', color = color_date['0210'],ms = ms)
for key_date in gps_results.keys():
    ax.plot(gps_results[key_date]['x'],gps_results[key_date]['D'],'^', color = color_date[key_date],ms = ms)

ax.set_xlabel(r'$x \; \mathrm{(m)}$')
ax.set_ylabel(r'$D \; \mathrm{(J)}$')

ax.set_xlim([0,630])

if scale == 'log':
    ax.set_ylim([4e6,3e8])
    ax.set_yscale(scale)
    ax.legend(loc = 'lower right')
else :
    ax.set_ylim([0.5e6,2e8])
    ax.set_yscale(scale)
    ax.legend(loc = 'upper right')


figname = f'{fig_folder}D_vs_x_multi_instrument_superposition_{scale}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

