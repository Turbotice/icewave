# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 17:14:33 2025

@author: sebas
"""


import os
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
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

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.drone.drone_projection as dp 
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.rw_data as rw
import icewave.tools.weather as weather
import icewave.geophone.gps_coordinates as geophone_gps
import icewave.field.gps as field_gps

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Set fig_folder path 

fig_folder = 'U:/Data/0211/DAS/Figures_article/Situation_picture/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

#%% Load DAS GPS position 

file_water_height = 'U:/Data/0211/DAS/fiber_water_height_GPS_structure_0211.h5'
DAS_water_height = rw.load_dict_from_h5(file_water_height)

#%% Load Geophone GPS position 

dates = ['0210','0211','0212']

GPS = {'Geophones':{},'thickness':{}}
date = '0210'
year = '2025'

path2logfile = f'U:/Data/{date}/Geophones'
geophones_table_path = 'C:/Users/sebas/git/icewave/sebastien/geophones/geophones_table'
UTC_timearea = pytz.timezone('UTC')

# load geophone GPS coordinates 
new_matrix = geophone_gps.get_GPS_coordinates(path2logfile, geophones_table_path, date, year)
    
#%% Compute GPS_logs for each deployment 

GPS_geoph = np.zeros((len(new_matrix.keys()),len(new_matrix['G01'].keys()),2))

for i,geo_key in enumerate(new_matrix.keys()):
    for j,acqu_key in enumerate(new_matrix[geo_key].keys()):
        GPS_logs = new_matrix[geo_key][acqu_key]
        
        if  GPS_logs['num_GPS_logs'] > 1:
            GPS_logs = geophone_gps.compute_avg_logs(GPS_logs)
            
            GPS_geoph[i,j,0] = GPS_logs['longitude']
            GPS_geoph[i,j,1] = GPS_logs['latitude']
        else:
            
            GPS_geoph[i,j,0] = GPS_logs['longitude'][0]
            GPS_geoph[i,j,1] = GPS_logs['latitude'][0]
                            
        print(GPS_logs)


        
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
        
#%% Maddie SMith GPS waypoints
MS_gps = {}
MS_gps['longitude'] = [-68.814822, -68.819197]
MS_gps['latitude'] = [48.347624, 48.347229]
MS_gps['h'] = [0.31,0.56]

    
#%% Load situation picture 

main_path = 'U/Data/'
date = '0205'
drone_ID = 'mesange'
exp_ID = '18-doc_010'

path2img = f'U:/Data/{date}/Drones/{drone_ID}/{exp_ID}/'
filelist = glob.glob(f'{path2img}*.jpg')
file2read = filelist[0]
print(file2read)

# fig_folder = f'U:/Data/{date}/Drones/{drone_ID}/Figures/{exp_ID}/'
# if not os.path.isdir(fig_folder):
#     os.mkdir(fig_folder)
    
img = cv.imread(file2read)
img = cv.cvtColor(img,cv.COLOR_BGR2RGB)

# img_name = filename.split('.')[0]

#%% Plot raw picture 

fig, ax = plt.subplots(figsize = (12,9))
ax.imshow(img)
ax.set_xlabel(r'$x_p$')
ax.set_ylabel(r'$y_p$')

# select pixels coord of beginning and end of fiber 
fiber_tips = np.array([[2395,2969],[2593,566]])
x_fiber = np.linspace(fiber_tips[0,0],fiber_tips[1,0],100)
y_fiber = np.linspace(fiber_tips[0,1],fiber_tips[1,1],100)

ax.plot(x_fiber,y_fiber,'-')

#%% Plot line with a given color


# create a parameter s along the line
s  = np.linspace(0,1,600)
x_fiber = (1-s)*fiber_tips[0,0] + s*fiber_tips[1,0]
y_fiber = (1-s)*fiber_tips[0,1] + s*fiber_tips[1,1]

# create line segments (pais of consecutive points)
points = np.array([x_fiber,y_fiber]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
# distance = np.sqrt((x_fiber - fiber_tips[0,0])**2 + (y_fiber - fiber_tips[0,1])**2)

# apparent distance
length_fiber = 600
app_distance = s*length_fiber

# create a norm and a line collection 
norm = colors.Normalize(vmin = app_distance.min(), vmax = app_distance.max())
lc = LineCollection(segments,cmap = 'viridis', norm = norm, linewidths = 4)
lc.set_array(app_distance) # values used for colormap

set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
ax.imshow(img)
ax.add_collection(lc)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$s \; \mathrm{(m)}$')

# hide labels
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# figname = f'{fig_folder}Situation_picture_0205_18_doc_010'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Set drone parameters 

# Set drone parameters
param_dict = {}
param_dict['H'] = 110.8 #110.8
param_dict['alpha_0'] = 22.4
param_dict['focale'] = 4100 # around 4100 (+- 100) for format 5280 x 2970 # initial guess 4100

# -68.8131039300668
# -68.81312780555555 #- from image data

param_dict['latitude'] =  48.347697722222215 #- from image data #48.3476967601294 - Column BN True 
param_dict['longitude'] =  -68.81312780555555 #- from image data #-68.8131401977678 - Column BN True 
param_dict['azimuth'] = 263


#%% Georectify picture 

[ny,nx,nc] = np.shape(img) 

x_edges = np.arange(0,nx + 1)
y_edges = np.arange(0,ny + 1)

x0 = (nx + 1) / 2
y0 = (ny + 1) / 2

Yedges,Xedges = np.meshgrid(y_edges,x_edges,indexing = 'ij')

# compute real coordinates for each pixels of the image 
Xreal,Yreal = dp.projection_real_space(Xedges,Yedges,x0,y0,param_dict['H'],param_dict['alpha_0']*np.pi/180,
                                       param_dict['focale'])
# shift image position
Xshift = 0
Yshift = -17

Xreal = Xreal + Xshift
Yreal = Yreal + Yshift

# compute object position
POS_realx,POS_realy = dp.projection_real_space(x_fiber,y_fiber,x0,y0,
                                               param_dict['H'],param_dict['alpha_0']*np.pi/180,param_dict['focale'])

#%% Plot georectified image 

fig, ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Xreal,Yreal,img[:,:,0],shading = 'auto', cmap = 'gray')
# ax.plot(POS_realx,POS_realy,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

c.set_rasterized(True)
# figname = fig_folder + 'Camera_coord_situation_picture_0205_18_doc_010' 
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute image GPS position 

# horizontal distance between center of metric coordinates and drone position 
dist2drone = param_dict['H']/np.tan(param_dict['alpha_0']*np.pi/180)
Lat0,Long0 = dp.LatLong_coords_from_referencepoint(param_dict['latitude'],param_dict['longitude'],
                                                param_dict['azimuth'],dist2drone)

rho,theta = dp.cart2pol(Xreal,Yreal)
theta = theta*180/np.pi # angle with respect to the drone orientation 
local_azimuth = param_dict['azimuth'] + 90 - theta 
# GPS coordinates of all pixels
Lat,Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                local_azimuth,rho)

#GPS coordinates of detected objects 
POS_rho,POS_theta = dp.cart2pol(POS_realx,POS_realy)
POS_theta = POS_theta*180/np.pi
POS_azimuth = param_dict['azimuth'] + 90 - POS_theta
POS_Lat,POS_Long = dp.LatLong_coords_from_referencepoint(Lat0,Long0,
                                                POS_azimuth,POS_rho)

#%% Plot figure using GPS coordinates of fiber 

set_graphs.set_matplotlib_param('single')
fig,ax = plt.subplots(figsize = (12,9))
c = ax.pcolormesh(Long,Lat,img[:,:,0],shading = 'auto',cmap = 'gray')
c.set_rasterized(True)
# ax.plot(POS_Long,POS_Lat,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 10)
ax.set_xlabel(r'Longitude $(^\circ)$',labelpad = 15)
ax.set_ylabel(r'Latitude $(^\circ)$',labelpad = 15)
ax.set_aspect(1/np.cos(Lat0*np.pi/180)) # scaling y/x

c.set_rasterized(True)

# plot DAS GPS position 
ax.plot(DAS_water_height['long'],DAS_water_height['lat'],'k')

# plot Geophones GPS coord
norm_acq = colors.Normalize(vmin = 0,vmax = 3)
full_cmap = mpl.colormaps['Oranges'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,1,256)))

for j in range(4):
    current_color = cmap(norm_acq(j))
    ax.plot(GPS_geoph[:,j,0],GPS_geoph[:,j,1],'^',color = current_color,mec = 'k',zorder = 3)
    
# plot Drillings GPS coords
for key_date in gps_results.keys():
    ax.plot(gps_results[key_date]['longitude'],gps_results[key_date]['latitude'],'p',
            color = 'tab:blue',mec = 'k',zorder = 2)
    
# ax.set_xlim([-225,380])
# ax.set_ylim([-120,120])

figname = fig_folder + 'Situation_picture_0205_18_doc_010_DAS_geoph_drillings_GPS_coords' 
# plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
# plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
# plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)


#%% Position fiber on image in (X,Y) coordinates 

# convert DAS GPS coordinates into (X,Y) coordinates 
X_DAS,Y_DAS = dp.GPS2XY(DAS_water_height['lat'], DAS_water_height['long'], Lat0, Long0, param_dict['azimuth'])

# convert Geophones GPS coord into (X,Y) coordinates
X_geoph, Y_geoph = dp.GPS2XY(GPS_geoph[:,:,1],GPS_geoph[:,:,0],Lat0,Long0,param_dict['azimuth'])

# convert thickness GPS into (X,Y) coordinates
X_gps = {}
Y_gps = {}
for key_date in gps_results.keys():
    X_gps[key_date],Y_gps[key_date] = dp.GPS2XY(gps_results[key_date]['latitude'],gps_results[key_date]['longitude'],
                                                Lat0, Long0, param_dict['azimuth'])

X_MS, Y_MS = dp.GPS2XY(MS_gps['latitude'],MS_gps['longitude'],
                                            Lat0, Long0, param_dict['azimuth'])

# interpolate new positions along the fiber 
s  = np.linspace(0,1,600)
X_interp = (1-s)*X_DAS[0] + s*X_DAS[-1] - X_DAS[0]
Y_interp = (1-s)*Y_DAS[0] + s*Y_DAS[-1] - Y_DAS[0]

# create line segments (pais of consecutive points)
points = np.array([Y_interp,X_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp)**2 + (Y_interp)**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param('single')
markersize = 6
fig, ax = plt.subplots(figsize = (12,6))
c = ax.pcolormesh(Yreal - Y_DAS[0],Xreal - X_DAS[0],img[:,:,0],shading = 'auto', cmap = 'gray')
c.set_rasterized(True)
# ax.plot(X_DAS,Y_DAS,'o',markerfacecolor = 'r',markeredgecolor = 'k',markersize = 5)
ax.set_xlabel(r'$X \; \mathrm{(m)}$',labelpad = 5)
ax.set_ylabel(r'$Y \; \mathrm{(m)}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

# Plot DAS in color
# ax.add_collection(lc)

# Plot DAS without colors
ax.plot(Y_interp ,X_interp,'k',lw = 2.5,label = 'DAS')

# plot geophones
norm_acq = colors.Normalize(vmin = 0,vmax = 3)
full_cmap = mpl.colormaps['Greens'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,0.8,256)))

for j in range(4):
    current_color = cmap(norm_acq(j))
    if j == 0:
        ax.plot(Y_geoph[:,j]-Y_DAS[0],X_geoph[:,j]-X_DAS[0],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3,label = 'Geophone')
    else :
        ax.plot(Y_geoph[:,j]-Y_DAS[0],X_geoph[:,j]-X_DAS[0],'^',color = current_color,mec = 'k',
            ms = 6,zorder = 3)
        

# plot drillings
norm_h = colors.Normalize(vmin = 0.3,vmax = 0.9)# create colormap for drilling thickness
cmap = mpl.colormaps['coolwarm'].resampled(256)

for i,key_date in enumerate(gps_results.keys()):
    for j in range(len(gps_results[key_date]['h'])):
        current_color = cmap(norm_h(gps_results[key_date]['h'][j]))
        if j == 0 and i == 0:
            scatter = ax.plot(Y_gps[key_date][j]-Y_DAS[0],X_gps[key_date][j]-X_DAS[0],'p',color = current_color,mec = 'k',
                    ms = 8,zorder = 2,label = 'Drilling')
        else:
            scatter = ax.plot(Y_gps[key_date][j]-Y_DAS[0],X_gps[key_date][j]-X_DAS[0],'p',color = current_color,mec = 'k',
                    ms = 8,zorder = 2)

for i in range(len(MS_gps['h'])):
    current_color = cmap(norm_h(MS_gps['h'][i]))
    ax.plot(Y_MS[i] - Y_DAS[0],X_MS[i] - X_DAS[0],'p',color = current_color,mec = 'k',
            ms = 8,zorder = 2)
    
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

sm = cm.ScalarMappable(cmap = cmap, norm = norm_h)
sm.set_array([])  # Only needed for the colorbar
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label(r'$h \; \mathrm{(m)}$')

# cbar = plt.colorbar(scatter,cax = cax)

ax.set_xlim(np.array([-225,390]) - Y_DAS[0]) # [-225,380]
ax.set_ylim(np.array([-120,120]) - X_DAS[0]) # [-120,120]

ax.legend(ncols = 3, bbox_to_anchor = (0.11, 1),
              loc='lower left')

figname = fig_folder + 'Camera_coord_YX_situation_picture_0205_18_doc_010__with_MS_core_colorbar_h' 
plt.savefig(figname + '.pdf', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)

#%% Plot initial figure with fiber position in pixel coordinates system 

# DAS pixel coordinates
xpix_interp,ypix_interp = dp.projection_pixel_space(X_interp, Y_interp, x0, y0, param_dict['H'], 
                                               param_dict['alpha_0']*np.pi/180, param_dict['focale'])

# geophones pixel coordinates 
xpix_geoph,ypix_geoph = dp.projection_pixel_space(X_geoph, Y_geoph, x0, y0, param_dict['H'], 
                                               param_dict['alpha_0']*np.pi/180, param_dict['focale'])

# gps results pixel coordinates
xpix_gps = {}
ypix_gps = {}

for key_date in X_gps.keys():
    xpix_gps[key_date],ypix_gps[key_date] = dp.projection_pixel_space(X_gps[key_date], Y_gps[key_date], 
                                                                      x0, y0, param_dict['H'], 
                                                   param_dict['alpha_0']*np.pi/180, param_dict['focale'])


# create line segments (pairs of consecutive points)
points = np.array([xpix_interp,ypix_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp - X_DAS[0])**2 + (Y_interp - Y_DAS[0])**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
c = ax.imshow(img)
ax.set_xlabel(r'$x_{pix}$',labelpad = 5)
ax.set_ylabel(r'$y_{pix}$',labelpad = 5)
ax.set_ylim([4210,0])
ax.set_aspect(1) # set aspect ratio to 1 

ax.add_collection(lc)

# plot geophones
norm_acq = colors.Normalize(vmin = 0,vmax = 3)
full_cmap = mpl.colormaps['Oranges'].resampled(256)
cmap = colors.ListedColormap(full_cmap(np.linspace(0.2,1,256)))

for j in range(4):
    current_color = cmap(norm_acq(j))
    ax.plot(xpix_geoph[:,j],ypix_geoph[:,j],'^',color = current_color,mec = 'k',zorder = 3)

# plot thickness measurements
for key_date in X_gps.keys():
    mask = ypix_gps[key_date] > 550
    ax.plot(xpix_gps[key_date][mask],ypix_gps[key_date][mask],'p',color = 'tab:blue',mec = 'k',zorder = 2)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$x \; \mathrm{(m)}$')

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

figname = f'{fig_folder}Pixel_coords_situation_picture_0205_18_doc_010_DAS_geoph_gps_GPS_no_axis_labels'
plt.savefig(figname + '.pdf', bbox_inches='tight', dpi = 600)
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)

#%% Create figure for presentations

xpix_interp,ypix_interp = dp.projection_pixel_space(X_interp, Y_interp, x0, y0, param_dict['H'], 
                                               param_dict['alpha_0']*np.pi/180, param_dict['focale'])


# create line segments (pairs of consecutive points)
points = np.array([xpix_interp,ypix_interp]).T.reshape(-1,1,2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# distance from fiber beginning
distance = np.sqrt((X_interp - X_DAS[0])**2 + (Y_interp - Y_DAS[0])**2)

# create a norm and a line collection 
norm = colors.Normalize(vmin = distance.min(), vmax = distance.max())
lc = LineCollection(segments,cmap = parula_map, norm = norm, linewidths = 4)
lc.set_array(distance) # values used for colormap


set_graphs.set_matplotlib_param(('single'))
fig, ax = plt.subplots()
c = ax.imshow(img)

# ax.set_xlabel(r'$x_{pix}$',labelpad = 5)
# ax.set_ylabel(r'$y_{pix}$',labelpad = 5)
ax.set_aspect(1) # set aspect ratio to 1 

ax.add_collection(lc)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(lc,cax = cax)
cbar.set_label(r'$x \; \mathrm{(m)}$')

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

figname = f'{fig_folder}Pixel_coords_situation_picture_0205_18_doc_010_DAS_GPS_powerpoint_version'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.svg', bbox_inches='tight',dpi = 600)
plt.savefig(figname + '.png', bbox_inches='tight',dpi = 600)