# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 09:35:18 2025

@author: sebas

Compute DAS GPS position and bathymetry using waypoints of active sources from Februar 11th 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob 
import pickle 

import os 
import scipy

import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.gps.gps_seb as gps_seb
import icewave.tools.weather as weather
import icewave.drone.drone_projection as dp

parula_map = matcmaps.parula()
plt.rcParams.update({
    "text.usetex": True}) # use latex

#%% Import GPS waypoints 

date = '0211'
path2DAS = f'U:/Data/{date}/DAS/'

path2data = f'U:/Data/{date}/GPS/'
gpx_filelist = glob.glob(path2data + '*.gpx')

gpx_file = gpx_filelist[1]
# get gpx
gpx = gps_seb.get_gpx(gpx_file)

# get waypoints position
wpts = {}
wpts['name'] = [waypoint.name for waypoint in gpx.waypoints]
Long,Lat = gps_seb.waypoints_loc(gpx)
wpts['lat'] = Lat
wpts['long'] = Long

# keep only waypoints of sources
name_pattern = 'S'
source_wpts = {'name':[],'lat':[],'long':[]}
# add waypoint of DAS interrogation unit
source_wpts['name'].append(wpts['name'][-1])
source_wpts['lat'].append(wpts['lat'][-1])
source_wpts['long'].append(wpts['long'][-1])

for i in range(len(wpts['name'])):
    if 'S' in wpts['name'][i]:
        source_wpts['name'].append(wpts['name'][i])
        source_wpts['lat'].append(wpts['lat'][i])
        source_wpts['long'].append(wpts['long'][i])



#%% Check relative distance between each source 

Lat0 = wpts['lat'][-1]
Long0 = wpts['long'][-1]
rho = []
for i in range(len(source_wpts['name'])):
    rho.append(gps_seb.distance_GPS(source_wpts['lat'][i], source_wpts['long'][i], Lat0, Long0))

rho = np.array(rho)

fig, ax = plt.subplots()
ax.plot(rho,'o')

#%% get GPS position of DAS 

# compute azimuth between first and last point of DAS waypoints 
Lat0 = source_wpts['lat'][0]
Long0 = source_wpts['long'][0]

Lat1 = source_wpts['lat'][-1]
Long1 = source_wpts['long'][-1]
psi = 360 + np.arctan2(np.cos(Lat0*np.pi/180)*(Long1 - Long0)*np.pi/180,(Lat1 - Lat0)*np.pi/180)*180/np.pi

# deduce GPS position of DAS 
fiber_length = 600
s = np.arange(0,fiber_length + 10,step = 10)

fiber = {}
fiber['s'] = s
fiber_lat,fiber_long = dp.LatLong_coords_from_referencepoint(Lat0, Long0, psi, s)

fiber['lat'] = fiber_lat
fiber['long'] = fiber_long

#%% Superpose DAS GPS position with bathymetry 
# Load array of bathymetry data 
file2load = 'U:/General/Bathymetrie/NONNA_CHS/Bathymetry/NONNA10_4830N06890W_array.pkl'
with open(file2load,'rb') as pf:
    data_bath = pickle.load(pf)

#%% Plot bathymetry and DAS GPS position 

set_graphs.set_matplotlib_param('single')
Haha_edges = [-68.8395,-68.8017,48.3355,48.3575]
#keep only points within this area 
mask = np.logical_and(data_bath[:,1] > Haha_edges[0],
                      data_bath[:,1] < Haha_edges[1])
mask = np.logical_and(mask,data_bath[:,0] > Haha_edges[2])
mask = np.logical_and(mask,data_bath[:,0] < Haha_edges[3])

Haha_coordinates = data_bath[mask,:]

fig, ax = plt.subplots()
scatter = ax.scatter(Haha_coordinates[:,1],Haha_coordinates[:,0],s = 10,c = Haha_coordinates[:,2],cmap = parula_map)
ax.plot(source_wpts['long'],source_wpts['lat'],'r^')
ax.plot(fiber['long'],fiber['lat'],'k-')

ax.set_xlabel(r'Longitude (°)')
ax.set_ylabel(r'Latitude (°)')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cbar = plt.colorbar(scatter,cax = cax)
cbar.set_label(r'$H \; \mathrm{(m)}$')

ax.set_xlim([Haha_edges[0], Haha_edges[1]])
ax.set_ylim([Haha_edges[2], Haha_edges[3]])

fig_folder = f'{path2DAS}Figures/DAS_water_height/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

figname = f'{fig_folder}DAS_GPS_bathymetry_{date}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Load bathymetry interpolator

path2bath_interp = 'U:/Data/Bathymetrie/linear_interpolator_bathymetry.pkl'
# load interpolator 
with open(path2bath_interp,'rb') as pf:
    interp_H = pickle.load(pf)

#%% Compute bathymetry for all waypoints 

H = np.zeros((len(fiber['s'])))
for i in range(len(fiber['s'])):
    H[i] = weather.get_bathymetry_GPS((fiber['lat'][i],fiber['long'][i]), interp_H)

fiber['H'] = H
fiber['units'] = {'s':'meter','lat':'deg','long':'deg','H':'meter'}


file2save = f'{path2DAS}fiber_GPS_bathymetry_{date}.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(fiber,pf)











