#%% import libraries
import numpy as np
import math
from datetime import datetime
from datetime import timedelta
import mplcursors
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy import read
import shutil
import os
import re
import pickle
from scipy.fftpack import fft, ifft
from scipy.linalg import svd
import warnings
from matplotlib.patches import PathPatch
import gpxpy
from IPython import get_ipython
import matplotlib
import scipy
import scipy.signal
# %% Inputs

#get_ipython().run_line_magic('matplotlib', 'qt') # cellule à commenter ou pas suivant l'éditeur utilisé
plt.close('all')

# ## Input

# donnees exemple :
year = '2024'
date = '0306' #date format, 'mmdd'
acq_nb = 4
acqu_numb = '000'+str(acq_nb)+'_new' #acquisition number # _new parce que le dossier 0004 ne veut pas s'ouvrir (contiennent les memes données) 

"""
year = '2025'
date = '0210' #date format, 'mmdd'
acq_nb = 1
acqu_numb = '000'+str(acq_nb) #acquisition number # _new parce que le dossier 0004 ne veut pas s'ouvrir (contiennent les memes données) 
"""
#direction = 1 # 1 ou 2 
#composante = 'N' #Z , E or N -> direction de la source
channel = 2  # 0 for E, 1 for N, 2 for Z. 
#flexure_wave = 1 # 1 to pick the dispersion curves of the flexure wave, 0 to pick those of the other 2 modes

#files need to be organised as: data/0210/Geophones/0001/minised files

#path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN' # arborescence sur ordi ludo
#geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

# windows:
#path2logfile = f'D:/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur ssd vasco
#geophones_table_path = 'D:/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur ssd vasco

#linux:
ordi = 'dell_vasco'


#donnees exemple:
if ordi=='adour':
    path2data = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN'
    path2logfile = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur storageshared
    geophones_table_path = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
    gpx_file_path = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/1000_Waypoints_2024-03-06.gpx'
elif ordi=='dell_vasco':
    path2data = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN'
    path2logfile = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur storageshared
    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
    gpx_file_path = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/1000_Waypoints_2024-03-06.gpx'

"""
path2data = f'S:/Data/{date}/'
path2logfile = f'S:/Data/{date}/Geophones/'
geophones_table_path = 
"""
# %% open the GPS files (containing the waypoints)

# Path to your .gpx file
# windows path : (sur disque ssd)
#gpx_file_path = 'D:/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/1000_Waypoints_2024-03-06.gpx'
# linux path :

# Open the .gpx file
with open(gpx_file_path, 'r') as gpx_file:
    gpx = gpxpy.parse(gpx_file)



# %% mettre toutes les données des gps, prises dans le fichier gpx, dans des listes
list_Waypoints_gpx = []
list_Latitudes_gpx = []
list_Longitudes_gpx = []
list_Elevation_gpx = []
list_Time_gpx = []
list_Symbol_gpx = []

for waypoint in gpx.waypoints:
    list_Waypoints_gpx.append(waypoint.name)
    list_Latitudes_gpx.append(waypoint.latitude)
    list_Longitudes_gpx.append(waypoint.longitude)
    list_Elevation_gpx.append(waypoint.elevation)
    list_Time_gpx.append(waypoint.time)
    list_Symbol_gpx.append(waypoint.symbol)


# %% Sélectionner les géophones, et la source qui nous intéressent, et extraire 
# l'instant des waypoints correspondants ainsi que les coordonnées spatiales (longitude,latitude)

# Define the path to the text file
#aparatus_numbers = [420,404,407,410,413] # numéros trouvés grace à la carte Tomo de Stéphane : file:///D:/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Carte_Tomo_2024_0306labeled.pdf
sources_numbers = [420,409,410] # numéros trouvés grace à la carte Tomo de Stéphane : file:///D:/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Carte_Tomo_2024_0306labeled.pdf
#aparatus_types = ['S','G','G','G','G']
list_idx_geophone = [4,7,10,13] # pareil, voir carte tomo de Stephane

#### pour voir toutes les sources utilisées pour la tomo:
list_idx_geophone = np.arange(1,17,1) # pareil, voir carte tomo de Stephane
sources_numbers = np.arange(401,423,1) # numéros trouvés grace à la carte Tomo de Stéphane : file:///D:/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Carte_Tomo_2024_0306labeled.pdf



#file_path = 'D:/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Map_Table.txt'
file_path = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Map_Table.txt'
file_path = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/2024_BICWIN/0306/GPS/Map_Table.txt'


# Read the file and store the data in a list of tuples
with open(file_path, 'r') as file:
    data = [line.strip().split() for line in file]

waypoints_numbers_sources = []
for i in range(len(data)):
    for j in range(len(sources_numbers)):
        if data[i][1]=='S_0'+str(sources_numbers[j]): # ici tous les appareils qui nous intéressent pour les waypoints sont les source (donc on prend que les 'S')
            print(data[i][1],int(data[i][0]))
            waypoints_numbers_sources.append(int(data[i][0]))
#chosen_numbers = [477,443,457,452,446] # trouvés à la main à la premiere ligne du fichier map_table
#%%
prefix = 'Sag240'
waypoints_names_sources = []
for i in range(len(waypoints_numbers_sources)):
    waypoints_names_sources.append(prefix+str(waypoints_numbers_sources[i]))

indices_sources = []
for i in range(len(list_Waypoints_gpx)):
    for name in waypoints_names_sources:
        if list_Waypoints_gpx[i]==name:
            indices_sources.append(i)
            print(list_Waypoints_gpx[i])
array_indices_sources = np.array(indices_sources)

array_Waypoints_gpx = np.array(list_Waypoints_gpx)
array_Latitudes_gpx = np.array(list_Latitudes_gpx)
array_Longitudes_gpx = np.array(list_Longitudes_gpx)
array_Elevation_gpx = np.array(list_Elevation_gpx)
array_Time_gpx = np.array(list_Time_gpx)
array_Symbol_gpx = np.array(list_Symbol_gpx)


# %% Plot longitude vs latitude
plt.figure()
plt.plot(list_Longitudes_gpx,list_Latitudes_gpx,'o')
plt.show()


# montrer seulement les appareils choisis
#colors_aparatus = convert_aparatus_to_color(aparatus_types)

plt.figure()
plt.plot(array_Longitudes_gpx[indices_sources],array_Latitudes_gpx[indices_sources],'bo',label='Sources')
#plt.plot(array_Longitudes[array_indices_sources][5],array_Latitudes[array_indices_sources][5],'ob',label='position source')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Coordinates of geophones and sources')
plt.legend()
plt.show()

print('time when the sources were disposed :',array_Time_gpx[indices_sources])

# %% dans cette partie on convertit pas vraiment en coordonnées (E,N) donc c'est pas super précis...
# suppose (to make simple) that all points are on a line
# let us compute distances with respect to the source

longitude_S = array_Longitudes_gpx[indices_sources][1]
latitude_S = array_Latitudes_gpx[indices_sources][1] # indice 0 car il peut y avoir plusieurs sources auquelles on s'intéresse (donc tableau)

#longitude_G = array_Longitudes[array_indices_sources][:5]
#latitude_G = array_Latitudes[array_indices_sources][:5]

# -----------------------------------------------------------------------------------------------------------------------------
#%% code qui ouvre les coordonnées GPS des géophones rangées dans un fichier par le programme GPS_coordinates_from_LOG_files.py
if ordi=='adour':
    path2coordfile = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/'+year+'_BICWIN/'+date+'/Geophones/GPS_coordinates_acq'+str(acq_nb)+'.txt'
    path2geophonetable = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'
elif ordi=='dell_vasco':
    path2coordfile = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/'+year+'_BICWIN/'+date+'/Geophones/GPS_coordinates_acq'+str(acq_nb)+'.txt'
    path2geophonetable = 'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'



def find_LongLat_from_idx_geophone(idx_geophone):
    geophones_table = np.loadtxt(path2geophonetable,skiprows=1)
    geo_SN = int(geophones_table[:,1][np.where(idx_geophone==geophones_table)[0][0]])
    geo_SN_str = str(geo_SN)
    geo_SN_last4str = geo_SN_str[-4:]
    geo_coord_table = np.loadtxt(path2coordfile,skiprows=1,dtype=str)
    indrow = np.where(geo_coord_table=='IGU_'+geo_SN_last4str)[0][0]
    latitude = float(geo_coord_table[indrow,1])
    longitude =  float(geo_coord_table[indrow,2])
    return longitude,latitude

def read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot):
    selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]

    current_stream = selected_streams[chan]
    print(current_stream[0])
    
    #ax[k].plot(np.array(datetime_values)[datetime_indices_to_plot],current_stream[0].data[datetime_indices_to_plot] / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) ),label = channel_dic[chan])        
    # ici il ne faut pas normaliser le signal (plus haut c'est fait seulement pour plotter tous les signaux)
    signal = current_stream[0].data[datetime_indices_to_plot]# / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) )
    return signal
def LongLat2deltaEdeltaN(coord1,coord2,latitude_approx,R=6371*1e3):

    coord1_rad = np.radians(np.array(coord1))
    coord2_rad = np.radians(np.array(coord2))
    latitude_approx_rad = np.radians(latitude_approx)
    
    deltaE = R * np.sin(np.pi/2-latitude_approx_rad) * (coord2_rad[0]-coord1_rad[0])
    deltaN = R * (coord2_rad[1]-coord1_rad[1])
    
    return deltaE,deltaN    

list_deltaE = [] 
list_deltaN = [] 
list_Long_geo = []
list_Lat_geo = []

for idx_geophone in list_idx_geophone: 
    Long,Lat = find_LongLat_from_idx_geophone(idx_geophone) 
    list_Long_geo.append(Long) 
    list_Lat_geo.append(Lat) 

for i in range(len(list_Long_geo)): 
     
    deltaE,deltaN = LongLat2deltaEdeltaN((longitude_S,latitude_S),(list_Long_geo[i],list_Lat_geo[i]),latitude_S) 
    list_deltaE.append(deltaE) 
    list_deltaN.append(deltaN) 

plt.figure()
plt.plot(list_Long_geo,list_Lat_geo,'s',color='tab:orange',label='geophones (from geophone log files)')
plt.plot(array_Longitudes_gpx[indices_sources],array_Latitudes_gpx[indices_sources],'ob',label='source (from gpx waypoints)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("GPS coord of sources and geophones")
plt.legend()
plt.show()                                  


plt.figure()
plt.plot(list_deltaE,list_deltaN,'s',color='tab:orange',label='geophones (from geophone log files)')
plt.plot(0,0,'ob',label='source (from gpx waypoints)')
plt.xlabel('"E" position (m)')
plt.ylabel('"N" position (m)')
plt.title("Positions converted in cartesian coord")
plt.legend()
plt.show()                                  

# %% code qui donne à la source la position (x,y) = (0,0) (en mètres), et qui déduit les coordonnées des géophones en cartésiennes
'''


longitude_S = array_Longitudes_gpx[indices_sources][0]
latitude_S = array_Latitudes_gpx[indices_sources][0]

E_S = 0
N_S = 0

E_G = np.zeros(len(longitude_G))
N_G = np.zeros(len(latitude_G))

for i in range(len(longitude_G)):
    E_G[i],N_G[i] = LongLat2deltaEdeltaN((longitude_S,latitude_S),(longitude_G[i],latitude_G[i]),latitude_S)
'''

# %% Fonctions (mais pas celles qui calculent les relations f vs k)


#----------------------------------------------------------------------------------------------
def read_data(path2data):
    miniseed_files = []

    print(path2data)
    # Iterate over files in the specified directory
    for filename in os.listdir(path2data):
        if filename.lower().endswith(".miniseed"):
#            file_path = os.path.join(path2data, filename)
            file_path = path2data + '/' + filename
            miniseed_files.append(file_path)

    # Read MiniSEED files
    streams = []
    for file_path in miniseed_files:
        stream = read(file_path)
        streams.append(stream)

    return streams, miniseed_files
#----------------------------------------------------------------------------------------------
def convert_to_utc_times(trace, time_vector):
    start_time_utc = UTCDateTime(trace.stats.starttime)
    return [start_time_utc + t for t in time_vector]
#----------------------------------------------------------------------------------------------
def rename_traces(stream, geophones_dict):
    for trace in stream:
        # Extract the last 5 digits after "SS."
        last_five_digits = trace.stats.station[-5:]

        # Check if the last five digits exist in the geophones dictionary
        if last_five_digits in geophones_dict:
            # Get the corresponding 2-digit number from the geophones dictionary
            new_station_code = geophones_dict[last_five_digits]

            # Replace the last 5 digits in the existing station code with the new ones
            trace.stats.station = f"{trace.stats.station[:-5]}{new_station_code}.{trace.stats.channel}"
        else:
            print(
                f"Warning: No entry found for {last_five_digits} in the geophones table.")

    # Sort traces based on the new station codes
    sorted_stream = sorted(stream, key=lambda trace: trace.stats.station)

    return sorted_stream
#----------------------------------------------------------------------------------------------

#%% Load geophone data
#table into a dictionary
geophones_dict = {}
with open(geophones_table_path, 'r') as table_file:
    for line in table_file:
        if line.startswith("Num"):
            continue  # Skip header line
        num, geo_sn = line.strip().split()
        # Assuming the last five digits are relevant for comparison
        last_five_digits = geo_sn[-5:]
        geophones_dict[last_five_digits] = num

# Read MiniSEED file directly
seismic_data_streams, miniseed_files = read_data (path2data +'/'+ date + '/Geophones/' + acqu_numb)

# Iterate over streams and rename traces
for i, stream in enumerate(seismic_data_streams):
    seismic_data_streams[i] = rename_traces(stream, geophones_dict)
    

def sort_key(trace):
    return trace[0].stats.station
# Sort the seismic_data_streams based on the custom sorting function
seismic_data_streams = sorted(seismic_data_streams, key=sort_key)

first_stream = seismic_data_streams[1]  #first stream is a obspy.core.trace object
first_trace = first_stream[0] #first stream is SS.01.SSN..SSN

# Find the largest start time (tstart) and the smallest end time (tend) among all streams
tstart = max([stream[0].stats.starttime for stream in seismic_data_streams])
tend = min([stream[-1].stats.endtime for stream in seismic_data_streams])
# Truncate all streams between tstart and tend
for i, stream in enumerate(seismic_data_streams):
    for trace in stream:
        if trace.stats.starttime < tstart:
            trace.trim(tstart, tend, pad=True, fill_value=0)
        elif trace.stats.endtime > tend:
            trace.trim(tstart, tend, pad=True, fill_value=0)


fs = first_trace.stats.sampling_rate # Needed, fn_svd
selected_indices = range(channel, len(seismic_data_streams), 3) # Needed for seismic_matrix


# If "Unknown format for file 0005\._01.0005.2024.02.11.20.06.55.000.E.miniseed" appears
# Create a NEW 000n file only with the miniseed files. 
# assign a string to the channel values 
signal_length = 1
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}
ch = channel_dic[channel]

####################### 
first_stream = seismic_data_streams[0]
first_trace = first_stream[0]
time_array = first_trace.times() # experimental time array 
time_vector = convert_to_utc_times(first_trace, first_trace.times())
datetime_values = [datetime.utcfromtimestamp(t.timestamp) for t in time_vector] # create datetime UTC array 
data_vector = first_trace.data
start_time_utc = UTCDateTime(first_trace.stats.starttime) 
fs = first_trace.stats.sampling_rate # acquisition frequency (Hz) 


#%% ---------------------- Plot all channels for the chosen geophones ---------------------------(vasco's modifications)
signal_length = 1 # qu'est ce que ce parametre?
channel_dic = {
    1: "N",
    2: "Z",
    0: "E",}
ch = channel_dic[channel]

#list_idx_geophone = [4]


#datetime_indices_to_plot = np.where((np.array(datetime_values)>datetime(2024,3,6,19,17,0))&(np.array(datetime_values)<datetime(2024,3,6,19,18,15)))
array_Time_gpx_without_tz = np.zeros_like(array_Time_gpx)
for i in range(len(array_Time_gpx)):
    array_Time_gpx_without_tz[i] = array_Time_gpx[i].replace(tzinfo=None)
time_to_add_after = timedelta(seconds=20)
datetime_indices_to_plot = np.where(((np.array(datetime_values)>array_Time_gpx_without_tz[indices_sources][0])&(np.array(datetime_values)<array_Time_gpx_without_tz[indices_sources][1]+time_to_add_after)))
#datetime_indices_to_plot = np.arange(len(np.array(datetime_values)))

fig, ax = plt.subplots(3,figsize = (10,10))

for i in range(len(list_idx_geophone)):
    idx_geophone = list_idx_geophone[i] # index of the geophone to be plotted

    selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]

    for k, chan in enumerate([2,0,1]):
            current_stream = selected_streams[chan]
            print(current_stream[0])
            #ax[k].plot(first_trace.times(),current_stream[0].data / max(np.abs(current_stream[0].data) ),label = channel_dic[chan])
            ax[k].plot(np.array(datetime_values)[datetime_indices_to_plot],current_stream[0].data[datetime_indices_to_plot] / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) ),label = channel_dic[chan])        
            ax[k].legend()

            ax[k].set_ylabel(r'$V/V_{max}$')
            ax[k].set_ylim([-1.1, 1.1])

            # pour afficher le moment où la source a été posée
            ax[k].vlines(array_Time_gpx[indices_sources][0],-1,1,color='red',linestyle='--')

    ax[2].set_xlabel('$t \; \mathrm{(s)}$')

fig.tight_layout()

# --------------------------------------------------------------------------------------------------------------------------------
# %% code qui donne à la source la position (x,y) = (0,0) (en mètres), et qui déduit les coordonnées des géophones en cartésiennes

def LongLat2deltaEdeltaN(coord1,coord2,latitude_approx,R=6371e3):
    coord1_rad = np.radians(np.array(coord1))
    coord2_rad = np.radians(np.array(coord2))
    latitude_approx_rad = np.radians(latitude_approx)
    
    deltaE = R * np.sin(np.pi/2-latitude_approx_rad) * (coord2_rad[0]-coord1_rad[0])
    deltaN = R * (coord2_rad[1]-coord1_rad[1])
    
    return deltaE,deltaN

#longitude_S = array_Longitudes_gpx[array_indices_sources][0]
#latitude_S = array_Latitudes_gpx[array_indices_sources][0]

#longitude_G = array_Longitudes_gpx[array_indices_sources][:5]
#latitude_G = array_Latitudes_gpx[array_indices_sources][:5]

#E_S = 0
#N_S = 0

#E_G = np.zeros(len(longitude_G))
#N_G = np.zeros(len(latitude_G))

#for i in range(len(longitude_G)):
#    E_G[i],N_G[i] = LongLat2deltaEdeltaN((longitude_S,latitude_S),(longitude_G[i],latitude_G[i]),latitude_S)
# %% dans cette partie on convertit pas vraiment en coordonnées (E,N) donc c'est pas super précis...
# suppose (to make simple) that all points are on a line
# let us compute distances with respect to the source

#longitude_S = array_Longitudes_gpx[indices_sources][0]
#latitude_S = array_Latitudes_gpx[indices_sources][0] # indice 0 car il peut y avoir plusieurs sources auquelles on s'intéresse (donc tableau)

#longitude_G = array_Longitudes[indices_to_study][:5]
#latitude_G = array_Latitudes[indices_to_study][:5]

# -----------------------------------------------------------------------------------------------------------------------------
#%% code qui ouvre les coordonnées GPS des géophones rangées dans un fichier par le programme GPS_coordinates_from_LOG_files.py
path2coordfile = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/'+year+'_BICWIN/'+date+'/Geophones/GPS_coordinates_acq'+str(acq_nb)+'.txt'
path2geophonetable = '/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table'

def find_LongLat_from_idx_geophone(idx_geophone):
    geophones_table = np.loadtxt(path2geophonetable,skiprows=1)
    geo_SN = int(geophones_table[:,1][np.where(idx_geophone==geophones_table)[0][0]])
    geo_SN_str = str(geo_SN)
    geo_SN_last4str = geo_SN_str[-4:]
    geo_coord_table = np.loadtxt(path2coordfile,skiprows=1,dtype=str)
    indrow = np.where(geo_coord_table=='IGU_'+geo_SN_last4str)[0][0]
    latitude = float(geo_coord_table[indrow,1])
    longitude =  float(geo_coord_table[indrow,2])
    return longitude,latitude

def read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot):
    selected_streams = seismic_data_streams[(idx_geophone - 1)*3 : (idx_geophone - 1)*3 + 3]

    current_stream = selected_streams[chan]
    print(current_stream[0])
    
    #ax[k].plot(np.array(datetime_values)[datetime_indices_to_plot],current_stream[0].data[datetime_indices_to_plot] / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) ),label = channel_dic[chan])        
    # ici il ne faut pas normaliser le signal (plus haut c'est fait seulement pour plotter tous les signaux)
    signal = current_stream[0].data[datetime_indices_to_plot]# / max(np.abs(current_stream[0].data[datetime_indices_to_plot]) )
    return signal
    

list_deltaE = [] 
list_deltaN = [] 
list_Long_geo = []
list_Lat_geo = []

for idx_geophone in list_idx_geophone: 
    Long,Lat = find_LongLat_from_idx_geophone(idx_geophone) 
    list_Long_geo.append(Long) 
    list_Lat_geo.append(Lat) 

for i in range(len(list_Long_geo)): 
     
    deltaE,deltaN = LongLat2deltaEdeltaN((longitude_S,latitude_S),(list_Long_geo[i],list_Lat_geo[i]),latitude_S) 
    list_deltaE.append(deltaE) 
    list_deltaN.append(deltaN) 

plt.figure()
plt.plot(list_deltaE,list_deltaN,'s',color='tab:orange',label='geophones (from geophone log files)')
plt.plot(E_G,N_G,'^g',label='geophone (from gpx waypoints)')
plt.plot(0,0,'ob',label='source (from gpx waypoints)')
plt.xlabel('"E" position (m)')
plt.ylabel('"N" position (m)')
plt.title("Positions converted in cartesian coord")
plt.legend()
plt.show()                                  

#%% afficher les signaux longitudinaux et transverses
# on commence par donner l'angle avec lequel la ligne de géophones est inclinée (dans le repère (E,N) choisi)
def compute_mean_angle(list_deltaE,list_deltaN):
    '''
    Fonction qui calcule l'angle moyen que l'on considère pour la tomographie "sur une ligne".
    Attention : il faut bien que tous les géophones ainsi que la source soient alignés !!!
    '''
    angles = np.zeros(len(list_deltaE))
    for i in range(len(list_deltaE)):
        angles[i] = np.arctan2(list_deltaN[i],list_deltaE[i])
    mean_angle = np.mean(angles) # on peut le calculer plus précisément...
    return mean_angle

mean_angle = compute_mean_angle(list_deltaE,list_deltaN)

# axial mode est appelé QS0
# la formule à utiliser sera : signal_QS0 = np.cos(mean_angle) * signal_E + np.sin(mean_angle) * signal_N

def extract_longi_transv_signal(idx_geophone,direction):
    '''
    prend en argument idx_geophone (le numéro du geophone considéré),
    et direction (qui peut être 'L' pour longi et 'T' pour transverse)
    '''
    chan = 0 # for E
    signal_E = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
    chan = 1 # for N
    signal_N = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
    if direction == 'L':
        signal_longi = np.cos(mean_angle) * signal_E + np.sin(mean_angle) * signal_N
        return signal_longi
    elif direction == 'T':
        signal_transv = np.sin(mean_angle) * signal_E - np.cos(mean_angle) * signal_N
        return signal_transv
# test
signal_longi4 = extract_longi_transv_signal(4,'L')
signal_transv4 = extract_longi_transv_signal(4,'T')
signal_longi7 = extract_longi_transv_signal(7,'L')
signal_transv7 = extract_longi_transv_signal(7,'T')

plt.figure()
chan = 0 # for E
signal_E = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
chan = 1 # for N
signal_N = read_data_choosing_geophone_channel(idx_geophone,chan,datetime_indices_to_plot)
plt.plot(signal_N,label='N')
plt.plot(signal_E,label='E')
plt.plot(signal_longi4,label='longi')
plt.plot(signal_transv4,label='transv')
plt.title('comparaison signalE, signalN, signal longi')
plt.legend()
plt.show()
plt.plot(signal_E[26000:29000],label='E')
plt.legend()
plt.show()
# shear horizontal mode est appelé QS0
# la formule à utiliser sera : signal_transv = np.sin(mean_angle) * signal_E - np.cos(mean_angle) * signal_N
plt.figure()
plt.plot(signal_longi4,label='longi4')
plt.plot(signal_longi7,label='longi7')
plt.legend()
plt.show()
plt.figure()
plt.plot(signal_longi4,label='transv4')
plt.plot(signal_longi7,label='transv7')
plt.legend()
plt.show()

#%% correler les signaux des geophones successifs pour obtenir le déphasage temporel
from scipy.signal import find_peaks
# pour commencer il faut faire un code qui détecte les moments où on tape, pour pouvoir ensuite zoomer si on veut faire des correlations
#typical_signal = signal_longi4[13090:13500] # signal typique dans la premiere tomo
typical_signal = signal_E[27000:28000]
typical_signal = typical_signal/np.max(np.abs(typical_signal)) # typical signal is normalized
typical_signal_withzeros = np.zeros(len(signal_longi4))
typical_signal_withzeros[:len(typical_signal)] = typical_signal

def find_9shocks(signal,typical_signal):# adjust height for find_peaks
    """
    function to find the 3*3 = 9 shocks in the desired signal
    inputs : 
    - signal (not normalized)
    - typical signal with which we want to correlate the signal to find the times of the 9 shocks (not normalized)
    Remark : returns indices (not times in UTC) that index the array
    """
    typical_signal = typical_signal/np.max(np.abs(typical_signal)) # typical signal is normalized
    typical_signal_withzeros = np.zeros(len(signal))
    typical_signal_withzeros[:len(typical_signal)] = typical_signal

    ycorr = scipy.signal.correlate(signal,typical_signal_withzeros)
    lags = scipy.signal.correlation_lags(len(signal),len(typical_signal_withzeros))

    pp = find_peaks(ycorr)
    #ppp = find_peaks(ycorr[pp[0]],width=2,distance=100,threshold=1e-2,height=80)
    ppp = find_peaks(ycorr[pp[0]],width=2,distance=100,threshold=1e-2,height=20)
    peak_heights = ppp[1]['peak_heights']
    argsort = np.argsort(peak_heights)
    nine_highest_peaks_idx = argsort[:9]
    plt.figure()
    plt.plot(lags,ycorr)
    plt.show()
    plt.plot(lags[pp[0]],ycorr[pp[0]])
    plt.vlines(lags[pp[0]][ppp[0]][nine_highest_peaks_idx],ymin=0,ymax=np.max(ycorr),colors='red',linestyles='--')
    print(lags[pp[0]][ppp[0]][nine_highest_peaks_idx])
    plt.show()
    if len(peak_heights)<9:
        print('didnt find all the shocks')
    return lags[pp[0]][ppp[0]]



# tests : 

find_9shocks(signal_longi4,typical_signal)
find_9shocks(signal_longi7,typical_signal)
#for i in range(3):
#    plt.plot(array_c_longi[i,:],'ob')
#plt.ylim(1500,1900)
find_9shocks(signal_transv4,typical_signal)
find_9shocks(signal_transv7,typical_signal)

#signalZ_4 = read_data_choosing_geophone_channel(4,2,datetime_indices_to_plot) # 2 correspond au channel 'Z'
#signalZ_7 = read_data_choosing_geophone_channel(7,2,datetime_indices_to_plot) # 2 correspond au channel 'Z'

#find_9shocks(signalZ_4,typical_signal)

# maintenant on fait une fonction qui précise suivant la direction du coup sur la glace qui renvoie les instants correspondants
def find_directed_shocks(indices_9shocks,direction):
    indices_9shocks = np.sort(indices_9shocks)
    if direction=='Z':# vertical component
        return indices_9shocks[:3]
    elif direction=='R':# 'R' for RADIAL with respect to the circle made by the sources
        return indices_9shocks[3:6]
    elif direction=='C':# 'C' for Colinear with respect to the circle made by the sources
        return indices_9shocks[6:9]

# maintenant on coupe deux signaux correspondant à deux géophones successifs pour focus seulement au voisinage du coup
def cut_2signals(signal1,signal2,idx_shock,nb_points_after_shock=2000):
    signal1_cut = signal1[idx_shock-nb_points_after_shock:idx_shock+nb_points_after_shock]
    signal2_cut = signal2[idx_shock-nb_points_after_shock:idx_shock+nb_points_after_shock]
    return signal1_cut,signal2_cut

# test :
from haversine import haversine

indices_9shocks = find_9shocks(signal_longi4,typical_signal)  # on prend comme indices de chocs de préférence ceux qui correspondent au géophone le plus près de la source
indices_shocks_longi = find_directed_shocks(indices_9shocks,'R')
signal_longi4_cut,signal_longi7_cut = cut_2signals(signal_longi4,signal_longi7,indices_shocks_longi[0])


# fonction qui fait tout ça :np.sqrt((8e9)/(2*1e3*(1+0.33)))
def detect_time_shift(signal1_cut,signal2_cut,num_geo1,num_geo2,freq_acq=1000):
    cross_corr = scipy.signal.correlate(signal2_cut,signal1_cut)
    lags = scipy.signal.correlation_lags(len(signal2_cut),len(signal1_cut))
    plt.plot(lags,cross_corr)
    plt.show()
    ind1 = np.where(np.array(list_idx_geophone)==num_geo1)[0][0]
    ind2 = np.where(np.array(list_idx_geophone)==num_geo2)[0][0]
    lat1,lon1 = list_Lat_geo[ind1],list_Long_geo[ind1]
    lat2,lon2 = list_Lat_geo[ind2],list_Long_geo[ind2]
    dis_1_2 = haversine((lat1,lon1),(lat2,lon2))*1e3 # suppose que la source et les geophones sont alignés!!!

    idx_max = np.argmax(cross_corr)

    # detection subpixellaire du max :
    x1 = lags[idx_max-1]
    y1 = cross_corr[idx_max-1]
    x2 = lags[idx_max]
    y2 = cross_corr[idx_max]
    x3 = lags[idx_max+1]
    y3 = cross_corr[idx_max+1]

    A = np.array([[x1**2,x1,1],[x2**2,x2,1],[x3**2,x3,1]])
    B = np.array([y1,y2,y3])
    Coefficients = np.linalg.solve(A,B)

    x_max_subpix = -Coefficients[1]/(2*Coefficients[0])
    delta_t_subpix = x_max_subpix/freq_acq
    return delta_t_subpix,dis_1_2

delta_t_subpix_4_7,dis_4_7 = detect_time_shift(signal_longi4_cut,signal_longi7_cut,4,7)

print("vitesse de l'onde entre les deux géophones 4 et 7 (détection subpixellaire): ",dis_4_7/delta_t_subpix_4_7)

# %% faire une fonction qui prend comme argument 2 numéros des géophones (cf carte Stéphane),
# la direction choisie (longitudinale ou transverse, sachant que source et géophones placés sur une ligne),
#  le signal typique choisi pour faire la csignal1,signal2,idx_shock,nb_points_after_shock=2000orrelation, et renvoie le décalage spatial et temporel
#  entre les géophones pour les trois différents chocs qui ont lieu sur la glace à l'emplacement de la source

def extract_3temporal_shifts_between_2geophones(num_geo1,num_geo2,direction): # ,typical_signal
    """
    fonction qui prend comme argument 2 numéros des géophones (cf carte Stéphane),
    la direction choisie (longitudinale ou transverse, sachant que source et géophones placés sur une ligne),
    le signal typique choisi pour faire la correlation, et renvoie le décalage spatial et temporel
    entre les géophones pour les trois différents chocs qui ont lieu sur la glace à l'emplacement de la source
    les 3 signaux comparés (car 3 coups par direction) sont indicés a,b,c respectivement
    Return : array_delta_t_subpix (len=3),distance entre les 2 géophones, array vitesses mesurée pour chaque coup (len=3) en m/s
    """
    # etape 1 : extraire les signaux désirés
    signal1 = extract_longi_transv_signal(num_geo1,direction)
    signal2 = extract_longi_transv_signal(num_geo2,direction)
    # etape 2 : les couper au bon endroit en fonction la direction choisie
    # indices_9shocks doit etre prealablement définie
    if direction=='L':
        signal1_cut_a , signal2_cut_a = cut_2signals(signal1,signal2,indices_9shocks[3])
        signal1_cut_b , signal2_cut_b = cut_2signals(signal1,signal2,indices_9shocks[4])
        signal1_cut_c , signal2_cut_c = cut_2signals(signal1,signal2,indices_9shocks[5])
    elif direction=='T':
        signal1_cut_a , signal2_cut_a = cut_2signals(signal1,signal2,indices_9shocks[6])
        signal1_cut_b , signal2_cut_b = cut_2signals(signal1,signal2,indices_9shocks[7])
        signal1_cut_c , signal2_cut_c = cut_2signals(signal1,signal2,indices_9shocks[8])
    print(len(signal1_cut_a),len(signal1_cut_b),len(signal1_cut_c))
    # etape 3 : calculer les shifts pour les 3 coups successifs
    # remarque : dans les 3 cas la distance reste la meme, c'est pourquoi on l'overwrite...
    delta_t_subpix_a,dis_1_2 = detect_time_shift(signal1_cut_a,signal2_cut_a,num_geo1,num_geo2)
    delta_t_subpix_b,dis_1_2 = detect_time_shift(signal1_cut_b,signal2_cut_b,num_geo1,num_geo2)
    delta_t_subpix_c,dis_1_2 = detect_time_shift(signal1_cut_c,signal2_cut_c,num_geo1,num_geo2)
    array_delta_t_subpix = np.array([delta_t_subpix_a,delta_t_subpix_b,delta_t_subpix_c])
    return (array_delta_t_subpix,dis_1_2,dis_1_2/np.array([delta_t_subpix_a,delta_t_subpix_b,delta_t_subpix_c]))

# tests : 
array_c_longi = np.zeros((3,3))
array_c_longi[:,0] = extract_3temporal_shifts_between_2geophones(4,7,'L')[2]
array_c_longi[:,1] = extract_3temporal_shifts_between_2geophones(7,10,'L')[2]
array_c_longi[:,2] = extract_3temporal_shifts_between_2geophones(10,13,'L')[2]
array_c_transv = np.zeros((3,3))
array_c_transv[:,0] = extract_3temporal_shifts_between_2geophones(4,7,'T')[2]
array_c_transv[:,1] = extract_3temporal_shifts_between_2geophones(7,10,'T')[2]
array_c_transv[:,2] = extract_3temporal_shifts_between_2geophones(10,13,'T')[2]

plt.figure()
for i in range(3):
    plt.plot(array_c_longi[i,:],'ob')
plt.title('vitesse des ondes longitudinales pour paires de geophones successifs')
plt.ylabel('c (m/s)')
#plt.ylim(1500,1900)
plt.show()
plt.figure()
for i in range(3):
    plt.plot(array_c_transv[i,:],'or')
plt.title('vitesse des ondes transverse pour paires de geophones successifs')
#plt.ylim(1500,1900)
plt.ylabel('c (m/s)')
plt.show()