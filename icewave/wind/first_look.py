"""
Usage
python3 first_look.py

February 2026 BicWin26 Hugo Jacquet
"""

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from read_trisonica import *
from compute_flux import *
from read_anemo import *


# Notes ____________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
> PREPROC
- enlever les spikes


> CONTEXT
- SAR images 29/02, 3/02, 10/02, 15/02
- get gps/accelerometers data
- récupérer les caractéristiques des anémo (height, location)
- récupérer le context (gobanos+notes)

> PHYSIC
- Convertir température sonique en température réelle
- calculer q à partir de HR (attention pas possible d'avoir les fluctuations
  ...)
- periode de moyenne: voir HDR JE-Sicart figure 1.23 et tester à posteriori si
  la moyenne est adaptée
- rotation du système pour avoir <V> = <W> = <v'w'> = 0

"""



# 
# INPUTS____________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾ 
# > CHOOSE SENSOR
sensor = 'trisonica'        # 'trisonica' or 'thies'
verbose = True              # speaks more
log_errors = True           # in the parsing
check_gaps = True           # in the original file
save_nc = True              # save dataset as nc in same location
force_save = False          # overwrite the .nc if already here

# > FILE PATHS
HOME = "/home/jacqhugo"
# file_path=f"{HOME}/BicWin26/terrain/0206/Wind/"
# file="20260205_1558_4.txt"
#
# file_path=f"{HOME}/BicWin26/terrain/0212/Mat_portatif/Trisonica/"
# file="serial_20260212_103030.txt"
# file="serial_20260212_103339.txt"
# file="serial_20260212_185053.txt"
file_path=f"{HOME}/BicWin26/terrain/0213/Mat_portatif/"
gps_file="Gobannos/2026-02-13T12_57_36-gps-1-610609801644-611509640581.csv"
tri_file="Trisonica/serial_20260213_125151.txt"

#thies_file=f"{HOME}/BicWin26/terrain/0213/Autruche/Thies/20260213_2244_thies.txt"
thies_file=f"{HOME}/BicWin26/terrain/0218/Autruche/Thies/20260218_0737_thies.txt"

# > SENSOR INFOS
dict_sensor = {
            'trisonica':{
                'file':file_path+tri_file,
                'kind':'trisonica',
                'sampling':5},
            'thies':{
                'file':thies_file,
                'kind':'thies',
                'sampling':20},
            }

# Preprocessing _____________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


# > Data Parsing
file_to_parse = dict_sensor[sensor]['file']
kind = dict_sensor[sensor]['kind']
sampling = dict_sensor[sensor]['sampling'

ds = parse_anemo_to_xarray(file_to_parse, 
                    verbose=True,
                    log_errors=True,
                    check_gaps=True,
                    kind=kind,
                    sampling_rate_hz=sampling,
                    save_nc=save_nc,
                    force_save=force_save)

# > GPS Data



# COMPUTES __________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

# > Mean

# > Fluctuations

# > Flux



# PLOTS _____________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

# > First look
if True:
    if sensor=='trisonica':
        fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
        ax.set_ylabel('Wind speed (m/s)')
        ax2= ax.twinx()
        ax2.plot(ds.time, ds.D, label='d', c='orange')
        ax.plot(ds.time, ds.S, label='S', c='blue')
        ax2.spines['left'].set_color('blue')
        ax2.spines['right'].set_color('orange')
        ax2.set_ylim([0,360])
        ax.set_ylim([0,12])
        ax2.set_ylabel('Direction (°) from north')
        ax.set_xlabel('time')
        ax.set_title('LI-550 portable')

        # Format the x-axis for time
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Every 5 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels

    if sensor=='thies':
        fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
        ax.set_ylabel('U (m/s)')
        ax.scatter(ds.time, ds.U,  c='b', marker='+')
        ax.set_xlabel('time')
        ax.set_title('Thies')

        # Format the x-axis for time
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=30))  # Every 5 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels







plt.show()

