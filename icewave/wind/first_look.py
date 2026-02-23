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
from sensors import *

# Notes ____________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
"""
> TRANSFER FILES

- find files:
  lz **/Mat_portatif/Trisonica/*.txt (on zsh, not bsh)




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
sensor = 'thies'        # 'trisonica' or 'thies'

mean_kwargs={'method':'block', # 'block' or 'moving'
             'period':600,      # (s), Time average period
             }
# avg_method = 'moving'        
# avg_period = 600            
verbose = True              # speaks more
log_errors = True           # in the parsing
check_gaps = True           # in the original file
save_nc = True              # save dataset as nc in same location
force_save = False          # overwrite the .nc if already here

# > FILE PATHS
# ToDo:
#  get all the file in 1 folder, for each sensor

HOME = "/home/jacqhugo"
# file_path=f"{HOME}/BicWin26/terrain/0206/Wind/"
# file="20260205_1558_4.txt"

# file_path=f"{HOME}/BicWin26/terrain/0212/Mat_portatif/Trisonica/"
# file="serial_20260212_103030.txt"
# file="serial_20260212_103339.txt"
# file="serial_20260212_185053.txt"

# file_path=f"{HOME}/BicWin26/terrain/0213/Mat_portatif/"
# gps_file="Gobannos/2026-02-13T12_57_36-gps-1-610609801644-611509640581.csv"
# tri_file="Trisonica/serial_20260213_125151.txt"

file_path=f"{HOME}/BicWin26/terrain/0220/Mat_portatif/"
tri_file="Trisonica/serial_20260220_124915.txt"

#thies_file=f"{HOME}/BicWin26/terrain/0213/Autruche/Thies/20260213_2244_thies.txt"
thies_file=f"{HOME}/BicWin26/terrain/0218/Autruche/Thies/20260218_0737_thies.txt"

# Preprocessing _____________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

# > initialize sensors

w_Thies = Anemometer(fs=20, 
                     files=[thies_file],
                     name='thies', 
                     description='Top of the mast', 
                     notes='')
w_Trisonica_portable = Anemometer(fs=5, 
                                  files=[file_path+tri_file],
                                  name='trisonica', 
                                  description='Front of the edge', 
                                  notes='')
T_IR120 = T_sensor(fs=1, 
                   files=[],
                   name='marmotte', 
                   kind='IR120')

if sensor=='thies':
    anemo = w_Thies
elif sensor=='trisonica':
    anemo = w_Trisonica_portable

# > Data Parsing
file_to_parse = anemo.files[0] 
kind = anemo.name
sampling = anemo.fs
mean_kwargs['sampling'] = sampling

file_nc = parse_anemo_to_xarray(file_to_parse, 
                    verbose=True,
                    log_errors=True,
                    check_gaps=True,
                    sampling_rate_hz=sampling,
                    save_nc=save_nc,
                    force_save=force_save,
                    anemo = anemo)

# > Motion correct of the mast

# > GPS Data

# > Footprint analysis (see Aubinet et al. 2012 Chapter 8 "Eddy Covariance: A
# Practical Guide to Measurement and Data Analysis")

# > Tilt correction 
compute_and_save(file_nc, 
                update_nc=True,
                operator=rotation_operator,
                 mean_kwargs=mean_kwargs)

# COMPUTES __________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


# > Mean
for var in ['U','V','W','T']:
    compute_and_save(file_nc,
                    update_nc=True,
                    operator=mean_operator,
                    X=var, 
                    period=mean_kwargs['period'], 
                    sampling=mean_kwargs['sampling'],
                    method=mean_kwargs['method'])


# > Flux (eddy covariance)



# > u* and z0
#
# See Foken et al. 2024 "Micrometeorology" part 4
#
#   - profile method (4.1)
#   - direct using eddy covariance (4.2)



# PLOTS _____________________________________________
# ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾

ds = xr.open_dataset(file_nc)

# > First look
if False:
    if sensor=='trisonica':
        fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
        ax.set_ylabel('Wind speed (m/s)')
        ax.plot(ds.time, ds.U1, label='U', c='blue')
        ax.plot(ds.time, ds.mean_U1, label='<U>', c='navy')
        ax.set_ylim([-5,10])
        ax.set_title('LI-550 portable')

        # Format the x-axis for time
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Every 5 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels

    if sensor=='thies':
        fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
        ax.set_ylabel('U (m/s)')
        ax.scatter(ds.time, ds.U1,  c='b', marker='+', label='U')
        ax.plot(ds.time, ds.mean_U1, label='<U>', c='navy')
        ax.set_title('Thies')

        # Format the x-axis for time
        ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=30))  # Every 5 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels
    

    date = pd.Timestamp(ds.time.values[0]).strftime('%d/%m/%Y')
    ax.set_xlabel(f'time {date}')
    plt.legend()

# > test the rotation part
if True:
    fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
    ax.plot(ds.time, ds.U1, label='U1', ls='-')
    ax.plot(ds.time, ds.U2, label='U2', ls='-')
    ax.plot(ds.time, ds.U, label='U',ls='--')
    ax.set_xlabel('time')
    ax.set_ylabel('')
    plt.legend()



plt.show()

