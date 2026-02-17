import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from read_trisonica import *
from compute_flux import *
from read_anemo import *


# Notes:
# - SAR images 29/02, 3/02, 10/02, 15/02


test_tri = False
test_thies = True


file_path="/home/jacqhugo/BicWin26/terrain/0206/Wind/"
file="20260205_1558_4.txt"

file_path="/home/jacqhugo/BicWin26/terrain/0212/Mat_portatif/Trisonica/"
file="serial_20260212_103030.txt"
file="serial_20260212_103339.txt"
file="serial_20260212_185053.txt"


file_path="/home/jacqhugo/BicWin26/terrain/0213/Mat_portatif/"
gps_file="Gobannos/2026-02-13T12_57_36-gps-1-610609801644-611509640581.csv"
tri_file="Trisonica/serial_20260213_125151.txt"

thies_file='/home/jacqhugo/BicWin26/terrain/0213/Autruche/Thies/20260213_2244_thies.txt'

date='12/02/2026'

#ds =  parse_trisonica_to_xarray(file_path+tri_file, reference_date=date, verbose=True)

if test_thies:
    ds = parse_anemo_to_xarray(thies_file, 
                        verbose=True,
                        log_errors=True,
                        check_gaps=True,
                        kind='thies',
                        sampling_rate_hz=20)
if test_tri:
    ds = parse_anemo_to_xarray(file_path+tri_file, 
                        verbose=True,
                        log_errors=True,
                        check_gaps=True,
                        kind='trisonica',
                        sampling_rate_hz=5)



# Thies: todo = remove '+' sign that can sometimes be in front of numerics


# get gps/accelerometers data

# ds = compute_mean_X(ds, 'U', period=600)

# Note:
# - enlever les spikes
#
#
# - Convertir température sonique en température réelle
# - calculer q à partir de HR (attention pas possible d'avoir les fluctuations ...)
# - periode de moyenne: voir HDR JE-Sicart figure 1.23
#    et tester à posteriori si la moyenne est adaptée
# - rotation du système pour avoir <V> = <W> = <v'w'> = 0

# - récupérer les caractéristiques des anémo (height, location)
#
# - récupérer le context (gobanos+notes)

# First look
if False:
    fig, ax = plt.subplots(2,1,figsize = (3,3), constrained_layout=True, dpi=100,
                        height_ratios=[3,1])
    ax[0].set_ylabel('Wind speed (m/s)')
    ax2= ax[0].twinx()
    ax2.plot(ds.time, ds.D, label='d', c='orange')
    ax[0].plot(ds.time, ds.S, label='S', c='blue')
    ax2.spines['left'].set_color('blue')
    ax2.spines['right'].set_color('orange')
    ax2.set_ylim([0,360])
    ax[0].set_ylim([0,12])
    ax2.set_ylabel('Direction (°) from north')
    ax[0].set_xlabel('time '+date)
    ax[0].set_title('LI-550 portable')

# Format the x-axis for time
    ax[0].xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Every 5 minutes
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels


fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
ax.set_ylabel('U (m/s)')
ax.scatter(ds.time, ds.U,  c='b', marker='+')
#ax.set_ylim([0,12.0])
ax.set_xlabel('time '+date)
ax.set_title('Thies')






plt.show()

