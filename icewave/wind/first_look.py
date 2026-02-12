import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from read_trisonica import *

file_path="/home/jacqhugo/BicWin26/terrain/0206/Wind/"
file="20260205_1558_4.txt"

file_path="/home/jacqhugo/BicWin26/terrain/0212/Mat_portatif/Trisonica/"
file="serial_20260212_103030.txt"
file="serial_20260212_103339.txt"
file="serial_20260212_185053.txt"

date='12/02/2026'

ds =  parse_log_to_xarray(file_path+file, reference_date=date, verbose=True)




# Note:
# - enlever les spikes
#
#
# - Convertir température sonique en température réelle
# - calculer q à partir de HR (attention pas possible d'avoir les fluctuations ...)
# - periode de moyenne: voir HDR JE-Sicart figure 1.23
#    et tester à posteriori si la moyenne est adaptée
# - rotation du système pour avoir <V> = <W> = <v'w'> = 0

# First look
fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=100)
ax.set_ylabel('Wind speed (m/s)')
ax2= ax.twinx()
ax2.plot(ds.time, ds.D, label='d', c='orange')
ax.plot(ds.time, ds.S, label='S', c='blue')
ax2.spines['left'].set_color('blue')
ax2.spines['right'].set_color('orange')
ax2.set_ylim([0,360])
ax.set_ylim([0,10])
ax2.set_ylabel('Direction (°) from north')
ax.set_xlabel('time '+date)
ax.set_title('LI-550 portable')

# Format the x-axis for time
ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Every 5 minutes
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))  # Format as HH:MM:SS
plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')  # Rotate labels





plt.show()

