import matplotlib.pyplot as plt
from read_marmotte import *


# Get the file with the command:
# scp -r -P 8022 u0_a236@192.168.1.184:/data/data/com.termux/files/home/git/pytruche/codes/marmotte .

# file
file_path = "/home/jacqhugo/BicWin26/icewave/icewave/wind/marmotte/"
filename = "20260217_141524.203_marmotte.txt"

file_path="/home/jacqhugo/BicWin26/terrain/calib_thermique/marmotte/"
filename="20260217_175931.335_marmotte.txt"


# .txt -> .nc
ds = parse_marmotte(file_path+filename)

ds.to_netcdf(file_path+filename[:-3]+"nc")

# plot
fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
ax.plot(ds.time, ds.IRSensor_T, label='Ts')
ax.plot(ds.time, ds.Airtemp, label='Airtemp')
ax.plot(ds.time, ds.IRTemp_C, label='Ts corrected' )
ax.set_ylabel('°C')
ax.set_xlabel('Time')
plt.legend()

plt.show()
