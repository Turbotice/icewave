"""
This script test the function in 'sensorIR120'.
Each step is detailled, see sensorIR120.py for more details
"""

import numpy as np
import matplotlib.pyplot as plt
from sensorIR120 import *

emissivity = 0.95

V = 400 - 512
V0 = 1024
N = 1024
Va = 740
R1 = 10e3 # Ohms

Alim = 5000 # mV
TensionV = V/N*Alim
TensionVa = Va/N*Alim
TensionV0 = V0/N*Alim

print(f'TensionV {TensionV} (mV)')
print(f'TensionVa {TensionVa} (mV)')

Rt = compute_Rt(Va, V0, R1)
print(f'Rt {Rt/1e3} kOhms')

IRSensorCan_temp = compute_Tt(Rt)
print(f'IRSensorCan_temp {IRSensorCan_temp} (°C)')

IRSensorVolt_TC = temperature_compensation(TensionV, IRSensorCan_temp)
print(f'IRSensorVolt_TC {IRSensorVolt_TC} (mV)')

IRSensor_E = compute_energy(IRSensorVolt_TC)
print(f'IRSensor_E {IRSensor_E} (W.m-2)')

LW_surf = add_body_sensor_E(IRSensor_E, IRSensorCan_temp)
print(f'IRSensor_Es {LW_surf} (W.m-2)')

IRSensor_T4 = LW_surf/(5.67e-8)


IRSensor_T = (IRSensor_T4**0.25) - 273.15 # °C
print(f'IRSensor_T {IRSensor_T} (°C)')

Airtemp = IRSensorCan_temp
Airtemp = -8
IRTemp_C = Corrected_temp(IRSensor_T4, Airtemp, emissivity) # °C
print(f'IRTemp_C {IRTemp_C} (°C)')

# Plot
R = np.arange(50, 200)
fig, ax = plt.subplots(1,1,figsize = (3,3), constrained_layout=True, dpi=100)
ax.plot(R, compute_Tt(R*1000))
ax.hlines(0, 0, 1000, color='gray', ls='--')
ax.set_ylabel('T (°C)')
ax.set_xlabel('R (kOhms)')
ax.set_xlim([50,200])
plt.show()
