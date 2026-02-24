"""
This scripts computes the surface temperature from the IR120 sensor


Note:
- We assume Air temperature = body temperature 
  but we can use other way to measure the real air temp*

TODO:
- fix the -512 that is added in the Arduino program
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime
from sensorIR120 import *


def extract_value(s, key):
    """Extract the number after 'key : ' in string s."""
    idx = s.find(key)
    if idx == -1:
        return None
    after = s[idx + len(key):]
    return float(after.split(",")[0].strip())


def parse_marmotte(file):
    print(file)
    # Values of the electric wiring
    R1 = 10e3 # Ohms
    V0 = 5e3 # mV
    N = 2**(10) # discretization

    # Physic constants
    sigma = 5.67e-8 # Stefan_Boltzman constant
    emissivity = 0.95


    data = {"time": [], "flux": [], "check": [], "R": []}
    with open(file) as f:
        for line in f:
            # Only process timestamped lines with all 3 variables
            if not line.startswith("[") or "Flux" not in line or "Check" not in line:
                continue
            
            # Parse timestamp
            timestamp_str = line[1:line.find("]")]
            timestamp = datetime.strptime(timestamp_str, "%Y-%m-%d %H:%M:%S.%f")
            
            # Parse values
            flux  = extract_value(line, "Flux (V) : ")
            check = extract_value(line, "Check (V) : ")
            r     = extract_value(line, "R (V) : ")
            
            if None in (flux, check, r):
                continue  # skip incomplete lines
            
            data["time"].append(timestamp)
            data["flux"].append(flux)
            data["check"].append(check)
            data["R"].append(r)

    ds = xr.Dataset(
        {
            "flux":  ("time", np.array(data["flux"])),
            "check": ("time", np.array(data["check"])),
            "R":     ("time", np.array(data["R"])),
        },
        coords={"time": np.array(data["time"], dtype="datetime64[ms]")},
    )

    # To fix: -512 
    ds['flux'] = ds['flux']- 512

    # Transform to real tensions
    ds['TensionFlux'] = ("time", vect_bits_to_tension(ds["flux"].values, N=N, Alim=V0))
    ds['TensionCheck'] = ("time", vect_bits_to_tension(ds["check"].values, N=N, Alim=V0))
    ds['TensionRt'] = ("time", vect_bits_to_tension(ds["R"].values, N=N, Alim=V0))

    # Transform to Ts, Rt
    ds['Rt'] = ("time", vect_compute_Rt(ds["TensionRt"].values, V0=V0, R1=R1))
    ds['IRSensorCan_temp'] = ("time", v_compute_Tt(ds['Rt'].values))
    ds['IRSensorVolt_TC'] = ("time", v_temperature_compensation(ds['TensionFlux'].values,
                                                                ds['IRSensorCan_temp'].values))
    ds['IRSensor_E'] = ("time", v_compute_energy(ds['IRSensorVolt_TC'].values))
    ds['LW_surf'] = ("time", v_add_body_sensor_E(ds['IRSensor_E'].values,
                                                ds['IRSensorCan_temp'].values))
    ds['IRSensor_T4'] = ds['LW_surf']/sigma
    ds['IRSensor_T'] = (ds['IRSensor_T4']**0.25) - 273.15

    ds['Airtemp'] = ds['IRSensorCan_temp'] # ASSUMPTION
    ds['IRTemp_C'] = ("time", v_Corrected_temp(ds['IRSensor_T4'].values,
                                            ds['Airtemp'].values,
                                            emissivity=emissivity))



    return ds


