"""
Functions to compte a surface temperature from a flux measured by the sensor IR120

This is following the documentation of the IR120 by Campbell
https://s.campbellsci.com/documents/eu/manuals/ir100_ir120%20-%20708.pdf

If you use your own IR120, you will need to adapt the callibration coefficient
"""


import numpy as np
import matplotlib.pyplot as plt

def compute_Rt(Va, V0, R1):
    return R1*Va/(V0-Va)

def compute_Tt(R): 
    """
    Thermistor temperature (Celsius)
    """

    # Calibration constants
    A = 9.325752E-4
    B = 2.212819E-4
    C = 1.333627E-7
    
    return 1/(A+B*np.log(R)+C*np.log(R)**3) - 273.15

def temperature_compensation(IRSensor_Volt, IRSensorCan_temp):
    return IRSensor_Volt * 1.0004**(IRSensorCan_temp - 25)

def compute_energy(IRSensorVolt_TC):
    # Calibration constants
    X = 8.015848E-6
    Y = 2.985651E-1
    Z = 3.024548E-1

    return (X*IRSensorVolt_TC**2 
            + Y*IRSensorVolt_TC
            + Z)

def add_body_sensor_E(IRSensor_E, IRSensorCan_Temp):
    """
    Takes the absolute energy measured, and correct it from body temp.
    Return emission from the surface (W/m2)
    """
    sig = 5.67e-8 # stefan boltzman constant
    return (IRSensor_E + sig * (IRSensorCan_Temp + 273.15)**4)


def Corrected_temp(IRSensor_T4, Airtemp, emissivity=0.95):
    """
    Computes a temperature in Celsius, correct from air temperature and emissivity
    """
    return ((IRSensor_T4 - ( (Airtemp + 273.15)**4*(1-emissivity) ))/emissivity )**0.25 - 273.15

def bits_to_tension(D, N, Alim):
    return D/N*Alim


# Vectorized versions
vect_bits_to_tension = np.vectorize(bits_to_tension, excluded=["arg1", "arg2"])
vect_compute_Rt = np.vectorize(compute_Rt, excluded=['arg1', 'arg2'])
v_compute_Tt = np.vectorize(compute_Tt)
v_temperature_compensation = np.vectorize(temperature_compensation)
v_compute_energy = np.vectorize(compute_energy)
v_add_body_sensor_E = np.vectorize(add_body_sensor_E)
v_Corrected_temp = np.vectorize(Corrected_temp, excluded=['arg2'])
