import pandas as pd
import os
import numpy as np
"""
gps_num = 1

date = '20260210'

waypoint = '049'
"""

def open_temperatures_csv(date, gps_num, waypoint): 
    path2data = f'F:/bicwin26/Summary/temperature_profiles/'
    filename = path2data + f'{date}_gps{gps_num}_waypoint{waypoint}.csv'

    df = pd.read_csv(filename,sep=';')

    return df