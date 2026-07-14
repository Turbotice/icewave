# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 15:33:52 2026

@author: sebas

Collect thickness drilling performed during BicWin25
"""


import icewave.field.gps as field_gps
import pickle

dates = ['0210','0211','0212']

records = {}
for date in dates :
    record = field_gps.get_records(date,year = '2025')
    records[date] = record['gps']['garminSP']
    
print(records)

path2save = '/media/turbots/DATA/thiou/labshared1/Banquise/Sebastien/'
file2save = f'{path2save}records_gps_BicWin25.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(records,pf)
    
    