# -*- coding: utf-8 -*-
"""
Created on Fri May 23 14:29:27 2025

@author: sebas
"""


import numpy as np
import re
import pickle

import glob 

table = {}

path2data = 'U:/Aurore_frasil/'
folderlist = glob.glob(f'{path2data}*mm_laser')

scale = 0.337 # scale in mm/pix
facq_pix = 1/(scale*1e-3) # scale in pixels / meter

for folder in folderlist:
    h = float(re.findall(r'e_(\d+\.\d+)mm',folder)[0]) # frasil thickness
    print(h)
    
    explist = glob.glob(f'{folder}/*Hz_*mm')
    for exp in explist:
        f_ex = float(re.findall(r'(\d+\.\d+)Hz',exp)[0]) # excitation_frequency 
        amplitude = float(re.findall(r'Hz_(\d+)mm',exp)[0])
        
        key = f'h_{h}mm_fex_{f_ex}Hz'
        print(key)
        param_dict = {}
        param_dict[key] = {'h':h,'f_ex': f_ex ,'theta':53.9, 'facq_t':37.215, 'scale':scale, 'facq_pix': facq_pix}
        print(param_dict)
        table.update(param_dict)

# Correct table of parameters for h = 2.5mm or h = 5.0mm and f_ex <= 4.2 Hz
for key in table.keys():
    h = table[key]['h']
    f_ex = table[key]['f_ex']
    
    test = (h == 2.5) or (h == 5.0 and f_ex < 4.2)
    if test:
        table[key]['theta'] = 44.17
        table[key]['facq_t'] = 67

#%% Save table 

file2save = f'{path2data}table_experimental_parameters.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(table,pf)


