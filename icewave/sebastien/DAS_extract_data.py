# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:41:26 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
import h5py 
import glob

import icewave.tools.matlab2python as mat2py
#%%

date = '0210'
path2data = f'E:/Data/{date}/DAS/'

filelist = glob.glob(path2data + '*.h5')
i = 1
file2load = filelist[i]
with h5py.File(file2load,'r') as f:
    print(list(f.keys()))
    
    a_group_key = list(f.keys())[0]
    data = mat2py.mat_to_dict(f[a_group_key], f[a_group_key])

#%% 

strain_rate = data['Source1']['Zone1']['Strain Rate [nStrain|s]']
t = data['Source1']['time']

