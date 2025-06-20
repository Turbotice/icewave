# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 17:07:17 2025

@author: sebas
"""

import numpy as np
import pickle

#%% Set up Febus DAS parameters 

dates = ['0210','0211','0212','0217']
param = {}

for date in dates:
    param[date] = {'fiber_length':[],'fs':[],'facq_x':[]}

# set parameters for 0210
param['0210']['fiber_length'] = 750
param['0210']['fs'] = 840
param['0210']['facq_x'] = 2

# set parameters for 0211
param['0211']['fiber_length'] = 700
param['0211']['fs'] = 84
param['0211']['facq_x'] = 0.5

# set parameters for 0212
param['0212']['fiber_length'] = 700
param['0212']['fs'] = 840
param['0212']['facq_x'] = 1.0

# set parameters for 0217
param['0217']['fiber_length'] = 700
param['0217']['fs'] = 840
param['0217']['facq_x'] = 5.0

param['units'] = {}
param['units']['fiber_length'] = 'meters'
param['units']['fs'] = 'Hz'
param['units']['facq_x'] = 'pix/meter'

path2save = 'U:/Data/'
file2save = f'{path2save}parameters_Febus_2025.pkl'

with open(file2save,'wb') as pf:
    pickle.dump(param,pf)






