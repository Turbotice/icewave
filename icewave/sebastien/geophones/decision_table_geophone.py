# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:05:30 2024

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os 

import seb
from seb import pickle_m

#%% Initial parameters 

g = 9.81
nu = 0.33
h_ice = np.arange(0.1,2.1,0.1)
rho_w = rho_w = 1027 # density of water 
rho_ice = 917
E = np.arange(2e9,11e9,1e9) # Young modulus 

# D = E*pow(h_ice,3)/12/(1-nu**2)
# l_D = pow(D/rho_w/g,0.25) # gravito-elastice length 
# lambda_c = 2*np.pi*l_D

#%% Create matrix with several values of h_ice and E 

row_h = np.reshape(h_ice,(np.size(h_ice),1))
col_E = np.reshape(E,(1,np.size(E)))

D = col_E*pow(row_h,3)/12/(1-nu**2)
l_D = pow(D/rho_w/g,0.25) # gravito-elastic length 
lambda_c = 2*np.pi*l_D

#%% Write table in a csv file 
path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Arctic_refuge_2024/Expedition_plan/'
csv_file = path + 'Decision_table_geophones_small_test.csv'

with open(csv_file, 'w', newline = '') as csvfile: 
    writer = csv.writer(csvfile) 

    row = ['h_ice\E']  + list(np.round(E*1e-9,2))
    writer.writerow(row)

    for i, h_ice_value in enumerate(h_ice):
        row = [round(h_ice_value,3)] + list(np.round(lambda_c[i,:],1))
        writer.writerow(row)

#%% 
a = np.arange(0,5)
b = np.arange(0,3)

row = np.reshape(a,(np.size(a),1))
col = np.reshape(b,(1,np.size(b)))

A = row*col 
#%%
fig, ax = plt.subplots()

ax.imshow()



