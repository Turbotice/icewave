# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:03:44 2024

@author: sebas

This script describes how to convert PIV data saved as .mat files using the function PIVmat2pkl

"""

import numpy as np
import pickle 
import h5py
import os 

import seb 
from seb import PIVmat2pkl
path = 'Y:/Banquise/Aurore/e_5mm/Data/'
filelist = os.listdir(path) # get all files/folders in the path 

for file in filelist: 
    folder = os.path.join(path,file)
    matfile_list = os.listdir(folder)
    for matfile in matfile_list: # loop over all files of the folder
        if matfile.endswith('processed.mat'):
            path2data = os.path.join(folder,matfile)
            print(path2data)
            PIVmat2pkl.PIVmat2pkl(path2data)
    

