# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:17:25 2023

@author: Banquise
"""


import numpy as np
import matplotlib.pyplot as plt

import baptiste.files.save as sv
import baptiste.display.display_lib as disp


import baptiste.math.fits as fit
import baptiste.math.RDD as rdd

import icewave.geophone.fonctions_extract_triangles_flexion as fct

params = {}

#%%



date = '20230313'
# folder ='/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/labshared2/Banquise/Rimouski_2023/Data/Geophones/'+date+'/'
# savefolder = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/labshared2/Banquise/Rimouski_2023/Traitements_donnees/baptsarahantonin'

folder = 'W:\Banquise/Rimouski_2023/Data/Geophones\\' +date+'\\'
savefolder = 'W:\Banquise\Rimouski_2023\Traitements_donnees\\baptsarahantonin\\'


data=fct.import_files(date,1, folder, savefolder)

#%%


params['voie']='Z'
params['num_geo'] = '4-6'


params = fct.direction(params)

#%%

plt.figure()
plt.plot(data[:,2])
plt.plot(data[:,5])
plt.plot(data[:,8])

