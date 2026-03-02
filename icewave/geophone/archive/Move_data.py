# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 00:08:40 2024

@author: moreaul
"""

#Bouge les datas des geophones au disque dur

import os
# import baptiste.files.file_management as fm
import shutil as sh
import pandas as pd
import numpy as np


path_DD = 'C:\\Canada_2024\\' #ORDI LUDO

geo = ['F:\\','E:\\','G:\\']

metadata = {}

def save_geo (path_geo, path_save, metadata, num_exp) :
    files = os.listdir(path_geo)
    for file in files :
        if ".MiniSeed" in file :
            if int(file[4:7]) == num_exp :
                sh.copy(path_geo + file, path_save + metadata['Serial_number'] + '.0' + file[4:7] + '.' + metadata['annee'] + '.' + 
                        metadata['mois'] + '.' + metadata['jour'] + '.' + metadata['heure'] + '.' + metadata['minute'] + '.' 
                        + metadata['seconde'] +'.000.' + XYZ_to_ZEN(file[7]) + '.Miniseed')
    
def XYZ_to_ZEN (a) :
    if a == 'X':
        return 'N'
    if a == 'Y':
        return 'E'
    if a == 'Z':
        return 'Z'


for path_geo in geo :
    files = os.listdir(path_geo)
    for file in files :
        if 'DigiSolo' in file :
            digi = pd.read_csv(path_geo + file, sep = '\t')
            digi = np.asarray(digi)
            num_exp = 0
            if np.shape(digi) == (0,1) :
                pass
            else :
                for i in range(len(digi)) :
                    
                    if 'Serial Number' in digi[i][0] :
                        metadata['Serial_number' ] = digi[i][0][-9:]
                    if '[Notify00001]' in digi[i][0] :
                        metadata['jour'] = digi[i+1][0][-12:-10]
                        metadata['mois'] = digi[i+1][0][-15:-13]
                        metadata['annee'] = digi[i+1][0][-20:-16]
                        metadata['heure'] = digi[i+1][0][-9:-7]
                        metadata['minute'] = digi[i+1][0][-6:-4]
                        metadata['seconde'] = digi[i+1][0][-3:-1]
                        
                        path_save = path_DD + metadata['annee'] + "\\" + metadata['mois'] + metadata['jour'] + "\\Geophones\\"
                        save_geo(path_geo, path_save + 'bon_format\\', metadata, num_exp)
                        num_exp += 1
                        
                sh.copy(path_geo + file, path_save + 'DigiSolo_' + metadata['Serial_number'][-4:] + '.txt')

        # if ".MiniSeed" in file :

        #     sh.copy(path_geo + file, path_DD + Serial_number + '.0' + file[4:7] + '.' + annee + '.4' + 
        #             mois + '.' + jour + '.' + heure + '.' + minute + '.' + seconde +'.000.' + file[7])
        
    
    

