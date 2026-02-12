# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:04:20 2024

@author: Banquise
"""

import pandas
import pickle 
import numpy as np
import shutil
import os

from pprint import pprint
import icewave.tools.datafolders as df
import glob

#%%


def mkdir_check(path) :
    if os.path.exists(path) :
        pass
        # print('Already exists')
    else :
        os.mkdir(path)
    

base = df.find_path(disk='Hublot24')
new_base = base.split('/')[0] + '/' +  base.split('/')[1] + '/' + 'Repository/'

paths_select = glob.glob(base + 'Summary/*.txt')

for path in paths_select :
    select = pandas.read_csv(path, header= None)
    select = np.asarray(select)[:,0]
    print(select[0].split('/')[1])
    for folder in select :
        # print(base + folder)
        date = folder.split('/')[1]
        drone = folder.split('/')[5]
        folder_name = folder.split('/')[7]
        mkdir_check(new_base + date)
        mkdir_check(new_base + date + '/Drones/')
        mkdir_check(new_base + date + '/Drones/' + drone)
        mkdir_check(new_base + date + '/Drones/' + drone + '/' + folder_name)
        list_files = os.listdir(base + date + '/Drones/' + drone + '/' + folder_name)
        for file in list_files :
            if True : #date == '0226' or date == '0306' :
                if len(file.split('.')) > 1 and file[0] != '.' :
                    # print(file)
                    shutil.copy(base + folder + '/' + file, new_base + folder + '/' + file )
                else : 
                    print(folder + ' ' + file + ' pas copié')
