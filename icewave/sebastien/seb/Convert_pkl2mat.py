# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:32:25 2023

@author: sebas
"""
#%%
import pickle
import scipy.io
import glob
import os

#%%
# def pkl2mat(file):
#     """ Converts a .pkl file to a .mat file and save it in the same 
#     directory as the initial pkl file """
#     file_pkl = open(file,'rb')
#     d

#%%
file = 'W:/SagWin2024/Data/0211/BoueeVague/Resultats/*.pkl'

# save_path = 'C:/Users/sebas/Stage_MIZ/Datas_mat/'
file_pkl = glob.glob(file)

# if os.path.isdir(save_path) == False:
#     os.mkdir(save_path)

#%%

a = file_pkl[0]
b = a.strip('pkl')
print(b)

#%%
for f in file_pkl :
    p = open(f,'rb')
    data = pickle.load(p)
    dict = {}
    dict['data'] = data
    
    f_mat = f.strip('pkl') + 'mat'
    print(f_mat)
    scipy.io.savemat(f_mat,dict)

#%%
name = 'probe_features'  #probe without 1
p = open('%s.pkl'%name,'rb')
data = pickle.load(p)
dict = {}
dict['data'] = data
scipy.io.savemat(name,dict)