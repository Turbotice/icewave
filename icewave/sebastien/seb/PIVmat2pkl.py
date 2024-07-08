# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:17:03 2024

@author: sebas
"""

import numpy as np
import pickle 
import h5py
# import seb
# from seb import pickle_m

def PIVmat2pkl(path,relevant_keys = ['Dt','Vx','Vy','W']):
    """ Convert PIV data obtained from PIVlab. Data saved in a .mat file is converted into a .pkl file 
    Input : - path : path to a .mat file
            - relevant_keys : (optional argument) keys that are kept from the .mat file (tables of parameters s and p are automatically saved)
    """
    dataset = {}
    with h5py.File(path, 'r') as f:

    # Create a dictionary from h5py file 

        for key in f.keys():
            print(key)
            if key in relevant_keys:
                dataset[key] = np.squeeze(np.array(f[key]))
            
        # Create new array of PIV parameters p and s 
        
        param = {}
        s_shape = np.shape(f['s'])
        param['s'] = {}
        s_headers = ['Int. area 1','Step size 1','Subpix. finder','Mask','ROI','Nr. of passes','Int. area 2','Int. area 3',
                     'Int. area 4','Window deformation','Repeated Correlation','Disable Autocorrelation','Correlation style',
                     'Repeat last pass','Last pass quality slope']
        
        for i in range(0,s_shape[1]):
            param['s'][s_headers[i]] = np.squeeze(f[f['s'][1,i]])
            
        param['s']['Window deformation'] = np.array('*spline') 
        
        p_shape = np.shape(f['p'])
        param['p'] = {}
        p_headers = ['ROI','CLAHE','CLAHE size','Highpass','Highpass size','Clipping','Wiener','Wiener size','Minimum intensity','Maximum intensity']
        
        for i in range(0,p_shape[1]):
            param['p'][p_headers[i]] = np.squeeze(f[f['p'][1,i]])


    dataset['s'] = param['s']
    dataset['p'] = param['p']

    dataset['Vx'] = np.transpose(dataset['Vx'],(1,2,0))
    dataset['Vy'] = np.transpose(dataset['Vy'],(1,2,0))
    
    pickle_file = path.replace('mat','pkl')
    with open(pickle_file,'wb') as file:
        pickle.dump(dataset,file)
    print('Data loaded as : ' + pickle_file)




# #%% Load raw matlab Data 

# # Load data 
# path = 'Y:/Banquise/Aurore/e_5mm/Data/3.5Hz_15mm_e_5mm_sans_laser/'
# filename = 'PIV_processed_i00_Dt1_b1_W32_xROI1_width1491_yROI1_height949post-processed.mat'

# print('Loading matlab data..')
# # mat = io.loadmat(path+filename)

# # f = h5py.File(path + filename,'r')
# dataset = {}
# with h5py.File(path + filename, 'r') as f:

# # Create a dictionary from h5py file 

#     relevant_keys = ['Dt','Vx','Vy','W']
#     for key in f.keys():
#         print(key)
#         if key in relevant_keys:
#             dataset[key] = np.squeeze(np.array(f[key]))
        
#     # Create new array of PIV parameters p and s 
    
#     param = {}
#     s_shape = np.shape(f['s'])
#     param['s'] = {}
#     s_headers = ['Int. area 1','Step size 1','Subpix. finder','Mask','ROI','Nr. of passes','Int. area 2','Int. area 3',
#                  'Int. area 4','Window deformation','Repeated Correlation','Disable Autocorrelation','Correlation style',
#                  'Repeat last pass','Last pass quality slope']
    
#     for i in range(0,s_shape[1]):
#         param['s'][s_headers[i]] = np.squeeze(f[f['s'][1,i]])
        
#     param['s']['Window deformation'] = np.array('*spline') 
    
#     p_shape = np.shape(f['p'])
#     param['p'] = {}
#     p_headers = ['ROI','CLAHE','CLAHE size','Highpass','Highpass size','Clipping','Wiener','Wiener size','Minimum intensity','Maximum intensity']
    
#     for i in range(0,p_shape[1]):
#         param['p'][p_headers[i]] = np.squeeze(f[f['p'][1,i]])


# dataset['s'] = param['s']
# dataset['p'] = param['p']

# dataset['Vx'] = np.transpose(dataset['Vx'],(1,2,0))
# dataset['Vy'] = np.transpose(dataset['Vy'],(1,2,0))

#%% Try to save data 

# pickle_file = path+filename
# pickle_file = pickle_file.replace('mat','pkl')
# # pickle_m.write(dataset,pickle_file)

# with open(pickle_file,'wb') as file:
#     pickle.dump(dataset,file)
#     print('Data loaded as : ' + pickle_file)
# #%% Try to load data 

# # new_data = pickle_m.read(pickle_file)

# with open(pickle_file,'rb') as file:
#     new_data = pickle.load(file)




