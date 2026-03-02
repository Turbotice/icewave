#%%

import numpy as np
from random import random
import matplotlib.pyplot as plt 
import pickle
import os , sys
import time
from scipy.interpolate import interp1d
import seaborn as sns


#%%-------------------------------LOADING DATA IN-PLANE--------------------------------------------------------------------------

year = '2025'
date = '0304' #date format, 'mmdd'
acqu_numb = '0001' #acquisition number 
ordi = 'dell_vasco'

def load_inplane_data(year='2025',date='0304',acqu_numb='0002',ordi='dell_vasco', rho_ice=800):

    if ordi=='adour':
        path2data = f'/media/turbots/Backup25/Data/{date}/Geophones/'
    elif ordi=='babasse':
        path2data = os.path.join('E:/Data/',date,'Geophones/')
    elif ordi=='dell_vasco':
        #path2data = f'D:/copie_BicWin25_geophones/Data/{date}/Geophones/'
        path2data = f'B:/Data/{date}/Geophones/'
    # path2data = 'C:/Users/sebas/icewave/icewave/sebastien/geophones/updatescriptspython/0211/Geophones/'
    # path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/' +year+'_BICWIN/'



    file2load2 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir.pkl'
    file2load3 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir.pkl'

    # Load phase velocity of in-plane waves 
    with open(file2load2, "rb") as f:
        data = pickle.load(f)
    cQS0 = data 
    with open(file2load3, "rb") as f:
        data = pickle.load(f)
    cSH0 = data 

    data_inplane = {
        'cQS0': cQS0,
        'cSH0': cSH0
    }


    nu = 1-2*(data_inplane['cSH0']/data_inplane['cQS0'])**2
    E = rho_ice*data_inplane['cQS0']**2*(1-nu**2)
    print(f'Young modulus, E = {E*1e-9} and Poisson coefficient, nu = {nu}')
    
    if len(E)==1:
        E=E[0]
        nu=nu[0]

    return data_inplane, (E, nu, rho_ice)
#%% Find a first guess of E and nu, assuming a given density for ice 
disk = 'B:'
path2fileinplaneacq = f'{disk}/General/Summary_geophone_lines/inversions/acquisitions_nb_SH0_QS0_only.txt'

datainplaneacq = np.loadtxt(path2fileinplaneacq,dtype=str)

all_data_inplane = []
E_arr = []
nu_arr = []
rho_ice_arr = []
date_arr = []
acqu_numb_arr = []
for i in range(len(datainplaneacq)):
    date = datainplaneacq[i][:4]
    acqu_numb = str(int(datainplaneacq[i][7:])).zfill(4)

    date_arr.append(date)
    acqu_numb_arr.append(acqu_numb)
    print(date, acqu_numb)
    result = load_inplane_data(year='2025', date=date, acqu_numb=acqu_numb)
    all_data_inplane.append(result[0])
    E_arr.append(result[1][0])
    nu_arr.append(result[1][1])
    rho_ice_arr.append(result[1][2])
    

# %% save csv (pas encore fait)
