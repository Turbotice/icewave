#%%

import numpy as np
from random import random
import matplotlib.pyplot as plt 
import pickle
import os , sys
import time
from scipy.interpolate import interp1d
import seaborn as sns


#%%-------------------------------LOADING DATA--------------------------------------------------------------------------
plt.close('all')

year = '2026'
date = '0203' #date format, 'mmdd'
acqu_numb = '0003' #acquisition number 
equation = 'stein'

ordi = 'dell_vasco'
# path to dispersion relation data 

if ordi=='adour':
    path2data = f'/media/turbots/Backup25/Data/{date}/Geophones/'
elif ordi=='babasse':
    path2data = os.path.join('E:/Data/',date,'Geophones/')
elif ordi=='dell_vasco':
    #path2data = f'D:/copie_BicWin25_geophones/Data/{date}/Geophones/'
    path2data = f'F:/bicwin26/data/{date}/Geophones/'
    #path2data = f'B:/Data/{date}/Geophones/'
# path2data = 'C:/Users/sebas/icewave/icewave/sebastien/geophones/updatescriptspython/0211/Geophones/'
# path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/' +year+'_BICWIN/'

# Load (f,k) points of QS dispersion relation 
direction = 1 # 1 ou 2 
filename = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
# file2load1 = path2data  + date + '/Geophones/' + filename
file2load1  = os.path.join(path2data,filename)
with open(file2load1, "rb") as f:
    data = pickle.load(f)
freq1 = data[0] 
kQS1 = data[1] 

# Load (f,k) points of QS dispersion relation
direction = 2 # 1 ou 2 
filename = year + '_' + date + '_acq'+acqu_numb+ 'disp_QS_dir' +str(direction) +'.pkl'
# file2load1 = path2data  + date + '/Geophones/' + filename
file2load1  = os.path.join(path2data,filename)
with open(file2load1, "rb") as f:
    data = pickle.load(f)
freq2 = data[0] 
kQS2 = data[1] 


mini = max(min(freq1),min(freq2))
maxi = min(max(freq1),max(freq2))
freq = np.linspace(mini,maxi , 15)

# Interpolate kQS1 at the new frequency points, without extrapolation
interpolator1 = interp1d(freq1, kQS1)
kQS1_interp = interpolator1(freq)

# Interpolate kQS2 at the new frequency points, without extrapolation
interpolator2 = interp1d(freq2, kQS2)
kQS2_interp = interpolator2(freq)

# Calculate the average values, ignoring NaNs
kQS = np.nanmean([kQS1_interp, kQS2_interp], axis=0)






##################
# rajouté par vasco pour fitte vite fait hydroelastic
def hydroelastic(k,D):
    g=9.81
    rho=900
    omega = np.sqrt(g*k+(D/rho)*(k**5))
    return (1/(2*np.pi))*omega
from scipy.optimize import curve_fit
popt,pcov= curve_fit(hydroelastic,kQS,freq)

##################




# Plot the results
plt.figure(figsize=(12, 6))
plt.plot(kQS1, freq1,  'bo', label='kQS1')
plt.plot( kQS2,freq2, 'go', label='kQS2')
plt.plot(kQS1_interp, freq,  'b--', label='kQS1 interpolated')
plt.plot(kQS2_interp, freq, 'g--', label='kQS2 interpolated')
plt.plot(kQS, freq, 'r--', label='Average kQS')
plt.xlabel(r'$k_{QS} \: \mathrm{(rad.m^{-1})}$')
plt.ylabel(r'$f \: \mathrm{(Hz)}$')


plt.plot(kQS,hydroelastic(kQS,popt[0]))

plt.legend()
plt.title('Comparison and Average of kQS1 and kQS2')
plt.show()

#%% Load data for inversion 
# file that will be saved 
file2save  = path2data + year + '_' + date + '_acq'+acqu_numb+'_bidir_inversion.pkl'
fig_save_path = path2data + '/Figures_inversion_MCMC_T_0p25/' 
if not os.path.isdir(fig_save_path):
    os.mkdir(fig_save_path)

file2load2 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cQS0_bidir.pkl'
file2load3 = path2data + year + '_' + date + '_acq'+acqu_numb+ '_cSH0_bidir.pkl'


# Load phase velocity of in-plane waves 
with open(file2load2, "rb") as f:
    data = pickle.load(f)
cQS0 = data 
with open(file2load3, "rb") as f:
    data = pickle.load(f)
cSH0 = data 

# load dictionnary of phase velocity 
# file2load2 = path2data +'Phase_velocity_dictionnary_acqu_0001_sig_length_0p3.pkl'
# with open(file2load2,'rb') as pfile:
#     dic_phase_velocity = pickle.load(pfile)
# print('Phase velocities loaded')

data = {
    'freq':freq,
    'kQS': kQS,
    'cQS0': cQS0,
    'cSH0': cSH0
}

"""data = {
    'cQS0': cQS0,
    'cSH0': cSH0}
"""
#%% Find a first guess of E and nu, assuming a given density for ice 

rho_ice = 917
nu = 1-2*(data['cSH0']/data['cQS0'])**2
E = rho_ice*data['cQS0']**2*(1-nu**2)


print(f'Young modulus, E = {E*1e-9} and Poisson coefficient, nu = {nu}')

# %%
print(popt,np.sqrt(pcov))
print(((12*(popt)*(1-nu)**2)/E)**(1/3)) # estimation epaisseur
# %%
