#%%
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append('../../..')
from tools import graphes

#%%
disk = "H"
date = '0223'
# Load the data
dir_path = f'{disk}:/data/{date}/Tests_Bresiliens/'
data_path = f'{dir_path}Summary_samples_tests.csv'

df = pd.read_csv(data_path, sep=';', encoding='latin-1')
#%%
def convert_to_float(value):
    if type(value) != str:
        return value
    if type(value) == str and value!='A NOTER':
        return float(value.replace(',', '.'))
    else:
        return np.nan
def convert_to_arrayoffloats(values):
    return np.array([convert_to_float(value) for value in values])


# Display the first few rows of the DataFrame
print(df.head())

arr_acqnum = df['acq_num'].values
arr_L_mm = df['L_mm'].values
arr_L_err_mm = df['L_err_mm'].values
arr_qualite_test = df['qualite_test'].values
arr_Fc_kN = convert_to_arrayoffloats(df['Fc_kN'].values)
arr_D_mm = df['D_mm'].values
arr_T_core_celsius = convert_to_arrayoffloats(df['T_core_celsius'].values)
arr_mass_kg = 1e-3 * df['mass_g'].values
arr_loc_core = df['loc_core'].values

#indices2plot = np.where((arr_qualite_test == 'ok')|(arr_qualite_test=='perfect'))[0]
mask_quality_perfect = (arr_qualite_test=='perfect')
mask_quality_ok = np.array(['ok' in x for x in arr_qualite_test])
mask_quality = mask_quality_ok|mask_quality_perfect

arr_sigma_c = 1e-6 * (2*arr_Fc_kN*1e3)/(np.pi*arr_D_mm*1e-3*arr_L_mm*1e-3)


print(arr_Fc_kN)
plt.figure(figsize=(10, 6))
plt.scatter(arr_L_mm[mask_quality], arr_sigma_c[mask_quality], color='blue', label='Qualité OK')
plt.xlabel('L (mm)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(0,np.max(arr_L_mm[mask_quality])*1.1)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs L')
plt.savefig(f'{dir_path}sigmac_vs_L.pdf', dpi=300)
plt.show()


plt.figure(figsize=(10, 6))
plt.scatter(arr_T_core_celsius[mask_quality], arr_sigma_c[mask_quality], color='orange', label='Qualité OK')
plt.xlabel('T (°C)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(-5, 0)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs T')
plt.savefig(f'{dir_path}sigmac_vs_T.pdf', dpi=300)
plt.show()

# %%


mask_HaHa = np.array(['HaHa' in x for x in arr_loc_core])
mask_Hatee = np.array(['Hatee' in x for x in arr_loc_core])
mask_Capelans_beach = np.array(['Capelans beach' in x for x in arr_loc_core])
mask_Capelans_floating_ice = np.array(['Capelans floating' in x for x in arr_loc_core])


plt.figure(figsize=(10, 6))
plt.scatter(arr_T_core_celsius[mask_quality&mask_HaHa], arr_sigma_c[mask_quality&mask_HaHa], color='tab:blue', label='Ha!Ha!')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Hatee], arr_sigma_c[mask_quality&mask_Hatee], color='tab:orange', label='Hatée')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Capelans_beach], arr_sigma_c[mask_quality&mask_Capelans_beach], color='tab:green', label='Capelans beach')
plt.scatter(arr_T_core_celsius[mask_quality&mask_Capelans_floating_ice], arr_sigma_c[mask_quality&mask_Capelans_floating_ice], color='tab:red', label='Capelans floating ice')

plt.xlabel('T (°C)')
plt.ylabel('Sigma_c (MPa)')
plt.xlim(-5, 0)
plt.ylim(0,np.nanmax(arr_sigma_c[mask_quality])*1.1)
plt.title('Sigma_c vs T')
plt.savefig(f'{dir_path}sigmac_vs_T.pdf', dpi=300)
plt.legend()
plt.show()


dict_fig1 = graphes.defaults_params_errorbar_categories()

dict_fig1['x'] = arr_T_core_celsius
dict_fig1['y'] = arr_sigma_c
dict_fig1['z'] = np.empty(len(arr_T_core_celsius),dtype=object)
dict_fig1['z'][mask_HaHa] = 'HaHa'
dict_fig1['z'][mask_Hatee] = 'Hatee'
dict_fig1['z'][mask_Capelans_beach] = 'Capelans beach'
dict_fig1['z'][mask_Capelans_floating_ice] = 'Capelans floating'
dict_fig1['ylabel'] = '$\sigma_c$ [MPa]'
dict_fig1['xlabel'] = 'T [°C]'
dict_fig1['savefig'] = False
dict_fig1['title'] = 'sigma_c vs T for cores in different locations'
#dict_fig1['symbols'] = {'HaHa':'o','Hatee':'o','Capelans beach':'o','Capelans floating':'o'}


graphes.errorbar_categories(dict_params=dict_fig1)




# %% load fits rigidity cores
with open(r"H:\data\0223\Tests_Bresiliens\results\fits_force_displacement.pkl", 'rb') as f:
    dict_allfits  = pickle.load(f)



# %%

acqnums_2 = np.empty(len(dict_allfits.keys()), dtype=object)
slopes_2 = np.empty(len(dict_allfits.keys()), dtype=object)
intercepts_2 = np.empty(len(dict_allfits.keys()), dtype=object)
slopeserr_2 = np.empty(len(dict_allfits.keys()), dtype=object)
interceptserr_2 = np.empty(len(dict_allfits.keys()), dtype=object)

for i in range(len(dict_allfits.keys())):
    acqnum_key = list(dict_allfits.keys())[i]
    acqnums_2[i] = acqnum_key
    slopes_2[i] = dict_allfits[acqnum_key]['fit_1']['a']
    intercepts_2[i] = dict_allfits[acqnum_key]['fit_1']['b']
    slopeserr_2[i] = np.diag(dict_allfits[acqnum_key]['fit_1']['cov'])[0]
    interceptserr_2[i] = np.diag(dict_allfits[acqnum_key]['fit_1']['cov'])[1]

slopes_2_usi = 1e6 * slopes_2 # units : N/m
slopeserr_2_usi = 1e6 * slopeserr_2
intercepts_2_usi = 1e3 * intercepts_2 # units : N
interceptserr_2_usi = 1e3 * interceptserr_2


common, a_ind, b_ind = np.intersect1d(arr_acqnum, acqnums_2, return_indices=True)

plt.figure()
plt.plot(arr_T_core_celsius[b_ind], slopes_2_usi*16/(3*np.pi*arr_L_mm[b_ind]*1e-3*arr_D_mm[b_ind]*1e-3), 'o')
plt.yscale('log')
# %%
def Weibull_cdf(x, x0=1, m=1):
    return 1 - np.exp(-(x/x0)**m)


plt.figure()
plt.hist(arr_sigma_c[mask_quality],bins=15,density=True,cumulative=True)
xvals = np.linspace(0,0.8,100)
x0=0.3
m=3.3
plt.plot(xvals, Weibull_cdf(xvals, x0=x0, m=m),label=f'Weibull, $\sigma_0$={x0}, m={m}')
x0=0.25
m=3
plt.plot(xvals, Weibull_cdf(xvals, x0=x0, m=m),label=f'Weibull, $\sigma_0$={x0}, m={m}')
x0=0.3
m=5
plt.plot(xvals, Weibull_cdf(xvals, x0=x0, m=m),label=f'Weibull, $\sigma_0$={x0}, m={m}')

plt.ylabel('Cumulative distribution function')
plt.xlabel('$\sigma_c$ [MPa]')
plt.legend()
plt.grid()