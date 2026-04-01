#%%
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import os
import pickle
# %%

disk = 'D:'

dates = ['1210','1211','1212','1215','1216','20260112','20260113','20260114','20260115','20260116','20260119','20260120','20260121','20260122','20260123']
# le 05/12 manip spéciale : changement de protocole pour fabrication échantillons :
# sortis de moules puis posés dans l'air froid du frigo, puis reposés dans l'eau à 0°
"""
temperature_samples = [0,0,0,0,0,-12,-12,-11,-11,-11,-1,-2,np.nan,np.nan,np.nan]
L_samples = np.array([8,8,8,8,8,8,8,8,8,8,8,8,8,8,8])*1e-2 # distance entre les 2 piliers en m (pas la taille réelle de l'échantillon)
w_samples = np.array([4,4,4,4,4,4,4,4,4,4,4,4,4,4,4])*1e-2 # largeur echantillon (m)
"""
temperature_samples = [0,0,0,0,0,-11,-11,-11,-1,-2,np.nan,np.nan,np.nan,np.nan,np.nan]
L_samples = np.array([8,8,8,8,8,8,8,8,8,8,8,8,8,8,8])*1e-2 # distance entre les 2 piliers en m (pas la taille réelle de l'échantillon)
w_samples = np.array([4,4,4,4,4,4,4,4,4,4,4,4,4,4,4])*1e-2 # largeur echantillon (m)

sigmac_dir_path = f'{disk}/manips_BsAs/Summary/sigma_c/'
mass_balance = 20 # masse de la balance de précision, en grammes
mass_pilar = 16 # masse du pilier, en grammes

dict_all_results = {}

thicknesses_avg = []
thicknesses_std = []
mc_avg = []
mc_err = []
array_L = []
array_temperatures_samples = []
array_defect_sizes = []

for i in range(len(dates)):
    date = dates[i]
    # Reading data from the pickle file
    with open(sigmac_dir_path+f'data_h_vs_mc_{date}.pkl', 'rb') as file:
        loaded_data = pickle.load(file)
    dict_all_results['dict_results_'+date] = loaded_data

    thicknesses_avg.append(dict_all_results['dict_results_'+date]['thicknesses_avg'])
    thicknesses_std.append(dict_all_results['dict_results_'+date]['thicknesses_std'])
    mc_avg.append(dict_all_results['dict_results_'+date]['mc_avg'])
    mc_err.append(dict_all_results['dict_results_'+date]['mc_err'])

    if np.isnan(L_samples[i]):
        array_L.append(dict_all_results['dict_results_'+date]['array_lengths'])
    else:
        array_L.append(np.ones(len(dict_all_results['dict_results_'+date]['thicknesses_avg'])) * L_samples[i])

    if 'defect_sizes' in dict_all_results['dict_results_'+date]:
        array_defect_sizes.append(dict_all_results['dict_results_'+date]['defect_sizes'])
    else:
        array_defect_sizes.append(np.ones(len(dict_all_results['dict_results_'+date]['thicknesses_avg'])) * np.nan)
    if np.isnan(temperature_samples[i]):
        array_temperatures_samples.append(dict_all_results['dict_results_'+date]['temperatures_samples'])
    else:
        array_temperatures_samples.append(temperature_samples[i] * np.ones(len(dict_all_results['dict_results_'+date]['thicknesses_avg'])))

thicknesses_avg = [item for sublist in thicknesses_avg for item in sublist]
thicknesses_std = [item for sublist in thicknesses_std for item in sublist]
mc_avg = [item for sublist in mc_avg for item in sublist]
mc_err = [item for sublist in mc_err for item in sublist]
array_L = [item for sublist in array_L for item in sublist]
array_temperatures_samples = [item for sublist in array_temperatures_samples for item in sublist]
array_defect_sizes = [item for sublist in array_defect_sizes for item in sublist]

thicknesses_avg = np.array(thicknesses_avg)
thicknesses_std = np.array(thicknesses_std)
mc_avg = np.array(mc_avg)
mc_err = np.array(mc_err)
array_L = np.array(array_L)
array_temperatures_samples = np.array(array_temperatures_samples)
array_defect_sizes = np.array(array_defect_sizes)

mask_nodefect = (np.isnan(array_defect_sizes)==True)



# %% plot results with labeled date

plt.figure()
for i in range(len(dates)):
    plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o',ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')

plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ (g)')
#²²plt.ylim(0,1500)
plt.xlim(0,15)
#plt.loglog()
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.show()

# %% plot results with labeled date : masse "effective"


plt.figure()
for i in range(len(dates)):
    dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o',ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')

A=3.5
plt.plot(np.linspace(3,14), A*np.linspace(3,14)**2,label='loi en h²')


plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ effective (g)')
#plt.ylim(100,4000)
#plt.xlim(1,8)
plt.loglog()
plt.ylim(10,6000)
plt.xlim(3,15)
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.show()

#%% mc eff vs h avec couleurs pour temperatures diff
%matplotlib qt
Lsample2plot = 8e-2

cmap = plt.cm.rainbow            # choix du colormap
norm = plt.Normalize(vmin=np.nanmin(array_temperatures_samples), vmax=np.nanmax(array_temperatures_samples))  # normalisation commune (peut être calculée dynamiquement)


plt.figure()

color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    plt.errorbar(thicknesses_avg[i], mc_avg[i] + mass_balance + mass_pilar, mc_err[i], thicknesses_std[i],linestyle='',marker='',ecolor='k')
    plt.scatter(thicknesses_avg[i], mc_avg[i] + mass_balance + mass_pilar, mc_err[i], color=color[i],norm=matplotlib.colors.LogNorm(),zorder=2)
plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ (g)')
plt.loglog()

plt.plot(np.linspace(3,16),50*np.linspace(3,16)**2)

plt.title('critical "force" (in grams) vs thickness (mm)')
#plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_granular_labeled_loglog.pdf',dpi=300)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

#%% mc eff vs h² avec couleurs pour temperatures diff
%matplotlib qt
Lsample2plot = 8e-2

cmap = plt.cm.rainbow            # choix du colormap
norm = plt.Normalize(vmin=-4, vmax=0)  # normalisation commune (peut être calculée dynamiquement)


plt.figure()

color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    plt.errorbar(thicknesses_avg[i]**2, mc_avg[i] + mass_balance + mass_pilar, mc_err[i], thicknesses_std[i]*thicknesses_avg[i]*2,linestyle='',marker='o', color=color[i],ecolor='k')

plt.xlabel('h² (mm²)')
plt.ylabel('$m_c$ effective (g)')
#plt.loglog()
plt.title('critical "force" (in grams) vs thickness (mm)')
#plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_granular_labeled_loglog.pdf',dpi=300)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

#%% contrainte vs h avec couleurs pour temperatures diff
%matplotlib qt
Lsample2plot = 8e-2

cmap = plt.cm.rainbow            # choix du colormap
norm = plt.Normalize(vmin=-10, vmax=0)  # normalisation commune (peut être calculée dynamiquement)


plt.figure()

color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    mceff = mc_avg[i] + mass_balance + mass_pilar
    w = w_samples[0] #pas de variation de w vs manips
    plt.errorbar(thicknesses_avg[i]*1e-3, 1.5*(mceff*1e-3*9.8*array_L[i])/(w*(thicknesses_avg[i]*1e-3)**2),  yerr=3 * ((mceff*1e-3*9.8 * array_L[i])/(w*((thicknesses_avg[i]*1e-3)**3)) * (thicknesses_std[i]*1e-3)), xerr=thicknesses_std[i]*1e-3, linestyle='',marker='o', color=color[i],ecolor='gray')

plt.xlabel('thickness (mm)')
plt.ylabel('sigma_c (Pa)')
#plt.loglog()
plt.title('critical stress (in grams) vs thickness (mm)')
#plt.legend()
#plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_granular_labeled_loglog.pdf',dpi=300)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

#%% contrainte vs T
%matplotlib qt
Lsample2plot = 8e-2

cmap = plt.cm.rainbow            # choix du colormap
norm = plt.Normalize(vmin=-4, vmax=0)  # normalisation commune (peut être calculée dynamiquement)


plt.figure()

color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    mceff = mc_avg[i] + mass_balance + mass_pilar
    w = w_samples[0] #pas de variation de w vs manips
    if (thicknesses_avg[i] > 0):
        plt.errorbar(array_temperatures_samples[i], 1.5*(mceff*1e-3*9.8*array_L[i])/(w*(thicknesses_avg[i]*1e-3)**2),  yerr=3 * ((mceff*1e-3*9.8 * array_L[i])/(w*((thicknesses_avg[i]*1e-3)**3)) * (thicknesses_std[i]*1e-3)), linestyle='',marker='o',color='b',ecolor='gray')

plt.xlabel('Temperature (°C)')
plt.ylabel('sigma_c (Pa)')

#plt.loglog()
plt.title('Contrainte à la rupture vs température')
#plt.legend()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.show()



# %%
