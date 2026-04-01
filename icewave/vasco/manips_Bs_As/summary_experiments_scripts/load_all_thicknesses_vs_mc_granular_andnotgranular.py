#%%
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import os
import pickle
# %%

disk = 'D:'

def plot_data_without_grains(Lsample2plot=8e-2):
    dates = ['1002', '1003','1006','1007','1009','1013','1014','1016','1020','1022','1023','1024','1110','1111','1113','1118','1119','1120','1125','1126','1127','1128','1202','1203','1205','20260126','20260127']
    temperature_samples = [0, 0, -15, -15, -13, 0, 0, 0,-15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-11,-10]
    L_samples = np.array([8,8,8,8,8,8,8,12,8,8,12,12,12,12,12,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,8,8])*1e-2 # distance entre les 2 piliers en m (pas la taille réelle de l'échantillon)
    w_samples = np.array([4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4])*1e-2 # largeur echantillon (m)

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

    # même plot mais avec differentes couleurs suivant les tailles d'échantillons

    mask_temperature = (array_temperatures_samples==0)
    mask_sizes = (array_L==0.08)
    
    plt.errorbar(thicknesses_avg[mask_temperature&mask_nodefect&mask_sizes], mc_avg[mask_temperature&mask_nodefect&mask_sizes]+mass_balance+mass_pilar, yerr=mc_err[mask_temperature&mask_nodefect&mask_sizes],xerr=thicknesses_std[mask_temperature&mask_nodefect&mask_sizes],linestyle='',marker='',ecolor='k',zorder=1)
    # Scatter
    sc = plt.scatter(thicknesses_avg[mask_temperature&mask_nodefect&mask_sizes], mc_avg[mask_temperature&mask_nodefect&mask_sizes]+mass_balance+mass_pilar, c=array_L[mask_temperature&mask_nodefect&mask_sizes], cmap='jet')

    # Valeurs uniques
    vals = np.unique(array_L)

    # Récupérer le cmap et normalisation
    cmap = sc.cmap
    norm = sc.norm

    # Créer une entrée de légende par valeur unique
    handles = [
        plt.Line2D(
            [], [], 
            marker="o", linestyle="", 
            color=cmap(norm(v)), 
            label=str(v)
        )
        for v in vals
    ]
    #plt.colorbar()

def plot_data_without_grains_temp(Lsample2plot=8e-2, plot=False):
    dates = ['1002', '1003','1006','1007','1009','1013','1014','1016','1020','1022','1023','1024','1110','1111','1113','1118','1119','1120','1125','1126','1127','1128','1202','1203','1205','20260126','20260127']
    temperature_samples = [0, 0, -15, -15, -13, 0, 0, 0,-15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-11,-10]
    L_samples = np.array([8,8,8,8,8,8,8,12,8,8,12,12,12,12,12,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,8,8])*1e-2 # distance entre les 2 piliers en m (pas la taille réelle de l'échantillon)
    w_samples = np.array([4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4])*1e-2 # largeur echantillon (m)

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

    if plot:
        for i in range(len(dates)):
            dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
            if temperature_samples[i] < 0:
                color = 'b'
            else:
                color = 'r'
            if L_samples[i]==Lsample2plot:
                plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o', color=color,ecolor='gray',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')
    return thicknesses_avg, thicknesses_std, array_temperatures_samples, mc_avg, mc_err

def plot_data_granular(Lsample2plot=8e-2):
    dates = ['1210','1211','1212','1215','1216','1218','1219','20260112','20260113','20260114','20260115','20260116','20260119','20260120','20260121','20260122','20260123']
    # le 05/12 manip spéciale : changement de protocole pour fabrication échantillons :
    # sortis de moules puis posés dans l'air froid du frigo, puis reposés dans l'eau à 0°
    temperature_samples = [0,0,0,0,0,-12,-12,-11,-11,-11,-1,-2,np.nan,np.nan,np.nan,np.nan,np.nan]
    L_samples = np.array([8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8])*1e-2 # distance entre les 2 piliers en m (pas la taille réelle de l'échantillon)
    w_samples = np.array([4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4])*1e-2 # largeur echantillon (m)

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

    # même plot mais avec differentes couleurs suivant les temperatures

    #Lsample2plot = 8e-2
    """
    for i in range(len(dates)):
        dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
        if temperature_samples[i] < -8:
            color = 'tab:blue'
        elif (temperature_samples[i]<0)&(temperature_samples[i]>-4):
            color='purple'
        else:
            color = 'tab:orange'

        if L_samples[i]==Lsample2plot:
            plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='^', color=color,ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')
    """
    return thicknesses_avg, thicknesses_std, array_temperatures_samples, mc_avg, mc_err

#%%
cmap = plt.cm.viridis            # choix du colormap
norm = plt.Normalize(vmin=-10, vmax=0)  # normalisation commune (peut être calculée dynamiquement)

plt.figure()
#plot_data_without_grains()
plot_data_without_grains_temp(plot=True)
#plot_data_granular()
mass_balance=20
mass_pilar=16
thicknesses_avg, thicknesses_std, array_temperatures_samples, mc_avg, mc_err = plot_data_granular()
color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    plt.errorbar(thicknesses_avg[i], mc_avg[i] + mass_balance + mass_pilar, mc_err[i], thicknesses_std[i],linestyle='',marker='^', color=color[i],ecolor='gray')

plt.loglog()
plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ (g)')
plt.title('critical "force" (in grams) vs thickness (mm)')
a=50
plt.plot(np.linspace(2,20), a*np.linspace(2,20)**2)
a=30
plt.plot(np.linspace(2,20), a*np.linspace(2,20)**2)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

#%%
cmap = plt.cm.jet            # choix du colormap
cmap_wg = plt.cm.jet
norm = plt.Normalize(vmin=-10, vmax=0)  # normalisation commune (peut être calculée dynamiquement)

plt.figure()
plt.style.use("seaborn-v0_8")
# wg for "without grains" :
thicknesses_avg_wg, thicknesses_std_wg, array_temperatures_samples_wg, mc_avg_wg, mc_err_wg = plot_data_without_grains_temp()
#plot_data_granular()
mass_balance=20
mass_pilar=16
thicknesses_avg, thicknesses_std, array_temperatures_samples, mc_avg, mc_err = plot_data_granular()
color_wg = cmap_wg(norm(array_temperatures_samples_wg))
color = cmap(norm(array_temperatures_samples))
for i in range(len(thicknesses_avg)):
    plt.errorbar(thicknesses_avg[i]**2, mc_avg[i] + mass_balance + mass_pilar, mc_err[i], thicknesses_std[i] * thicknesses_avg[i] * 2,linestyle='',marker='^', color=color[i],ecolor='gray')
for i in range(len(thicknesses_avg_wg)):
    plt.errorbar(thicknesses_avg_wg[i]**2, mc_avg_wg[i] + mass_balance + mass_pilar, mc_err_wg[i], thicknesses_std_wg[i] * thicknesses_avg_wg[i] * 2,linestyle='',marker='o', color=color_wg[i],ecolor='gray',alpha=0.4)


#plt.loglog()
plt.xlabel('h² (mm²)')
plt.ylabel('$m_c$ effective (g)')
plt.title('critical "force" (in grams) vs thickness (mm)')
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

# %%
cmap = plt.cm.jet            # choix du colormap
cmap_wg = plt.cm.jet
norm = plt.Normalize(vmin=-10, vmax=0)  # normalisation commune (peut être calculée dynamiquement)

plt.figure()
# wg for "without grains" :
thicknesses_avg_wg, thicknesses_std_wg, array_temperatures_samples_wg, mc_avg_wg, mc_err_wg = plot_data_without_grains_temp()
#plot_data_granular()
mass_balance=20
mass_pilar=16
thicknesses_avg, thicknesses_std, array_temperatures_samples, mc_avg, mc_err = plot_data_granular()
color_wg = cmap_wg(norm(array_temperatures_samples_wg))
color = cmap(norm(array_temperatures_samples))

L = 0.08
w = 0.04


for i in range(len(thicknesses_avg_wg)):
    Fc_wg = (mc_avg_wg[i] + mass_balance + mass_pilar)*9.81 * 1e-3
    plt.errorbar(thicknesses_avg_wg[i]*1e-3, (3/2)*(Fc_wg*L)/(w*(thicknesses_avg_wg[i]*1e-3)**2),yerr= 3 * ((mc_avg_wg[i]*1e-3*9.8 * L)/(w*((thicknesses_avg_wg[i]*1e-3)**3)) * (thicknesses_std_wg[i]*1e-3)),xerr= thicknesses_std_wg[i]*1e-3,linestyle='',marker='o', color=color_wg[i],ecolor='gray',alpha=0.6)

for i in range(len(thicknesses_avg)):
    Fc = (mc_avg[i] + mass_balance + mass_pilar)*9.81 * 1e-3
    plt.errorbar(thicknesses_avg[i]*1e-3, (3/2)*(Fc*L)/(w*(thicknesses_avg[i]*1e-3)**2), yerr=3 * ((mc_avg[i]*1e-3*9.8 * L)/(w*((thicknesses_avg[i]*1e-3)**3)) * (thicknesses_std[i]*1e-3)),xerr= thicknesses_std[i]*1e-3,linestyle='',marker='^', color=color[i],ecolor='gray')

#plt.loglog()
plt.xlabel('h (m)')
plt.ylabel('$\sigma_c$ (Pa)')
plt.title('Critical stress (Pa) vs thickness (m)')
#plt.xlim(0,1.8e-2)
#plt.ylim(0,5e6)
#plt.yscale('log')
plt.loglog()
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Temperature of sample (°C)')
plt.show()

# %%
