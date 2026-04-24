#%%
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import os
import pickle

matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
# %%

disk = 'D:'

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

# %% plot all results

mask = (thicknesses_std/thicknesses_avg)<1/15


plt.figure()
plt.errorbar(thicknesses_avg, mc_avg, yerr=mc_err, xerr=thicknesses_std,linestyle='',marker='o',color='b',ecolor='k')
plt.errorbar(thicknesses_avg[mask], mc_avg[mask], yerr=mc_err[mask], xerr=thicknesses_std[mask],linestyle='',marker='o',color='r',ecolor='k')

plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ (g)')
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.savefig(sigmac_dir_path + f'mc_vs_h_all_results.pdf',dpi=300)
plt.plot(np.linspace(1,10),50 * np.linspace(1,10)**2,'-')
#plt.xlim(0,)
plt.loglog()
plt.show()


# %% plot results with labeled date

plt.figure()
for i in range(len(dates)):
    plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o',ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')

plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ (g)')
plt.ylim(100,4000)
plt.xlim(1,8)
plt.loglog()
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.show()

# %% plot results with labeled date : masse "effective"


plt.figure()
for i in range(len(dates)):
    dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o',ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')

plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ effective (g)')
plt.ylim(100,6000)
plt.xlim(1,10)
plt.loglog()
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.show()
#%% même plot mais avec differentes couleurs suivant les temperatures

Lsample2plot = 8e-2
legendeachserie=False

plt.figure()
for i in range(len(dates)):
    dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    if temperature_samples[i] < 0:
        color = 'b'
    else:
        color = 'r'
    if L_samples[i]==Lsample2plot:
        if legendeachserie:
            plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o', color=color,ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C')
        else:
            plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o', color=color,ecolor='k')

if legendeachserie==False:
    plt.errorbar([],[],linestyle='',marker='o',color='b',ecolor='k',label='T$\simeq$-10°C')
    plt.errorbar([],[],linestyle='',marker='o',color='r',ecolor='k',label='T = 0°C')



plt.xlabel('$h$ [mm]', fontsize=15)
plt.ylabel('$m_c$ [g]', fontsize=15)
plt.ylim(100,8000)
plt.xlim(1,10)
plt.loglog()
plt.title('Critical force at fracture vs thickness', fontsize=15)
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.savefig('C:/Users/Vasco Zanchi/Desktop/presentations/interfreeze_2026/mc_eff_vs_h_all_results_labeled_loglog.pdf',dpi=300)

plt.show()

#%% même plot mais avec differentes couleurs suivant les tailles d'échantillons


plt.figure()
for i in range(len(dates)):
    dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'] = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    if temperature_samples[i] < 0:
        plot=False
    elif (L_samples[i]==8e-2)&(temperature_samples[i]==0):
        plot=True
        color='tab:orange'
    elif (L_samples[i]==12e-2)&(temperature_samples[i]==0):
        plot=True
        color='tab:green'
    if plot:
        plt.errorbar(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg'], dict_all_results['dict_results_'+dates[i]]['mc_avg_eff'], yerr=dict_all_results['dict_results_'+dates[i]]['mc_err'], xerr=dict_all_results['dict_results_'+dates[i]]['thicknesses_std'],linestyle='',marker='o', color=color,ecolor='k',label=dates[i]+' , T = '+str(temperature_samples[i])+' °C ; L='+str(L_samples[i]))

plt.xlabel('thickness (mm)')
plt.ylabel('$m_c$ effective (g)')
plt.ylim(100,4000)
plt.xlim(1,8)
plt.loglog()
plt.title('critical "force" (in grams) vs thickness (mm)')
plt.legend()
plt.savefig(sigmac_dir_path + f'mc_eff_vs_h_all_results_labeled_loglog.pdf',dpi=300)
plt.show()

#%% même plot mais avec differentes couleurs suivant les tailles d'échantillons

mask_temperature = (array_temperatures_samples==0)

plt.figure()
plt.errorbar(thicknesses_avg[mask_temperature&mask_nodefect], mc_avg[mask_temperature&mask_nodefect]+mass_balance+mass_pilar, yerr=mc_err[mask_temperature&mask_nodefect],xerr=thicknesses_std[mask_temperature&mask_nodefect],linestyle='',marker='', color=color,ecolor='k',zorder=1)
# Scatter
sc = plt.scatter(thicknesses_avg[mask_temperature&mask_nodefect], mc_avg[mask_temperature&mask_nodefect]+mass_balance+mass_pilar, c=array_L[mask_temperature&mask_nodefect], cmap='jet')

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
plt.loglog()
#plt.xlim(0.9,14)
#plt.ylim(0, 6000)
plt.ylabel('$m_c$ [g]', fontsize=15)
plt.xlabel('$h$ [mm]', fontsize=15)
plt.legend(handles=handles, title="Sample length [m]")
plt.title('Critical force at fracture vs thickness', fontsize=15)
plt.savefig('C:/Users/Vasco Zanchi/Desktop/presentations/interfreeze_2026/mc_vs_h_allsamples.pdf', dpi=300)

plt.show()

#%% contrainte vs h avec differentes couleurs suivant les tailles d'échantillons

w = 4e-2 

mask_temperature = (array_temperatures_samples==0)

#sigmac_array = (3/2) * (mc_avg*1e-3*9.81*array_L)/(w*(thicknesses_avg*1e-3)**2)
#sigmac_err_array = 3 * (mc_avg*1e-3*9.81*array_L) * (thicknesses_std*1e-3)/(w*(thicknesses_avg*1e-3)**3)
sigmac_array = (3/2) * ((mc_avg+mass_balance+mass_pilar)*1e-3*9.81*array_L)/(w*(thicknesses_avg*1e-3)**2)
sigmac_err_array = 3 * ((mc_avg+mass_balance+mass_pilar)*1e-3*9.81*array_L) * (thicknesses_std*1e-3)/(w*(thicknesses_avg*1e-3)**3)


plt.figure()
plt.errorbar(thicknesses_avg[mask_temperature&mask_nodefect], sigmac_array[mask_temperature&mask_nodefect],xerr=thicknesses_std[mask_temperature&mask_nodefect],yerr=sigmac_err_array[mask_temperature&mask_nodefect],linestyle='',marker='', color=color,ecolor='k',zorder=1)
# Scatter
sc = plt.scatter(thicknesses_avg[mask_temperature&mask_nodefect], sigmac_array[mask_temperature&mask_nodefect], c=array_L[mask_temperature&mask_nodefect], cmap='jet')

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
plt.xlim(0, np.nanmax(thicknesses_avg)*1.1)
plt.ylim(0, np.nanmax(sigmac_array)*1.1)
plt.xlabel('thickness (mm)')
plt.ylabel('$\sigma_c$ (Pa)')
plt.legend(handles=handles, title="Sample length [m]")
plt.savefig('C:/Users/Vasco Zanchi/Desktop/presentations/interfreeze_2026/sigmac_vs_h_variousL_allexp.pdf', dpi=300)
plt.show()

#%% contrainte vs h avec differentes couleurs suivant les tailles d'échantillons


thicknesses_avg_normalized = thicknesses_avg*1e-3/array_L
thicknesses_std_normalized = thicknesses_std*1e-3/array_L

plt.figure()
plt.errorbar(thicknesses_avg_normalized[mask_temperature&mask_nodefect], sigmac_array[mask_temperature&mask_nodefect],xerr=thicknesses_std_normalized[mask_temperature&mask_nodefect],yerr=sigmac_err_array[mask_temperature&mask_nodefect],linestyle='',marker='', color=color,ecolor='k',zorder=1)
# Scatter
sc = plt.scatter(thicknesses_avg_normalized[mask_temperature&mask_nodefect], sigmac_array[mask_temperature&mask_nodefect], c=array_L[mask_temperature&mask_nodefect], cmap='jet')

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
plt.xlim(0, np.nanmax(thicknesses_avg_normalized)*1.1)
plt.ylim(0, np.nanmax(sigmac_array)*1.1)
plt.xlabel('normalized thickness h/L')
plt.ylabel('$\sigma_c$ (Pa)')
plt.legend(handles=handles, title="Sample length [m]")
plt.savefig('C:/Users/Vasco Zanchi/Desktop/presentations/interfreeze_2026/sigmac_vs_hsurL_variousL_allexp.pdf', dpi=300)

plt.show()

#%%

#%% voir si la longueur de defaut "equivalente" dépend de l'épaisseur

# fracture toughness ice : between 50 and 150 kPa * m^1/2
K_Ic = 1.5 * 1e5

plt.figure()
#plt.errorbar(thicknesses_avg_normalized[mask_temperature&mask_nodefect], sigmac_array[mask_temperature&mask_nodefect],xerr=thicknesses_std_normalized[mask_temperature&mask_nodefect],yerr=sigmac_err_array[mask_temperature&mask_nodefect],linestyle='',marker='', color=color,ecolor='k',zorder=1)
# Scatter
sc = plt.scatter(thicknesses_avg_normalized[mask_temperature&mask_nodefect], ((1/(2*np.pi))*(K_Ic**2/sigmac_array[mask_temperature&mask_nodefect]**2))/array_L[mask_temperature&mask_nodefect], c=array_L[mask_temperature&mask_nodefect], cmap='jet')

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
plt.xlim(0, np.nanmax(thicknesses_avg_normalized)*1.1)
plt.ylim(0, np.nanmax(((1/(2*np.pi))*(K_Ic**2/sigmac_array[mask_temperature&mask_nodefect]**2))/array_L[mask_temperature&mask_nodefect])*0.3)
plt.plot(np.linspace(0, np.nanmax(thicknesses_avg_normalized)*1.1), np.linspace(0, np.nanmax(thicknesses_avg_normalized)*1.1), '--', color='k')
#plt.loglog()
plt.xlabel('normalized thickness h/L')
plt.ylabel('normalized equivalent defect length $a_{eq}/L$')
plt.legend(handles=handles, title="array_L")
plt.show()

# %% conversion en contrainte vs h/L

plt.figure()
for i in range(10,len(dates)):
    mc_avg_eff = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    force_eff_newtons = mc_avg_eff * 1e-3 * 9.81
    thicknesses_avg_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg']) * 1e-3
    thicknesses_std_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_std']) * 1e-3
    plt.errorbar(thicknesses_avg_meters/L_samples[i], (3/2) * (force_eff_newtons * L_samples[i])/(w_samples[i]*thicknesses_avg_meters**2), xerr=thicknesses_std_meters/L_samples[i], yerr=3 * ((force_eff_newtons * L_samples[i])/(w_samples[i]*thicknesses_avg_meters**3)) * thicknesses_std_meters,linestyle='',marker='o',ecolor='k',label=dates[i]+', L = '+str(L_samples[i])+' , T = '+str(temperature_samples[i])+' °C')
    plt.legend()

plt.xlabel('h/L')
plt.ylabel('$\sigma_c$ (Pa)')
plt.title('critical stress vs normalized thickness h/L')
plt.savefig(sigmac_dir_path + f'sigmac_vs_hsurL_all_results_labeled_loglog.pdf',dpi=300)
plt.show()

# %% conversion en contrainte vs h/L avec les temperatures

plt.figure()
for i in range(len(dates)):
    mc_avg_eff = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar
    force_eff_newtons = mc_avg_eff * 1e-3 * 9.81
    thicknesses_avg_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg']) * 1e-3
    thicknesses_std_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_std']) * 1e-3
    if temperature_samples[i] < 0:
        color = 'tab:blue'
    else:
        color = 'tab:orange'
    if L_samples[i]==0.08:
        marker='o'
        #color='tab:blue'
    else:
        marker='^'
        #color='tab:orange'
    plt.errorbar(thicknesses_avg_meters/L_samples[i], (3/2) * (force_eff_newtons * L_samples[i])/(w_samples[i]*thicknesses_avg_meters**2), xerr=thicknesses_std_meters/L_samples[i], yerr=3 * ((force_eff_newtons * L_samples[i])/(w_samples[i]*thicknesses_avg_meters**3)) * thicknesses_std_meters,marker=marker, color=color, linestyle='',ecolor='k',label=dates[i]+', L = '+str(L_samples[i])+' , T = '+str(temperature_samples[i])+' °C')

plt.xlabel('h/L')
plt.ylabel('$\sigma_c$ (Pa)')
plt.title('critical stress vs normalized thickness h/L')
#plt.legend()
plt.savefig(sigmac_dir_path + f'sigmac_vs_hsurL_all_results_labeled_loglog.pdf',dpi=300)
plt.show()
# %% plot de la contrainte critique qui prend en compte la taille des grains
# load data grain size

def load_grain_sizes_byhand(date,sample_name,disk=disk):
    path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{date}/{sample_name}_grains_sizes_mm.txt"
    datagrainsizes = np.loadtxt(path2grainsizes) 

    return datagrainsizes, np.mean(datagrainsizes), np.std(datagrainsizes)


cmap = plt.cm.viridis             # choix du colormap
norm = plt.Normalize(vmin=1, vmax=1.8)  # normalisation commune (peut être calculée dynamiquement)

plt.figure()
for i in range(len(dates)):
    print(i)
    for j in range(len(dict_all_results['dict_results_'+dates[i]]['list_samples'])): # boucle for sur toutes les acquisitions
        sample_name = dict_all_results['dict_results_'+dates[i]]['list_samples'][j][:3]
        path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{dates[i]}/{sample_name}_grains_sizes_mm.txt"
        #print(path2grainsizes)
        if os.path.exists(path2grainsizes):
            print('exists')
            datagrainsizes,grainsizeavg,_ = load_grain_sizes_byhand(dates[i], sample_name)
            color = cmap(norm(grainsizeavg))  # transformation valeur → couleur
            #color = cmap(norm(np.max(datagrainsizes)))
            #color = cmap(norm(np.min(datagrainsizes)))
            mc_avg_eff = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar # on veut afficher les points 1 par 1
            force_eff_newtons = mc_avg_eff * 1e-3 * 9.81
            thicknesses_avg_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg']) * 1e-3
            thicknesses_std_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_std']) * 1e-3
            plt.errorbar(thicknesses_avg_meters[j]/L_samples[i], (3/2) * (force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**2), xerr=thicknesses_std_meters[j]/L_samples[i], yerr=3 * ((force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**3)) * thicknesses_std_meters[j],linestyle='',marker='o',ecolor='k', color=color)
            #plt.errorbar(thicknesses_avg_meters[j]/L_samples[i], (3/2) * (force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**2), xerr=thicknesses_std_meters[j]/L_samples[i], yerr=3 * ((force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**3)) * thicknesses_std_meters[j],linestyle='',marker='o',ecolor='k', color=color)
        else:
            pass

# create a ScalarMappable so colorbar knows the mapping
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
plt.colorbar(sm, ax=plt.gca(), label='Average grain size (mm)')

plt.xlabel('h/L')
plt.ylabel('$\sigma_c$ (Pa)')
plt.title('critical stress vs normalized thickness h/L')
plt.legend()
#plt.savefig(sigmac_dir_path + f'sigmac_vs_hsurL_all_results_labeled_loglog.pdf',dpi=300)
plt.show()


# %% plot de la contrainte critique qui prend en compte la taille des grains
# load data grain size

cmap = plt.cm.viridis             # choix du colormap
norm = plt.Normalize(vmin=1, vmax=1.8)  # normalisation commune (peut être calculée dynamiquement)

plt.figure()
for i in range(len(dates)):
    print(i)
    for j in range(len(dict_all_results['dict_results_'+dates[i]]['list_samples'])): # boucle for sur toutes les acquisitions
        sample_name = dict_all_results['dict_results_'+dates[i]]['list_samples'][j][:3]
        path2grainsizes = f"{disk}/manips_BsAs/Summary/polariseurs/{dates[i]}/{sample_name}_grains_sizes_mm.txt"
        print(path2grainsizes)
        if os.path.exists(path2grainsizes):
            print(dates[i],sample_name)
            print('exists')
            datagrainsizes,grainsizeavg,grainsizestd = load_grain_sizes_byhand(dates[i], sample_name)
            color = cmap(norm(grainsizeavg))  # transformation valeur → couleur
            #color = cmap(norm(np.max(datagrainsizes)))
            #color = cmap(norm(np.min(datagrainsizes)))
            mc_avg_eff = dict_all_results['dict_results_'+dates[i]]['mc_avg'] + mass_balance + mass_pilar # on veut afficher les points 1 par 1
            force_eff_newtons = mc_avg_eff * 1e-3 * 9.81
            thicknesses_avg_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_avg']) * 1e-3
            thicknesses_std_meters = np.array(dict_all_results['dict_results_'+dates[i]]['thicknesses_std']) * 1e-3
            #plt.errorbar(thicknesses_avg_meters[j]/L_samples[i], (3/2) * (force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**2), xerr=thicknesses_std_meters[j]/L_samples[i], yerr=3 * ((force_eff_newtons[j] * L_samples[i])/(w_samples[i]*thicknesses_avg_meters[j]**3)) * thicknesses_std_meters[j],linestyle='',marker='o',ecolor='k', color=color)
            if np.isnan(L_samples[i]): 
                L_sample = dict_all_results['dict_results_'+dates[i]]['array_lengths'][j]
            else:
                L_sample = L_samples[i]
            plt.errorbar(grainsizeavg, (3/2) * (force_eff_newtons[j] * L_sample)/(w_samples[i]*thicknesses_avg_meters[j]**2), xerr=grainsizestd, yerr=3 * ((force_eff_newtons[j] * L_sample)/(w_samples[i]*thicknesses_avg_meters[j]**3)) * thicknesses_std_meters[j],linestyle='',marker='o',ecolor='k',color='b')

        else:
            pass

# create a ScalarMappable so colorbar knows the mapping
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])   # nécessaire pour certaines versions de matplotlib
#plt.colorbar(sm, ax=plt.gca(), label='Average grain size (mm)')

plt.xlabel('grain size (mm)')
plt.ylabel('$\sigma_c$ (Pa)')
plt.title('critical stress vs grain size')
plt.legend()
#plt.xlim(0, 10)
#plt.ylim(0,4e6)
plt.loglog()
plt.show()

# %%


dict_path = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/dict_all_results_1darrays.pkl"

dict2save = {}

dict2save['thicknesses_avg'] = thicknesses_avg
dict2save['thicknesses_std'] = thicknesses_std
dict2save['mc_avg'] = mc_avg
dict2save['mc_err'] = mc_avg
dict2save['mass_balance'] = mass_balance
dict2save['mass_pilar'] = mass_pilar
dict2save['sigmac_array'] = sigmac_array
dict2save['sigmac_err_array'] = sigmac_err_array
dict2save['array_L'] = array_L
dict2save['array_temperatures_samples'] = array_temperatures_samples
dict2save['array_defect_sizes'] = array_defect_sizes
dict2save['mask_nodefect'] = mask_nodefect

with open(dict_path, 'wb') as f:
    pickle.dump(dict2save, f)
