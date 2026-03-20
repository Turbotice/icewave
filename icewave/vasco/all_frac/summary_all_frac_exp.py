#%%
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"

#%%

disk_Gre = 'D:/Grenoble'
path2dict_results_fracwaves_Gre = f'{disk_Gre}/Gre25/Summary/fracture_postprocessing/resultats/dict_fracwaves_Gre.pkl'

with open(path2dict_results_fracwaves_Gre, 'rb') as f:
    dict_results_fracwaves_Gre = pickle.load(f)


disk = 'D:'
dict_dir = f"{disk}/manips_BsAs/Summary/dictionaries_alldata/"

path2data_BsAs_Gre25 = f"{dict_dir}dict_BsAs_Gre25_data.pkl"

with open(path2data_BsAs_Gre25, 'rb') as f:
    dict_BsAs_Gre25 = pickle.load(f)


dir2savefig = f'{disk}/General/figures'

# %%
dict_results_fracwaves_Gre

# %%
dict_BsAs_Gre25
#%%
# premierement on peut afficher Gre25 et BsAs sur 
# le même graphe pour epsilonc vs h par exemple
#####################################################
exp_success_Gre = dict_results_fracwaves_Gre['tab_plot']
maskerr = dict_results_fracwaves_Gre['tab_epsilonc_err']/dict_results_fracwaves_Gre['tab_epsilonc_avg'] < 1/5

maskplotGre = maskerr&exp_success_Gre
indices2plotGre = np.where(maskplotGre)[0]
#####################################################
mask_samples8cm = dict_BsAs_Gre25['L']==0.08
mask_samples12cm = dict_BsAs_Gre25['L']==0.12

mask_samples_size = mask_samples8cm | mask_samples12cm
#mask_temperature = dict_BsAs_Gre25['temperatures_samples']==0
mask_positive = dict_BsAs_Gre25['epsilonc_avg']>0
mask_errBsAs = dict_BsAs_Gre25['epsilonc_err']/dict_BsAs_Gre25['epsilonc_avg']<1/5
indices2plotBsAs = np.where(mask_samples_size&mask_positive&mask_errBsAs)[0]

colors = np.empty(len(mask_samples_size), dtype=object)
for i in range(len(colors)):
    if dict_BsAs_Gre25['temperatures_samples'][i]==0:
        colors[i] = 'red'
    elif dict_BsAs_Gre25['temperatures_samples'][i]<0:
        colors[i] = 'blue'
    else:
        print("positive temperatures not possible!")

symbols = np.empty(len(mask_samples_size), dtype=object)
for i in range(len(symbols)):
    if dict_BsAs_Gre25['L'][i]==0.08:
        symbols[i] = 's'
    elif dict_BsAs_Gre25['L'][i]==0.12:
        symbols[i] = '^'
    else:
        symbols[i] = 'o'
        print("different size than 8cm or 12cm found !")

#####################################################

plt.figure()
plt.errorbar(dict_results_fracwaves_Gre['tab_h_avg'][indices2plotGre],
             dict_results_fracwaves_Gre['tab_epsilonc_avg'][indices2plotGre],
             yerr=dict_results_fracwaves_Gre['tab_epsilonc_err'][indices2plotGre],
             xerr=dict_results_fracwaves_Gre['tab_h_std'][indices2plotGre],
             linestyle='',ecolor='gray')
plt.scatter(dict_results_fracwaves_Gre['tab_h_avg'][indices2plotGre],
            dict_results_fracwaves_Gre['tab_epsilonc_avg'][indices2plotGre],
            c='green', zorder=3, label='Wave tank, Grenoble')

plt.errorbar(dict_BsAs_Gre25['h_avg'][indices2plotBsAs],
              dict_BsAs_Gre25['epsilonc_avg'][indices2plotBsAs],
                yerr=dict_BsAs_Gre25['epsilonc_err'][indices2plotBsAs],
                  xerr=dict_BsAs_Gre25['h_std'][indices2plotBsAs],
                    linestyle='',ecolor='gray')
for i in range(len(colors)):
    if i in indices2plotBsAs:
        plt.scatter(dict_BsAs_Gre25['h_avg'][i],
                    dict_BsAs_Gre25['epsilonc_avg'][i],zorder=2,
                    edgecolors=colors[i],c='white', marker=symbols[i])

plt.scatter([], [], edgecolors='red', facecolors='white',
            marker='s', label='BsAs, L=8cm, T=0°C')
plt.scatter([], [], edgecolors='blue', facecolors='white',
            marker='s', label='BsAs, L=8cm, T=-10°C')
plt.scatter([], [], edgecolors='red', facecolors='white',
            marker='^', label='BsAs, L=12cm, T=0°C')

plt.xlabel('Thickness (m)', fontsize=13)
plt.ylabel('$\epsilon_c$', fontsize=13)
plt.legend()
plt.loglog()
plt.title("$\epsilon_c$ vs $h$, for wave tank and flexural 3 points tests", fontsize=15)
plt.savefig(f'{dir2savefig}/epsilonc_vs_hGre25_BsAs_loglog.pdf', dpi=300)

#plt.loglog()
# %%
mask_Gre25samples = dict_BsAs_Gre25['Gre25_samples']['Critical_stress_MPa']>1.5
indices2plotGre25samples = np.where(mask_Gre25samples)[0]

E_eff_moy = 5e9
E_eff_err = 2e9

plt.figure()
plt.title("Data from Grenoble (wave tank and 3pts tests)")
plt.errorbar(dict_BsAs_Gre25['Gre25_samples']['Thickness'][indices2plotGre25samples]*1e-3, 
         dict_BsAs_Gre25['Gre25_samples']['Critical_stress_MPa'][indices2plotGre25samples]*1e6/E_eff_moy,
         yerr=dict_BsAs_Gre25['Gre25_samples']['Critical_stress_MPa'][indices2plotGre25samples]*1e6*E_eff_err/(E_eff_moy**2),
         linestyle='',
         marker='o',
             label='3pts bending tests from Gre25, T=-10°C')

plt.errorbar(dict_results_fracwaves_Gre['tab_h_avg'][indices2plotGre],
             dict_results_fracwaves_Gre['tab_epsilonc_avg'][indices2plotGre],
             yerr=dict_results_fracwaves_Gre['tab_epsilonc_err'][indices2plotGre],
             xerr=dict_results_fracwaves_Gre['tab_h_std'][indices2plotGre],
             linestyle='',ecolor='gray')
plt.scatter(dict_results_fracwaves_Gre['tab_h_avg'][indices2plotGre],
            dict_results_fracwaves_Gre['tab_epsilonc_avg'][indices2plotGre],
            c='green', zorder=3, label='Wave tank, Grenoble')


plt.xlabel('Thickness (m)', fontsize=13)
plt.ylabel('$\epsilon_c$', fontsize=13)
plt.legend()
plt.loglog()
plt.savefig(f'{dir2savefig}/epsilonc_vs_h_Gre25_wavetank_3pts_loglog.pdf', dpi=300)

# %%

plt.errorbar(dict_BsAs_Gre25['h_avg'][indices2plotBsAs],
              dict_BsAs_Gre25['sigma_c_avg'][indices2plotBsAs],
                yerr=dict_BsAs_Gre25['sigma_c_std'][indices2plotBsAs],
                  xerr=dict_BsAs_Gre25['h_std'][indices2plotBsAs],
                    linestyle='',ecolor='gray')
for i in range(len(colors)):
    if i in indices2plotBsAs:
        plt.scatter(dict_BsAs_Gre25['h_avg'][i],
                    dict_BsAs_Gre25['sigma_c_avg'][i],zorder=2,
                    edgecolors=colors[i],c='white', marker=symbols[i])