# -*- coding: utf-8 -*-
"""
Created on Sat May 30 10:39:06 2026

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import seaborn as sns
from scipy.stats import gaussian_kde

import icewave.sebastien.set_graphs as set_graphs

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')
#%% Load data 

year = '2025'
date = '0210' #date format, 'mmdd'
acqu_numb = '0003' #acquisition number 

path2data = os.path.join('F:/Rimouski_2025/Data/',date,'Geophones/')

# set fig_folder
fig_folder = f'{path2data}/MCMC_results/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)

file2load = f'{path2data}{year}_{date}_acq{acqu_numb}_MCMC_inversion_bidir_results.pkl'
with open(file2load,'rb') as pf:
    data = pickle.load(pf)
    
T_critic = data[0]
X = data[1]
misfit_accepted = data[2]


#%% Plot MCMC results as histograms

# Keep only the correct number of tests
nonzeroidx = np.nonzero(X[0, :])[0][-1]

data_to_plot = [X[i, :nonzeroidx] for i in range(4)]
dict_data = {'h':data_to_plot[0],'E':data_to_plot[1],
             'nu':data_to_plot[2],'rho':data_to_plot[3]}

set_graphs.set_matplotlib_param('square')
fig, axes = plt.subplots(nrows=2, ncols=2,constrained_layout = True)
axes = axes.flatten() # Transforme la matrice 2x2 en une liste 1D pour itérer facilement

# Dictionnaire pour stocker nos résultats mathématiques
results = {}

# 3. Boucle de traitement pour chaque variable
for i, (variable, valeurs) in enumerate(dict_data.items()):
    ax = axes[i]
    
    # a) Tracer l'histogramme de la densité de probabilité
    # density=True permet d'avoir une aire totale égale à 1 (probabilité)
    if i == 1:
        valeurs = valeurs/1e9
    elif i == 3:
        valeurs = valeurs/1017
        
    ax.hist(valeurs, bins=30, density=True, alpha=0.5, 
            color='skyblue', edgecolor='black', label='Histogramme')
    
    # b) Calculer le Kernel Density Estimate (KDE)
    kde = gaussian_kde(valeurs)
    x_eval = np.linspace(min(valeurs), max(valeurs), 1000) # Grille de 1000 points pour lisser la courbe
    kde_eval = kde(x_eval)
    
    # Tracer la courbe du KDE
    ax.plot(x_eval, kde_eval, color='darkblue', linewidth=2, label='KDE')
    
    # c) Déterminer la valeur ayant la plus forte probabilité (le pic du KDE)
    # On cherche l'index de la valeur maximale dans les résultats du KDE
    index_max_prob = np.argmax(kde_eval)
    valeur_optimale = x_eval[index_max_prob]
    
    # d) Déterminer un ordre de grandeur de l'incertitude (écart-type)
    ecart_type = np.std(valeurs)
    
    # Sauvegarde en mémoire
    results[variable] = {
        'max': valeur_optimale,
        'std': ecart_type
    }
    
    # e) Ajout d'indications visuelles sur le graphique
    # Ligne verticale pour la valeur optimale
    ax.axvline(valeur_optimale, color='red', linestyle='--', linewidth=2, 
               label=f'Valeur la plus probable : {valeur_optimale:.2f}')
    
    # Zone ombrée pour représenter l'incertitude (± 1 écart-type autour de la moyenne/mode)
    ax.axvspan(valeur_optimale - ecart_type, valeur_optimale + ecart_type, 
               color='orange', alpha=0.15, label=f'Incertitude (±{ecart_type:.2f})')
    
    # Habillage du graphique
    ax.set_ylabel('PDF')
    
axes[0].set_xlabel(r'$h_{ice} \; \mathrm{(m)}$')
axes[1].set_xlabel(r'$E \; \mathrm{(GPa)}$')
axes[2].set_xlabel(r'$\nu$')
axes[3].set_xlabel(r'$\rho_{ice}/\rho_w$')

figname = f'{fig_folder}{year}_{date}_acq{acqu_numb}_Histogramm_MCMC_results'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

file2save = f'{path2data}{year}_{date}_acq{acqu_numb}_MCMC_distribution_results.pkl'
with open(file2save,'wb') as pf:
    pickle.dump(results,pf)






