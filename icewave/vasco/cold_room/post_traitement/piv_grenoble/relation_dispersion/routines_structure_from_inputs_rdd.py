#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import csv
import re
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import copy  # pour éviter de modifier les données partagées
#%%
system_loc = 'windows_server'


##Windows
if system_loc=='windows_server':
    sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/')
    #sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/cold_room/post_traitement/piv_grenoble/relation_dispersion')
    
##Linux
elif system_loc=='linux_server':
    sys.path.append('/media/vasco/OS/Users/Vasco Zanchi/Documents/git_turbotice/')
from vasco.tools import open_csv
from vasco.cold_room.post_traitement.piv_grenoble.relation_dispersion.functions.disprel import *
# %% creation d'un dictionnaire pour ranger les résultats obtenus

##windows
if system_loc=='windows_server':
    inputs_file_path = 'R:/Gre25/Summary/dispersion_relations/inputs_RDD.csv'
##linux_server
elif system_loc=='linux_server':
    inputs_file_path = '/media/turbots/GreDisk/Gre25/Summary/dispersion_relations/inputs_RDD.csv'
## linux local
#inputs_file_path = '/media/vasco/Samsung_T5/Grenoble/Gre25/Summary/dispersion_relations/inputs_RDD.csv'


# inputs :
(dict_lists , list_categories , list_types) = open_csv.open_csv_table_experiments(file_path=inputs_file_path)



for i in range(len(list_categories)):
    if list_categories[i]=='\ufefflist_dates':
        list_categories[i] = 'list_dates'

print(list_categories)
print(list_types)

##windows
#general_folder = f'R:/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
##Linux
#general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'


if '\ufefflist_dates' in dict_lists:
    dict_lists['list_dates'] = dict_lists.pop('\ufefflist_dates')

dict_results = {}
for i in range(len(dict_lists['list_dates'])):
    dict_results[dict_lists['list_dates'][i]] = {}
for i in range(len(dict_lists['list_dates'])):
    dict_results[dict_lists['list_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])] = {}    
for i in range(len(dict_lists['list_dates'])):
    dict_results[dict_lists['list_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])]['camera_'+str(dict_lists['list_camera_SN'][i])] = {}


#%%
##############################################################################
####### Routine pour plein d'acquisitions, soit exécuté en parallele, ########
######  soit normalement #####################################################
##############################################################################

parallel=False

def process_single_case(i):
    date = dict_lists['list_dates'][i]
    acq = dict_lists['list_acq_num'][i]
    cam = dict_lists['list_camera_SN'][i]
    
    dd = copy.deepcopy(dict_results[date]['acq_' + str(acq)]['camera_' + str(cam)])
    dd['W'] = dict_lists['list_W'][i]
    dd['Dt'] = dict_lists['list_Dt'][i]

    ## Dossier général (Linux, sur serveur)
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date[-4:]}/cameras/manip_relation_dispersion/Acquisition_{acq}/camera_{cam}/'
    ##windows, local:
    #general_folder = f'D:/Grenoble/Gre25/Data/{date[-4:]}/cameras/manip_relation_dispersion/Acquisition_{acq}/camera_{cam}/'
    ##windows, serveur:
    #general_folder = f'R:/Gre25/Data/{date[-4:]}/cameras/manip_relation_dispersion/Acquisition_{acq}/camera_{cam}/'
    
    ##linux, local:
    #general_folder = f'/media/vasco/Samsung_T5/Grenoble/Gre25/Data/{date[-4:]}/cameras/manip_relation_dispersion/Acquisition_{acq}/camera_{cam}/'

    tab_f_exc, tab_freq_acq = list_all_freq(general_folder=general_folder)
    dd['tab_f_exc'] = tab_f_exc
    dd['tab_freq_acq'] = tab_freq_acq

    tab_v_phase = np.zeros(len(tab_f_exc))
    tab_v_phase_err = np.zeros(len(tab_f_exc))

    for j in range(len(tab_f_exc)):
        tab_v_phase[j], tab_v_phase_err[j] = routine_vitesse_phase(
            f_exc=tab_f_exc[j],
            freq_acq=tab_freq_acq[j],
            general_folder=general_folder,
            W=dd['W'],
            Dt=dd['Dt'],
            direction=-1,
            index_profile_line=dict_lists['list_indices_profile_line'][i],
            xlim_fit=(dict_lists['x_inf_fit'][i], dict_lists['x_sup_fit'][i]),
            dcm=dict_lists['list_dcm'][i],
            dpx=dict_lists['list_dpx'][i],
            plot=False
        )

    dd['tab_v_phase'] = tab_v_phase
    dd['tab_v_phase_err'] = tab_v_phase_err

    # Retourner l'entrée mise à jour avec sa clé d'origine
    return (date, acq, cam, dd)


if parallel: # dans le cas où on veut faire les calculs en parallel (beaucoup plus rapide)


    # Lancement du parallélisme
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_single_case, i) for i in range(len(dict_lists['list_dates']))]

        for future in futures:
            date, acq, cam, dd = future.result()
            dict_results[date]['acq_' + str(acq)]['camera_' + str(cam)] = dd

########### Dans le cas où on veut faire avec un seul coeur :

elif parallel==False:

    for i in range(len(dict_lists['list_dates'])):
        dd = dict_results[dict_lists['list_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])]['camera_'+str(dict_lists['list_camera_SN'][i])]
        dd['W'] = dict_lists['list_W'][i]
        dd['Dt'] = dict_lists['list_Dt'][i]
        #dd['tab_v_phase'] , dd['tab_v_phase_err'] = routine_vitesse_phase(f_exc=tab_f_exc[j],freq_acq=tab_freq_acq[j],general_folder=general_folder,W=dd['W'],Dt=dd['Dt'],index_profile_line=dict_lists['list_indices_profile_line'][i],xlim_fit=(dict_lists['list_xinf'][i],dict_lists['list_xsup'][i]),dcm=dict_lists['list_dcm'][i],dpx=dict_lists['list_dpx'])
        ##windows
        general_folder = f'R:/Gre25/Data/{dict_lists['list_dates'][i][-4:]}/cameras/manip_relation_dispersion/Acquisition_{str(dict_lists['list_acq_num'][i])}/camera_{dict_lists['list_camera_SN'][i]}/'
        ##Linux
        #general_folder = f'/media/turbots/GreDisk/Gre25/Data/{dict_lists['list_dates'][i][-4:]}/cameras/manip_relation_dispersion/Acquisition_{str(dict_lists['list_acq_num'][i])}/camera_{dict_lists['list_camera_SN'][i]}/'
        print(general_folder)

        tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)
        dd['tab_f_exc'] = tab_f_exc
        dd['tab_freq_acq'] = tab_freq_acq
        tab_v_phase = np.zeros(len(tab_f_exc))
        tab_v_phase_err = np.zeros(len(tab_f_exc)) 
        for j in range(len(tab_f_exc)):
            tab_v_phase[j],tab_v_phase_err[j] = routine_vitesse_phase(f_exc=tab_f_exc[j],freq_acq=tab_freq_acq[j],general_folder=general_folder,W=dd['W'],Dt=dd['Dt'],direction=-1,index_profile_line=dict_lists['list_indices_profile_line'][i],xlim_fit=(dict_lists['x_inf_fit'][i],dict_lists['x_sup_fit'][i]),dcm=dict_lists['list_dcm'][i],dpx=dict_lists['list_dpx'][i],plot=False)
        dd['tab_v_phase'] = tab_v_phase
        dd['tab_v_phase_err'] = tab_v_phase_err

        # enlever le try except (il faut avoir tous les fichiers piv prêts...)



#%% enregistrer les données


if system_loc=='linux_server':
    pkl_file_lists = "/media/turbots/GreDisk/Gre25/Summary/dispersion_relations/resultats/dict_lists.pkl"
elif system_loc=='windows_server':
    pkl_file_lists = "R:/Gre25/Summary/dispersion_relations/resultats/dict_lists.pkl"

if system_loc=='linux_server':
    pkl_file_results = "/media/turbots/GreDisk/Gre25/Summary/dispersion_relations/resultats/dict_results.pkl"
elif system_loc=='windows_server':
    pkl_file_results = "R:/Gre25/Summary/dispersion_relations/resultats/dict_results.pkl"



#pkl_file_lists = "/home/vasco/Bureau/donnees_temp/dict_lists.pkl"

# save pickle file
with open(pkl_file_lists, "wb") as f:
    dict_lists = pickle.dump(dict_lists,f)

# save pickle file
with open(pkl_file_results, "wb") as f:
    dict_lists = pickle.dump(dict_results,f)


# %%
