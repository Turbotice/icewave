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

from functions.disprel import *

#matplotlib.use('TkAgg')
#%% changer generalfolder si besoin

#date = '20241203'
date = '0528'

acq_num = 4
#camera_SN = '40300722'
#camera_SN = '40437120'
camera_SN = 40307970

#W = 64
#Dt = 20
W = 64
Dt = 50
direction = -1

system_loc = 'windows_server'

if system_loc=='windows_local':
    general_folder = f'D:/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
elif system_loc=='linux_local':
    general_folder = f'/media/vasco/Samsung_T5/Grenoble/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
elif system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'



#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_1/camera_{camera_SN}/'

# pour la manip plexi :
#camera_SN = '40307970'
#camera_SN = '40300722'
#date = '20250130'
#acq_num = 5
#general_folder = f'/run/user/1003/gvfs/afp-volume:host=thiou.local,user=vasco,volume=labshared1/Banquise/Vasco/plexi/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#W = 64
#Dt = 1
#direction = 1

#%% lister toutes les frequences
tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)

#%% lancer routine
tab_v_phase = np.zeros(len(tab_f_exc))
tab_v_phase_err = np.zeros(len(tab_f_exc))
for i in range(len(tab_f_exc)):
    tab_v_phase[i],tab_v_phase_err[i] = routine_vitesse_phase(f_exc=tab_f_exc[i],freq_acq=tab_freq_acq[i],general_folder=general_folder,W=W,Dt=Dt,direction=direction,index_profile_line=3,xlim_fit=(0,50),dcm=16,dpx=1452)#,camera_SN='40437120')#'40300722')#)

#%%  relation dispersion
def v_phase_flexural(f,D):
    rho=1e3
    omega = 2*np.pi*f
    return (D/rho)**(1/5) * omega**(3/5)


E = 1.5e9
h = 2.4e-3
nu = 0.4
D = (E*h**3)/(12*(1-nu**2))

mask = (tab_v_phase_err < 1/30*tab_v_phase) & (tab_v_phase>0)# & (tab_f_exc<100)


popt,pcov = curve_fit(v_phase_flexural, tab_f_exc[mask],tab_v_phase[mask])

plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='.',linestyle='')
plt.plot(tab_f_exc,v_phase_flexural(tab_f_exc,D))
plt.plot(tab_f_exc,v_phase_flexural(tab_f_exc,popt[0]),label='fit : D='+str(np.round(popt[0],1))+' +- '+str(np.round(np.sqrt(pcov[0][0])))+' J')
plt.fill_between(tab_f_exc,v_phase_flexural(tab_f_exc,popt[0]-np.sqrt(pcov[0][0])),v_phase_flexural(tab_f_exc,popt[0]+np.sqrt(pcov[0][0])),color='red',alpha=0.2)

#plt.ylim(0,60)
#plt.xlim(0,200)
plt.xlabel('frequency (Hz)',fontsize=15)
plt.ylabel('$v_{\phi} (m/s)$',fontsize=15)
plt.legend()
plt.show()



# %% creation d'un dictionnaire pour ranger les résultats obtenus

# inputs :

list_categories,dict_lists,list_types = load_csv_input_RDD(file_path='R:/Gre25/Summary/dispersion_relations/inputs_RDD.csv')# le file path est par defaut celui dans hublot

def assign_input_param2global_variables(dict_lists=dict):
    """
    cree les variables dont les noms et les valeurs sont stockées dans dict_lists en tant que variables globales
    """
    for key in dict_lists:
        globals()[key] = dict_lists[key]

assign_input_param2global_variables(dict_lists=dict_lists)
list_xlim_fit = []
for i in range(len(x_inf_fit)):
    list_xlim_fit.append((x_inf_fit[i],x_sup_fit[i]))
# function to run and store the routine for a chosen acquisition and desired values of params
def routine_et_structure(idx_routine):
    date = list_dates[idx_routine]
    acq_num = list_acq_num[idx_routine]
    W = list_W[idx_routine]
    Dt = list_Dt[idx_routine]
    camera_SN = list_camera_SN[idx_routine]
    index_profile_line = list_indices_profile_line[idx_routine]
    xlim_fit = list_xlim_fit[idx_routine]
    dcm = list_dcm[idx_routine]
    dpx = list_dpx[idx_routine]
    

    if (date in dico) == False:
        dico[date] = {}
    if ('acq'+str(acq_num) in dico[date])==False:
        dico[date]['acq'+str(acq_num)] = {}
    dico_acq = dico[date]['acq'+str(acq_num)]
    if ('camera_'+camera_SN in dico_acq)==False:
        dico_acq['camera_'+camera_SN] = {}
    dico_camera = dico_acq['camera_'+camera_SN]
    if (f'W{str(W)}_Dt{str(Dt)}' in dico_camera)==False:
        dico_camera[f'W{str(W)}_Dt{str(Dt)}'] = {}
    dico_piv = dico_camera[f'W{str(W)}_Dt{str(Dt)}']
    if (f'index_profile_line{str(index_profile_line)}' in dico_piv)==False:
        dico_piv[f'index_profile_line{str(index_profile_line)}'] = {}
    dico_piv_profile = dico_piv[f'index_profile_line{str(index_profile_line)}']
    if (f'_xlim_fit{str(xlim_fit)}' in dico_piv_profile)==False:
        dico_piv_profile[f'_xlim_fit{str(xlim_fit)}'] = {}
    
    dico_pivprof_xlimfit = dico_piv_profile[f'_xlim_fit{str(xlim_fit)}']


    # run the routine if data is not in the dictionary:
    if ('tab_v_phase' in dico_pivprof_xlimfit)==False:
        dico_pivprof_xlimfit['dcm'] = list_dcm[idx_routine]
        dico_pivprof_xlimfit['dpx'] = list_dpx[idx_routine]
        if computer=='DellVasco':
            general_folder = f'K:/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
        elif computer=='Leyre':
            general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

        tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)
        dico_pivprof_xlimfit['tab_f_exc'] = tab_f_exc
        dico_pivprof_xlimfit['tab_freq_acq'] = tab_freq_acq
        tab_v_phase = np.zeros(len(tab_f_exc))
        tab_v_phase_err = np.zeros(len(tab_f_exc))
        for i in range(len(tab_f_exc)):
            tab_v_phase[i],tab_v_phase_err[i] = routine_vitesse_phase(f_exc=tab_f_exc[i],freq_acq=tab_freq_acq[i],general_folder=general_folder,W=W,Dt=Dt,index_profile_line=index_profile_line,xlim_fit=xlim_fit,dcm=dcm,dpx=dpx,plot=False)#,camera_SN='40437120')#'40300722')#)

        dico_pivprof_xlimfit['tab_v_phase'] = tab_v_phase
        dico_pivprof_xlimfit['tab_v_phase_err'] = tab_v_phase_err

#############################################################
# run for several acqusitions
## first, if required, create an empty dict
if ('dico' in globals())==False:
    dico = {}
## then fill it with the results
for idx_routine in range(len(list_dates)):
    routine_et_structure(idx_routine)

################################################################
# FIN DE LA PARTIE OÙ ON RUN POUR PLEIN D'ACQUISITION ##########
################################################################
#%%
plt.figure(figsize=(15,10))
for idx_routine in range(len(list_dates)):
    tab_f_exc = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_f_exc']
    tab_freq_acq = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_freq_acq']
    tab_v_phase = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase']
    tab_v_phase_err = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase_err']
    #print(tab_v_phase)
    mask = (tab_v_phase_err < 1/20*tab_v_phase) & (tab_v_phase>0)
    plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='o',linestyle='',label=list_dates[idx_routine])
    #plt.errorbar(tab_f_exc,tab_v_phase,tab_v_phase_err,marker='.',linestyle='')
plt.ylim(1,40)
plt.xlabel('frequency (Hz)')
plt.ylabel('phase velocity (m/s)')
plt.legend()
#plt.loglog()
plt.show()
#%% fitter en regime flexion pour avoir modul de flexion
def compute_omega(k,D):
    return np.sqrt(((D/rho)*k**5 + (T/rho)*k**3 + g*k)*np.tanh(k*H))
def compute_freq(k,D):
    return (1/(2*np.pi))*compute_omega(k,D)





"""
e = 2e-3 
H = 14.5e-2 - e 
E = 9e9
nu = 0.4
D_value = (E*(e**3)/(12*(1-nu**2)))
rho = 1e3
T = 0
lambdamin = 5e-2
lambdamax = 1
g = 9.81
tab_lambda = np.linspace(lambdamin,lambdamax,10000)
"""
########################## ICI CHOISIR "L'INDICE ROUTINE" A AFFICHER ET FITTER ##############
list_idx_routine=[0,1,15,16]
#############################################################################################

plt.figure(figsize=(15,10))

for idx_routine in list_idx_routine:
    tab_f_exc = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_f_exc']
    tab_freq_acq = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_freq_acq']
    tab_v_phase = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase']
    tab_v_phase_err = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase_err']


    mask = (tab_v_phase_err < 1/10*tab_v_phase) & (tab_v_phase>0)
    plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='o',linestyle='',label=list_dates[idx_routine]+' ; acq'+str(list_acq_num[idx_routine])+' ; cam '+str(list_camera_SN[idx_routine]))


    popt,pcov = curve_fit(v_phase_flexural,tab_f_exc[mask],tab_v_phase[mask],sigma=tab_v_phase_err[mask])
    #plt.figure()
    tab_f_fit = np.linspace(10,200,1000)

    plt.plot(tab_f_fit,v_phase_flexural(tab_f_fit,popt[0]),label='Fit (if only flexural) : D = ('+str(np.round(popt[0],decimals=2))+ '+-'+ str(np.round(np.sqrt(pcov[0][0]),decimals=2)) + ') J')

    #plt.errorbar(np.array(dico_dataset['tab_f_exc'])[mask],np.array(dico_dataset['v_phase'])[mask],np.array(dico_dataset['v_phase_err'])[mask],marker='o')
    #plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='.',linestyle='')
    plt.ylim(0,50)
    plt.xlim(0,300)
    plt.xlabel('frequency (Hz)',fontsize=15)
    plt.ylabel('$v_{\phi}$',fontsize=15)
plt.legend()
plt.show()

