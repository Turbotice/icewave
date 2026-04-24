#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd

#%%
def Kice_from_Keff_Ksticks(Keff, Kstick, nsticks=2):
    return 1/((1/Keff)-(nsticks/Kstick))

def Kstick_from_reflength(L_mm, Kstick_ref, reflength_mm=200):
    return (Kstick_ref * L_mm)/reflength_mm

def Eice_approx_from_Kice_dimsample(Kice, L_mm, D_mm):
    return (16/(3*np.pi*D_mm*1e-3*L_mm*1e-3))*Kice

def Eice_approx_from_Keff_dimsample(Keff, L_mm, D_mm, Kstick_ref=63*1e6, reflength_mm=200, nsticks=2):
    Kstick = Kstick_from_reflength(L_mm, Kstick_ref, reflength_mm=reflength_mm)
    Kice = Kice_from_Keff_Ksticks(Keff, Kstick, nsticks=nsticks)
    Eice = Eice_approx_from_Kice_dimsample(Kice, L_mm, D_mm)
    return Eice

#%%
# test
# load for one acquisition the corresponding dataset (pkl file)
disk = 'H:'
date = '0223'
dictdir = f'{disk}/data/{date}/Tests_Bresiliens/results'
inputfile = f'{dictdir}/fits_force_displacement.pkl'
outputfile = f'{dictdir}/fits_force_displacement_correcwood.pkl'

with open(inputfile, 'rb') as f:
    inputdict = pickle.load(f)

outputdict = inputdict.copy()
print(outputdict)

# %%
def dimsamples_from_csv_summary(csvfilepath=f'{disk}/data/{date}/Tests_Bresiliens/Summary_samples_tests.csv'):
    df = pd.read_csv(csvfilepath, sep=';', encoding='latin1')
    array_L_mm = df['L_mm']
    array_L_err_mm = df['L_err_mm']
    array_D_mm = df['D_mm']
    array_acq = df['acq_num']
    dd = {}
    for i in range(len(array_acq)):
        dd[array_acq[i]] = {}
        dd[array_acq[i]]['L_mm'] = array_L_mm[i]
        dd[array_acq[i]]['L_err_mm'] = array_L_err_mm[i]
        dd[array_acq[i]]['D_mm'] = array_D_mm[i]
    return dd

def correct_fit_params_1acq(dimsamples_infos, acq='2', L_mm=185, subdict_acq=outputdict['2']):
    
    original_keys = [k for k in subdict_acq.keys() if 'fit' in k and 'correc_wood' not in k]
    for k in original_keys:    
        if 'fit' in k:
            a_kN_per_mm = subdict_acq[k]['a']
            Keff = a_kN_per_mm * 1e6
            D_mm = dimsamples_infos[acq]['D_mm']
            L_mm = dimsamples_infos[acq]['L_mm']
            L_err_mm = dimsamples_infos[acq]['L_err_mm']
            Eice_approx = Eice_approx_from_Keff_dimsample(Keff=Keff, L_mm=L_mm, D_mm=D_mm)

            subdict_acq['L_mm'] = L_mm
            subdict_acq['L_err_mm'] = L_err_mm
            subdict_acq['D_mm'] = D_mm

            subdict_acq[k+'_correc_wood'] = {}
            subdict_acq[k+'_correc_wood']['Eice_approx'] = Eice_approx
            # pour l'instant on ne sembête pas avec les propagations d'incertitude...


def correct_fit_params_allacqs(csvfilepath=f'{disk}/data/{date}/Tests_Bresiliens/Summary_samples_tests.csv'):
    dimsamples_infos = dimsamples_from_csv_summary(csvfilepath=csvfilepath)
    for i in range(len(dimsamples_infos.keys())):
        acq = list(dimsamples_infos.keys())[i]
        L_mm = dimsamples_infos[acq]['L_mm']
        if acq in outputdict:
            correct_fit_params_1acq(dimsamples_infos=dimsamples_infos, acq=acq, L_mm=L_mm, subdict_acq=outputdict[acq])
        else:
            pass

correct_fit_params_allacqs()

with open(outputfile, 'wb') as f:
    pickle.dump(outputdict, f)

# %%
