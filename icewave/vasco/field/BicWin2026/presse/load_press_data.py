#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_single_test_curve(acqnum='4', date='0223', disk='H:',plot=False):
    # load csv data
    df = pd.read_csv(f"{disk}/data/{date}/Tests_Bresiliens/presse_data/{acqnum}.csv",
                sep=',',
                skiprows=6,
                header=[0,1],
                    encoding='latin-1')

    dict_single_test = {
        'date' : date,
        'acqnum' : acqnum,
        'units' : "(Actuator : mm , Load : kN , Time : sec)",
        }

    for i in range(len(df.keys())):
        key = df.keys()[i]
        new_key = key[0].replace(' ','')
        dict_single_test[new_key] = df[key]
    if plot:
        plt.figure()
        plt.plot(dict_single_test['Actuator'].values, dict_single_test['Load'].values)
        plt.xlabel('Displacement (mm)', fontsize=13)
        plt.ylabel('Load (kN)', fontsize=13)
        plt.title('Applied force vs time', fontsize=13)
        plt.show()
    return dict_single_test 

#%%
# exemple pour afficher les données
"""%matplotlib qt


"""
