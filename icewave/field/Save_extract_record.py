# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:04:03 2024

@author: Banquise
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:25:04 2024

@author: Banquise
"""
import pandas
import pickle 
import numpy as np
import os


import icewave.tools.datafolders as df




#%%  Select data and save it 


def open_dico(path):
    a_file = open(path, "rb")
    dico = pickle.load(a_file)
    a_file.close
    return dico


def save_dico(dico, path):
    a_file = open(path, "wb")
    pickle.dump(dico, a_file)
    a_file.close()

def PM_to_normal(time_fr):
    if time_fr.split()[-1] == 'PM' :
        if int(time_fr.split(':')[0]) > 11 :
            hour = str(int(time_fr.split(':')[0]))
        else :
            hour = str(int(time_fr.split(':')[0]) + 12)
        hh = ''
        for i in time_fr.split(':')[1:] :
            hh += ':' +  i
        time = hour + hh[:-3]
    
    elif time_fr.split()[-1] == 'AM' :
        time = time_fr.split()[0]
        
    return time


def create_timeline(flight_record) :
    timeline = []
    for i in flight_record['CUSTOM.updateTime [local]'] :
        timeline += [PM_to_normal(i)[:-3]]
    return timeline



def extract_time_f_r(flight_record, t0, tf) :
    """
    

    Parameters
    ----------
    flight_record : Dict
        Flight record 
    t0 : "12:46:22"
        DESCRIPTION.
    tf : "12:48:22"
        DESCRIPTION.

    Returns
    -------
    new_fr : dict
        flight record between t0 and tf.

    """
    timeline = np.array(create_timeline(flight_record)) # cree une timeline de tous les temps existants avec le bon format (pas avec les PM)
    new_fr = {}
    for key in flight_record.keys() :
        new_fr[key] = []
    for i in range (len(timeline)) :
        if timeline[i] <= tf and timeline[i] >= t0 :
            for key in flight_record.keys() :
                new_fr[key] += [flight_record[key][i]]
            
    return new_fr
    


def extract_from_fr(date, drone, selection = True, disque = 'K') :
    #avec une date et un drone, on cree record avec toutes les valeurs interessantes et au temps d'enregistrement des cameras depuis les flight record
    #base = disque + ':/Share_hublot
    base = df.find_path(disk='Hublot24')

    path_SP = base + date +'/Summary/records_' + date + '.pkl'
    path_fr = base + date +'/Drones/'  + drone+'/flightrecords/Flightrecord_dict.pkl'    
    flight_record = open_dico(path_fr)
    summary = open_dico(path_SP)
    path_selection = base+'/Summary/' + date + '_path_drone.txt'
    select = pandas.read_csv(path_selection, header= None)
    select = 'K' + np.asarray(select)[:,0]
    datas = []
    for i in select :
        # if i.split('//')[-2] == drone : #on garde les noms de dossiers du bon drone, ATTENTION  écrire les drones en minuscule
        datas += [i.split('/')[-1]]
    
    
    record = {}
    
    for folder in summary['drones'][drone].keys() :
        if folder in datas or selection == False :
            record[folder] = {}
            for j in range( len(summary['drones'][drone][folder])) :
            
                t0 = min(summary['drones'][drone][folder][j]['time'])
                tf = max(summary['drones'][drone][folder][j]['time'])
                
                key_name = summary['drones'][drone][folder][j]['name']
                
                record[folder][key_name] = extract_time_f_r(flight_record, t0, tf)
                


    
    return record
                    


def save_record(date, drone, record, disque = 'K') :
    #save le record là où il y a les données
    base = df.find_path(disk='Hublot24')
    for key in record.keys() :
        path = base + date + '/Drones/' + drone + '/' + key + "/record_" + date + '_' + drone + '_' + key + '.pkl'
        print(path)
        save_dico(record[key], path )
        print('saved')
        
def MAIN(date, drone, selection = True, disque = 'K', save = False) :
    record = extract_from_fr(date, drone, selection = selection, disque = disque)
    
    if save :
        save_record(date, drone, record, disque = disque)
        
    return record


