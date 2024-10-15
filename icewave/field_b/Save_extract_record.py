# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:25:04 2024

@author: Banquise
"""
import pandas
import pickle 
import numpy as np

 
#%%


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


def extract_from_fr(date, drone, liste = True, selection = True, disque = 'K') :
    
    #avec une date et un drone, on cree record avec toutes les valeurs interessantes et au temps d'enregistrement des cameras depuis les flight record
    
    path_SP = disque + ':\Share_hublot\\Data\\' + date +'\\Summary\\records_' + date + '.pkl'
    path_fr = disque + ':\Share_hublot\\Data\\' + date +'\\Drones\\'  + drone+'\\flightrecords\\Flightrecord_dict.pkl'    
    flight_record = open_dico(path_fr)
    summary = open_dico(path_SP)

    path_selection = disque + ':\Share_hublot\\Summary\\' + date + '_path_drone.txt'
    select = pandas.read_csv(path_selection, header= None)
    select = 'K' + np.asarray(select)[:,0]
    datas = []
    for i in select :
        # if i.split('\\')[-2] == drone : #on garde les noms de dossiers du bon drone, ATTENTION  écrire les drones en minuscule
        datas += [i.split('\\')[-1]]
    

    timeline = np.array(create_timeline(flight_record)) # cree une timeline de tous les temps existants avec le bon format (pas avec les PM)
    
    keys_date = ['CUSTOM.date [local]', 'CUSTOM.updateTime [local]']
    keys_bool = ['CAMERA.isPhoto', 'CAMERA.isVideo']
    #, 'OSD.flyTime', 'OSD.flyTime [s]', 'OSD.latitude', 'OSD.longitude', 'OSD.height [ft]']
    keys1 = ['OSD.altitude [ft]', 'OSD.mileage [ft]', 'OSD.hSpeed [MPH]', 'OSD.xSpeed [MPH]', 'OSD.ySpeed [MPH]', 'OSD.zSpeed [MPH]']
    keys2 = ['OSD.pitch', 'OSD.roll', 'OSD.yaw', 'OSD.yaw [360]', 'OSD.gpsNum']
    keys3 = ['GIMBAL.pitch', 'GIMBAL.roll', 'GIMBAL.yaw', 'GIMBAL.yaw [360]']            
    keys_float = keys1+keys2+keys3
    
    record = {}
    
    
    
    for j in summary['drones'][drone].keys() :
        if j in datas or selection == False :
            record[j] = {}
            print('Dossier :', j)
            if liste : # On peut ajouter des keys dans les keys ou enlever j (nom du dossier),  rajouter le path pour sauvegarder au bon endroit
                for key in keys_float:
                    record[j][key] = []
                for key in keys_date:
                    record[j][key] = []
                for key in keys_bool:
                    record[j][key] = []
                    
            else : #Fait tout en array
                for key in keys_float:
                    record[j][key] = np.array([])
                for key in keys_date:
                    record[j][key] = np.array([])
                for key in keys_bool:
                    record[j][key] = np.array([])
            for k in summary['drones'][drone][j]['time'] :

                if k in timeline : #si le temps du summary est dans la timeline, on cherche les indices et on garde toutes les keys interessantes
                    indices = [np.where(timeline == k)[0]]
    
                    if liste :
                        
                        for key in keys_float:
                            record[j][key] += np.array(flight_record[key])[indices].tolist()[0]
                        for key in keys_date:
                            record[j][key] += np.array(flight_record[key])[indices].tolist()[0]
                        for key in keys_bool:
                            record[j][key] += np.array(flight_record[key])[indices].tolist()[0]
                    
                    else : #Fait tout en array
                        for key in keys_float:
                            record[j][key] = np.append(record[j][key], np.array(flight_record[key])[indices])
                        for key in keys_date:
                            record[j][key] = np.append(record[j][key], np.array(flight_record[key])[indices])
                        for key in keys_bool:
                            record[j][key] = np.append(record[j][key], np.array(flight_record[key])[indices])
                            
    return record
                    

def save_record(date, drone, record, disque = 'K') :
    #save le record là où il y a les données
    for key in record.keys() :
        path = disque + ':\Share_hublot\\Data\\' + date + '\\Drones\\' + drone + '\\' + key + "\\record_" + date + '_' + drone + '_' + key + '.pkl'
        print(path)
        save_dico(record[key], path )
        print('saved')
        
    

def MAIN(date, drone, liste = True, selection = True, disque = 'K', save = False) :
    record = extract_from_fr(date, drone, liste = liste, selection = selection, disque = disque)
    
    if save :
        save_record(date, drone, record, disque = disque)
        
    return record
        
date = '0223'
drone = 'mesange' #Bernache ou bernache ?????

record = MAIN(date, drone, save = False)




