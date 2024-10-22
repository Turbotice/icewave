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


from pprint import pprint
import icewave.tools.datafolders as df
import glob

import icewave.tools.rw_data as rw_data


#%%  Select data and save it 

global drones
drones = ['mesange','Bernache','Fulmar']
global base
base = df.find_path(disk='Hublot24')

def get_jpgfiles(date):
    base = df.find_path(disk='Hublot24')
    #drones = ['mesange', 'Bernache',]
    jpgfiles = {}
    print(base)
    for key in drones:  
        jpg = glob.glob(base+date+'/Drones/'+key+'/*/*.jpg')#/*/*.srt')
        if len(jpg)>0:
            jpgfiles[key] = jpg
        else:
            print(f"No data for {key} on {date}")
    return jpgfiles

def get_jpg_record(jpgfile,drone='mesange'):
    base = df.find_path(disk='Hublot24')
    nbase = len(base)
    record={}
    if drone=='mesange':
        h0 = -1
    elif drone=='bernache':
        h0 = 5
    elif drone=='Fulmar':
        h0 = 5
    else:
        h0 = 0
        print(key)
        print('Drone unknown')
    record = {}
    name = jpgfile.split('/')[-2]#.split('.')[0]
            # print(i,files,name)
    record['name']=jpgfile.split('/')[-1].split('.')[0]
    record['path']=jpgfile[nbase:].split('.')[0]
    record['time'] = time_BA_to_SP(int(jpgfile[-17:-11]) + h0 * 10000)
    record['folder'] = name
    return record


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


def time_BA_to_SP(time) :
    return str(int(time))[:2] + ':' + str(int(time))[2:4] + ':' + str(int(time))[4:]

def create_timeline(flight_record) :
    timeline = []
    for i in flight_record['CUSTOM.updateTime [local]'] :
        timeline += [PM_to_normal(i)[:-3]]
    return timeline

def extract_time_fr(flight_record, t0, tf) :
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
    
def extract_all_fr(date, drone, selection = True) :
    #avec une date et un drone, on cree record avec toutes les valeurs interessantes et au temps d'enregistrement des cameras depuis les flight record
    #base = disque + ':/Share_hublot
    base = df.find_path(disk='Hublot24')

    path_SP = base + date +'/Summary/records_' + date + '.pkl'
    path_fr = base + date +'/Drones/'  + drone+'/flightrecords/Flightrecord_dict.pkl'  
    path_selection = base+'/Summary/' + date + '_path_drone.txt'
    
    flight_record = open_dico(path_fr)
    summary_SP = open_dico(path_SP)
    select = pandas.read_csv(path_selection, header= None)
    select = base + np.asarray(select)[:,0]
    
    datas = []
    
    for i in select :
        datas += [i.split('/')[-1]]
    
    record = {}
    
    for folder in summary_SP['drones'][drone].keys() :
        if folder in datas or selection == False :
            record[folder] = []
            for j in range( len(summary_SP['drones'][drone][folder])) :
            
                t0 = min(summary_SP['drones'][drone][folder][j]['time'])
                tf = max(summary_SP['drones'][drone][folder][j]['time'])
                
                name = summary_SP['drones'][drone][folder][j]['name']
                
                record[folder] += [extract_time_fr(flight_record, t0, tf)]
                record[folder][j]['name'] = name
    return record
                    
def save_record(date, drone, record) :
    #save le record là où il y a les données
    base = df.find_path(disk='Hublot24')
    for key in record.keys() :
        path = base + date + '/Drones/' + drone + '/' + key + "/record_" + date + '_' + drone + '_' + key + '.pkl'
        print(path)
        save_dico(record[key], path )
        print('saved')
        
def MAIN(date, drone, selection = True, save = False) :
    record = extract_all_fr(date, drone, selection = selection)
    
    if save :
        save_record(date, drone, record)
        
    return record


