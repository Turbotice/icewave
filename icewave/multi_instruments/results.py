import pylab as plt
import numpy as np

#import icewave.field.multi_instruments as multi
#import icewave.geometry.display as display
#import icewave.display.graphes as graphes

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

def get_base(disk='Elements',year='2024'):
    return df.find_path(disk=disk,year=year)

def get_record(date,year='2024'):
    base = get_base(year=year)
    filename =  base+f'{date}/Summary/records_{date}.pkl'
    return rw_data.load_pkl(filename)

def make_result(date,typ,instrument,name,var,value,index=0,time=None,records=None,year='2024',comments=''):
    """ Create a dictionnary of results, where each result is associated to a single key build as followed :
        key = year_dateTtime_typ_instrument_name_var_index
        
        Inputs : - date : str, date of experiment, exemple : '0226' for 26th February
                 - typ : str, type of instrument, exemple : drone, gps, phone, buoy, etc...
                 - instrument : str, name of instrument, exemple : bernache, mesange, garmin_SP, etc...
                 - name : str, name of sampling, exemple : H_01, 10-waves_12, etc...
                 - var : str, name of variable measured : H, A, f0, etc...
        
        Optional inputs : - index : int, index of measurement of the corresponding variable
                          - time : str, time of measurement, can be extracted from records files. 
                                   Should be at the format : hh:mm:ss
                          - records : dictionnary, table of records parameters 
                          - year : str
                          - comments : str, additional comments to this variable 
                          
        Outputs : - results : a dictionnary containing a single key (described above) and several 'sub_keys' which are :
                    + year, date, time, 
                    + latitude, longitude, 
                    + var, value of the measured variable, can also be an array 
                    + path, path to raw data 
                    + name_instrument, type_instrument and comments """
    
    results={}
    if records==None:
        print('retrieve associated records')
        records = get_record(date,year=year)
    try:
        record = records[typ][instrument][name]
        if type(record)==dict:
            pass
            #print(record.keys())
        elif type(record)==list:
            record = record[index]
    except:
        print('key not found in records !')
        print(typ,instrument,name)

    if time==None:
        print('No time specified !')
        time = record['time'][index]
        
    tstamp = year+'_'+date+'T'+str(time)
    key = '_'.join([tstamp,typ,instrument,name,var,str(index)])

    results[key] = {}
    results[key]['year']=year
    results[key]['date']=date
    results[key]['time']=time
    results[key]['latitude']=record['latitude'][index]
    results[key]['longitude']=record['longitude'][index]

    if value=='auto':
        #try:
        #print(record.keys())
        results[key][var]=float(record['params'][var])
        #except:
        #    print('key var does not exist in record')
    else:
        results[key][var]=value

    results[key]['name_instrument']=instrument
    results[key]['type_instrument']=typ

    if 'path' in record.keys():
        results[key]['path']=record['path']
    else:
        print(f'No path available for {key}')
    results[key]['comments']=comments    
    #results[key]['path']=None

#    record = records[typ][instrument][name]

    return results

def save_result(result):
    key = list(result.keys())[0]

    typ = result[key]['type_instrument']
    year = result[key]['year']
    date = result[key]['date']
    
    base = get_base(year=year)
    filename =  base+f'{date}/Summary/results_{typ}_{date}.pkl'

    rw_data.write_pkl(filename,result)
