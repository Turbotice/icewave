# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:15:20 2015

@author: stephane
"""

import os.path
import numpy as np
import h5py
import csv

#to use to read data files in ASCII formats (or others ?)
#to use to write data files in ASCII formats (or others ?) with label on the top :
# in particular to create a catalog of the parameters of all the existing data (!)


def read_xml():
    pass

def read_csv(filename,delimiter=','):
    rows = []
    with open(filename,'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=delimiter, quotechar='|')
        for row in spamreader:
            rows.append(row)
    return rows

def write_csv(filename,data):
    folder = os.path.dirname(filename)
    if not os.path.exists(folder):
        print(f"Creating folder {folder}")
        os.makedirs(folder)

    with open(filename, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',quotechar='|')
        keys = list(data.keys())
        spamwriter.writerow(keys)
        n = len(data[keys[0]])

        for i in range(n):
            line = [data[key][i] for key in keys]
            spamwriter.writerow(line)

def writedict_csv(filename,data,symbol='#'):
    with open(filename, 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',quotechar='|')#, quoting=csv.QUOTE_MINIMAL)

        keys = list(data.keys())
        print(keys)
        header = [symbol]+list(data[keys[0]].keys())
        spamwriter.writerow(header)
        for key in data.keys():
            row = [key]+[data[key][k] for k in data[key].keys()]
            spamwriter.writerow(row)

def write_h5(filename,data):
    hf = h5py.File(filename, 'w')
    for key in data.keys():
        if type(data[key])==dict:
            hf.create_group(key)
            for k in data[key].keys():
                hf[key].create_dataset(k,data=data[key][k])
        else:
            hf.create_dataset(key, data=np.asarray(data[key]))
    hf.close()

def save_h5_rec(filename,data):
    with h5py.File(filename, 'w') as hf:
        hf =  write_h5_rec(hf,data)
    
def write_h5_rec(hf,data):
    for key in data.keys():
        if type(data[key])==dict:
            hf.create_group(key)
            print(key)
            hf = write_h5_rec(hf[key],data[key])
        else:
            print(data[key])
            hf.create_dataset(key, data=np.asarray(data[key]))
    return hf

def read_h5(filename):
    hf = h5py.File(filename,'r')
    return hf

def write_pkl(filename,data):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_pkl(filename):
    import pickle
    with open(filename, 'rb') as handle:
        handle.seek(0)
        data = pickle.load(handle)#, protocol=pickle.HIGHEST_PROTOCOL)
    return data

def csv2dict(table,headerindex=0,symbol='#'):
    data = {}
    if table[0][0]==symbol:
        keys = table[0][1:]
        #print(keys)
        for tab in table[1:]:
            try:
                #try to convert to int the key
                tab[0]=int(tab[0])
            except:
                pass #do nothing
            #print(tab)
            data[tab[0]]={}
            for (t,key) in zip(tab[1:],keys):
                if '.' in t:
                    data[tab[0]][key]=float(t)
                else:
                    try:
                        data[tab[0]][key]=int(t)
                    except:
                        data[tab[0]][key]=str(t)
    else:
        header = table[headerindex]
        data = {}
        for key in header:
            data[key]=[]
        for i in range(headerindex+1,len(table)):
            for j,key in enumerate(header):
                data[key].append(table[i][j])
    return data


def save_dict_to_h5(dict_obj,filename):
    """ Save dictionnary in HDF5 file """
    with h5py.File(filename,'w') as hf:
        write_dict_h5_rec(hf,dict_obj)
        
def write_dict_h5_rec(h5group,dict_obj):
    """ Write a dictionary (dict_obj) to an open h5py.File object. This recursive function converts
    numpy arrays, arrays of strings, list of strings and dictionnaries to equivalents h5 structures """

    for key,value in dict_obj.items():
        
        if isinstance(value,dict): # if value is a dictionnary
            subgroup = h5group.create_group(key)
            write_dict_h5_rec(subgroup,value)
            
        elif isinstance(value, list):# if value is a list
            list_group = h5group.create_group(key)
            for i, item in enumerate(value):
                item_key = str(i)

                # Nested dictionary
                if isinstance(item, dict):
                    subgroup = list_group.create_group(item_key)
                    write_dict_h5_rec(subgroup, item)

                # Nested numpy array
                elif isinstance(item, np.ndarray):
                    if item.dtype.kind in {'U','O'} and np.all(np.vectorize(lambda x: isinstance(x,str))(item)):
                        dt = h5py.string_dtype(encoding='utf-8')
                        list_group.create_dataset(item_key, data=item.astype(object), dtype=dt)
                    else:
                        list_group.create_dataset(item_key, data=item)

                # Nested string
                elif isinstance(item, str):
                    dt = h5py.string_dtype(encoding='utf-8')
                    list_group.create_dataset(item_key, data=item, dtype=dt)

                # Nested scalar
                elif isinstance(item, (int, float, np.integer, np.floating)):
                    list_group.create_dataset(item_key, data=item)

                # Nested list => recurse
                elif isinstance(item, list):
                    sublist_group = list_group.create_group(item_key)
                    # recursive list handling
                    for j, subitem in enumerate(item):
                        subkey = str(j)
                        if isinstance(subitem, dict):
                            g = sublist_group.create_group(subkey)
                            write_dict_h5_rec(g, subitem)
                        else:
                            # call the same handler: wrap the single element as dict
                            write_dict_h5_rec(sublist_group, {subkey: subitem})    
            
        elif isinstance(value,np.ndarray): # if value is a numpy array
            # case 1 : numpy array of strings
            if value.dtype.kind in {'U','O'} and np.all(np.vectorize(lambda x: isinstance(x,str))(value)):
                dt = h5py.string_dtype(encoding = 'utf-8')
                h5group.create_dataset(key, data = value.astype(object),dtype = dt)
                
            # case 2 : any other numpy array 
            else :
                h5group.create_dataset(key, data = value)
                
        elif isinstance(value,str): # if value is a string 
            dt = h5py.string_dtype(encoding = 'utf-8')
            h5group.create_dataset(key, data = value,dtype = dt)
            
        elif isinstance(value,list) and all(isinstance(v,str) for v in value): # if value is a list of strings
            dt = h5py.string_dtype(encoding = 'utf-8')
            h5group.create_dataset(key, data = np.array(value,dtype = object),dtype = dt)
            
        elif isinstance(value, (int, float, np.integer, np.floating)):
                # store scalars directly
                h5group.create_dataset(key, data=value)
            
        else:
            raise TypeError(f"Unsupported data type for key '{key}': {type(value)}")
            
def load_dict_from_h5(filename):
    """ Load dictionnary saved in a HDF5 file """
    with h5py.File(filename,'r') as hf:
        out = load_dict_h5_rec(hf)
        return out
        
def load_dict_h5_rec(h5group):
    """ Load dictionnary from an open h5py.file object. """
    out = {}
    
    # First detect if this group represents a LIST
    # HDF5 groups used to store lists have numeric keys: "0", "1", "2", ...
    # if all(k.isdigit() for k in h5group.keys()) and len(h5group.keys()) > 0:
    #     # → This group is a list
    #     print('This group is a list')
    #     out_list = []
    #     for idx in sorted(h5group.keys(), key=lambda x: int(x)):
    #         item = h5group[idx]

    #         if isinstance(item, h5py.Group):
    #             out_list.append(load_dict_h5_rec(item))

    #         elif isinstance(item, h5py.Dataset):
    #             data = item[()]

    #             if isinstance(data, bytes):
    #                 out_list.append(data.decode("utf-8"))

    #             elif isinstance(data, np.ndarray) and data.dtype.kind in {"U", "S", "O"}:
    #                 list_string = [
    #                     x.decode("utf-8") if isinstance(x, bytes) else str(x)
    #                     for x in data
    #                 ]
    #                 out_list.append(np.array(list_string))

    #             else:
    #                 out_list.append(data)

    #     return out_list
    
    for key,item in h5group.items():
        if isinstance(item,h5py.Group): # if item is a h5 group 
            out[key] = load_dict_h5_rec(item)
            
        elif isinstance(item,h5py.Dataset): # if item is a h5 dataset 
            data = item[()]
            
            if isinstance(data,bytes): # decode single string
                out[key] = data.decode('utf-8')
            
            elif isinstance(data,np.ndarray) and data.dtype.kind in {'U','S','O'}: # decode arrays of strings
                list_string = [x.decode('utf-8') if isinstance(x,bytes) else str(x) for x in data]
                out[key] = np.array(list_string)
            
            else : 
                out[key] = data
                
    return out
