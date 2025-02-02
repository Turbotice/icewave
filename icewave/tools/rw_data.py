# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:15:20 2015

@author: stephane
"""

import os.path
import numpy as np
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

def write_pkl(filename,data):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_pkl(filename):
    import pickle
    with open(filename, 'rb') as handle:
        data = pickle.load(handle)#, protocol=pickle.HIGHEST_PROTOCOL)
    return data

def csv2dict(table,headerindex=0,symbol='#'):
    data = {}
    if table[0][0]==symbol:
        keys = table[0][1:]
        print(keys)
        for tab in table[1:]:
            print(tab)
            data[tab[0]]={}
            for (t,key) in zip(tab[1:],keys):
                if '.' in t:
                    data[tab[0]][key]=float(t)
                else:
                    data[tab[0]][key]=int(t)
    else:
        header = table[headerindex]
        data = {}
        for key in header:
            data[key]=[]
        for i in range(headerindex+1,len(table)):
            for j,key in enumerate(header):
                data[key].append(table[i][j])
    return data
