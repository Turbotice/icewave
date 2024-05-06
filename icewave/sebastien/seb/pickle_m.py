# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:59:44 2023

@author: sebas
"""
import pickle
import os
import sys

def write(obj,filename):  
    
    try:
        f=open(filename,'wb')
    except 'EOFError':
        print('Empty file')
        return None
    p = pickle.Pickler(f,2)
#        S_str=self.decode()
    try:
        p.dump(obj)    
    except '_pickle.PicklingError':
        print('Sdata class has been modified')
    f.close()    
    
def read(filename):
    if os.path.isfile(filename):
       # print("Reading Sdata")# from "+filename)
        #has to be secured with a unic identifier.
        #bad method : the filename is used to compare with the attributes previously used to generate it
        v=sys.version_info

        f=open(filename,'rb')
        if v[0]==3:
            buf=f.read()
        
            S=pickle.loads(buf,encoding='latin1')
        else:
            S=pickle.load(f)
            
        f.close()
        
        return S
    else:
        print(filename+ "does not exist")
        return None