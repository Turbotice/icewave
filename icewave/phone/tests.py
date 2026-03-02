import glob
import os
import numpy as np
import time
import datetime 
import icewave.phone.load as load

from pprint import pprint
global path
path = 'Test_data_FP3/'
global path_android
path_android = '/storage/self/primary/Download/Gobannos/'

def data_stats():
    #load all the files in the path folder
    data = load.load_folder(path)
    data = load.sync_time(data)
    load.stat(data)

    return data

def parse_files(date=None):
    if date==None:
        date = str(datetime.date.today())

    filelist = glob.glob(path_android+date+'*.csv')

    experiments={}
    for filename in filelist:
        get_start_time()
    pprint(filelist)
    #parse  per start time

def get_start_time(filename):
    filename = os.path.basename(filename)
    filename = '20'+filename.split('20')[1]
    '-'.join(filename.split('-')[:3])
             
#parse_files()
#data_stats()




