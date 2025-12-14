import glob
import os
import numpy as np
import time

import icewave.phone.load as load

global path
path = 'Test_data_FP3/'

def data_stats():
    #load all the files in the path folder
    data = load.load_folder(path)
    data = load.sync_time(data)
    load.stat(data)

    return data
