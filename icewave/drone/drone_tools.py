# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:19:04 2026

@author: sebas

Gathers all functions used to drone analysis such as
- collecting initial time of videos, 


"""

import glob
from datetime import datetime
import pytz

def get_UTC0_from_SRT(path2SRT,drone_ID,exp_ID,year = '2024'):
    """ Compute UTC_t0, beginning of movie from SRT file"""
    
    filelist = glob.glob(f'{path2SRT}*.SRT')
    print(filelist)
    file_SRT = filelist[0].split(exp_ID)[1]
    date_string = file_SRT.split('_')[1]

    format_date = '%Y%m%d%H%M%S'
    current_time = datetime.strptime(date_string,format_date)
    print(current_time)

    if drone_ID == 'mesange':
        UTC_t0 = current_time.astimezone(pytz.utc)
    elif drone_ID == 'bernache':
        current_time = current_time.replace(tzinfo = pytz.timezone('America/Montreal'))
        UTC_t0 = current_time.astimezone(pytz.utc)
        
    return UTC_t0