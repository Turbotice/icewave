# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 10:03:32 2025

@author: sebas

Module Python gathering all useful functions to deal with time, time stamps, UTC time from field work 

"""


import numpy as np 
from time import strftime, localtime
from datetime import datetime
import pytz

def epoch2datetime(t,timezone = None,format_date = '%Y-%m-%d %H:%M:%S.f'):
    """ Convert time since epoch to UTC time 
    Inputs: - t, array or list of times since epoch
            - timezone, pytz timezone in which we want datetime 
            - format_date, string, format in which datetime will be written """  
    
    t_datetime = []
    for t_epoch in t:
        current_string = strftime(format_date, localtime(t_epoch))
        current_time = datetime.strptime(current_string, format_date)
        if timezone is None:
            adjusted_time = current_time
        else:
            adjusted_time = current_time.astimezone(timezone)
            
        t_datetime.append(adjusted_time)
    
    t_datetime = np.array(t_datetime)
    t_datetime = np.reshape(t_datetime,(t.shape))
    return t_datetime


def generate_datetime_txt(UTC_0):
    """ Generate a string from a datetime UTC_0"""
    
    txt = f'{UTC_0.year}_{UTC_0.month}_{UTC_0.day}_{UTC_0.hour}_{UTC_0.minute}_{UTC_0.second}'
    return txt