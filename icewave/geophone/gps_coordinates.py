# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 17:22:26 2025

@author: sebas
"""

import os
import re
import numpy as np 
from datetime import datetime, time , timedelta
import pytz
import glob 


# Define the find_coordinates function
def find_coordinates(file_path, date, year):
    """ Collect coordinates of all acquisition in a given DigiSolo.txt file """
    
    num_snippets = 0
    target_date_found = False
    target_date_pattern = year + '/' + date[:2] + '/' + date[2:]
    GPS_coordinates_matrices = []
    current_snippet = ""

    with open(file_path, 'r') as file:
        for line in file:
            if target_date_found and current_snippet:
                current_snippet += line
                if line.startswith("[Notify00001]"):
                    num_snippets += 1
                    latitudes = [float(value) for value in re.findall(r'Latitude\s*=\s*(-?\d+\.\d+)', current_snippet)]
                    longitudes = [float(value) for value in re.findall(r'Longitude\s*=\s*(-?\d+\.\d+)', current_snippet)]
                    num_GPS_logs = min(len(latitudes), len(longitudes))
                    UTC_times = re.findall(r'Synchronization\nUTC Time = "(\d+/\d+/\d+,\d+:\d+:\d+)"',current_snippet)
                    UTC_t = [datetime.strptime(t,'%Y/%m/%d,%H:%M:%S') for t in UTC_times]
                    # lat_avg = sum(latitudes) / len(latitudes) if latitudes else None
                    # long_avg = sum(longitudes) / len(longitudes) if longitudes else None
                    GPS_coordinates_matrices.append({
                        "acquisition": num_snippets,
                        "num_GPS_logs": num_GPS_logs,
                        "latitude": latitudes,
                        "longitude": longitudes,
                        "UTC_t":UTC_t
                    })
                    current_snippet = ""
            elif target_date_found and line.startswith("Start Acquisition FileName ="):
                current_snippet += line.partition("Start Acquisition FileName =")[2]
            elif line.startswith("UTC Time ="):
                match = re.search(r'UTC Time = "(\d{4}/\d{2}/\d{2})', line)
                if match and match.group(1) >= target_date_pattern:
                    target_date_found = True

        # If the file ends with an acquisition snippet
        if current_snippet:
            num_snippets += 1
            latitudes = [float(value) for value in re.findall(r'Latitude\s*=\s*(-?\d+\.\d+)', current_snippet)]
            longitudes = [float(value) for value in re.findall(r'Longitude\s*=\s*(-?\d+\.\d+)', current_snippet)]
            num_GPS_logs = min(len(latitudes), len(longitudes))
            UTC_times = re.findall(r'Synchronization\nUTC Time = "(\d+/\d+/\d+,\d+:\d+:\d+)"',current_snippet)
            UTC_t = [datetime.strptime(t,'%Y/%m/%d,%H:%M:%S') for t in UTC_times]
            # lat_avg = sum(latitudes) / len(latitudes) if latitudes else None
            # long_avg = sum(longitudes) / len(longitudes) if longitudes else None
            GPS_coordinates_matrices.append({
                "acquisition": num_snippets,
                "num_GPS_logs": num_GPS_logs,
                "latitude": latitudes,
                "longitude": longitudes,
                "UTC_t":UTC_t
            })

    return num_snippets, GPS_coordinates_matrices

# Main function
def get_GPS_coordinates(path2logfile,geophones_table_path,date,year):
    """ Build a dictionnary of GPS coordinates for all geophones and all acquisitions for a given date. 
    Geophones are refered using keys of the geophones_table_path : 'G01' to 'G16'. 
    Acquisition number is indexed from 1 """
    log_files = [file for file in os.listdir(path2logfile) if file.startswith("DigiSolo") and file.endswith(".txt")]
    print(log_files)
    # Sort the log files based on DigiSolo number
    log_files.sort(key=lambda x: int(re.search(r'\d+', x).group()))

    all_matrices = {}

    for log_file in log_files:
        log_file_path = os.path.join(path2logfile, log_file)
        num_snippets, GPS_coordinates_matrices = find_coordinates(log_file_path, date, year)
        
        acquisitions_dict = {entry['acquisition']: {
            'num_GPS_logs': entry['num_GPS_logs'],
            'latitude': entry['latitude'],
            'longitude': entry['longitude'],
            'UTC_t' : entry['UTC_t']
        } for entry in GPS_coordinates_matrices}

        all_matrices[log_file] = acquisitions_dict
        
    # rename dictionnary using geophone labels 'G01 to G16' thanks to geophones table path
    new_matrix = {}

    geophones_table = {}
    with open(geophones_table_path, 'r') as table_file:
        # Skip the header
        next(table_file)
        for line in table_file:
            num, geo_sn = line.strip().split('\t')
            geophones_table[geo_sn] = num

    # create a dictionnary with ordered geophones location 
    for key in all_matrices.keys():
        file_number_match = re.search(r'DigiSolo_(\d{4})',key)
        file_number = file_number_match.group(1)
        geo_sn = "45302" + file_number.zfill(4)
        
        # Look up the corresponding two-digit number from the geophones table
        if geo_sn in geophones_table:
            two_digit_number = geophones_table[geo_sn]
            label = f"G{two_digit_number}"
        else:
            label = f"Unknown"
            
        new_matrix[label] = all_matrices[key]
      
    new_matrix = dict(sorted(new_matrix.items()))

    return new_matrix

def compute_avg_logs(GPS_logs):
    """ Compute averaged GPS position from a dictionnary of GPS positions and time 
    Input : - GPS_logs, dictionnary containing keys : 'latitude', 'longitude', 'UTC_t' (list of datetime object at which
                                                                                        GPS position has been recorded)
    Output : - avg_logs, dictionnary containing keys : 'latitude', 'longitude', 'UTC_t' and a single value for each key """
    
    avg_logs = {}
    # Compute averaged position of geophones 
    for key in ['longitude','latitude']:
        avg_logs[key] = np.mean(GPS_logs[key])   
        
    time_epoch_list = [date_obj.replace(tzinfo = pytz.timezone('UTC')).timestamp()
                     for date_obj in GPS_logs['UTC_t']]
    mean_epoch = np.mean(time_epoch_list)
    avg_logs['UTC_t'] = datetime.utcfromtimestamp(mean_epoch)
    
    return avg_logs

