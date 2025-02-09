#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 18:35:07 2024

@author: moreaul
"""
import os
import shutil

#%% 
#####################################################################
## Change miniseed files name, using geophone correspondance table
#####################################################################

#path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0206/Geophones'
path2data = 'E:Data/0206/Geophones/'
geophones_table_path = 'D:/PC Seb/These PMMH/Rimouski_2024/Geophones/geophones_table'


# Read geophones_table into a dictionary
geophones_dict = {}
with open(geophones_table_path, 'r') as file:
    for line in file:
        num, geo_SN = line.strip().split()
        geophones_dict[geo_SN] = num

# Specify the directory where the MiniSEED files are located
miniseed_directory = path2data

# Iterate through each MiniSEED file in the directory
for filename in os.listdir(miniseed_directory):
    if filename.endswith(".miniseed"):
        # Extract the geophone serial number from the filename
        geophone_serial_number = filename[:9]
        
        # Check if the geophone serial number is in the geophones_table
        if geophone_serial_number in geophones_dict:
            # Replace the first 9 characters with the associated 2-digit number
            new_filename = geophones_dict[geophone_serial_number] + filename[9:]
            
            # Rename the file
            old_path = os.path.join(miniseed_directory, filename)
            new_path = os.path.join(miniseed_directory, new_filename)
            os.rename(old_path, new_path)
            print(f'Renamed: {filename} -> {new_filename}')
        else:
            print(f'Geophone serial number {geophone_serial_number} not found in geophones_table.')
            
#%% 
####################################
## Sort miniseed files into folders
####################################

#path2data = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Rimouski_2024/Data/0210/geophones/'

for filename in os.listdir(path2data):
    if filename.endswith(".miniseed"):
        new_filename = filename.replace('..', '.', 2)
        full_path_old = os.path.join(path2data, filename)
        full_path_new = os.path.join(path2data, new_filename)
        os.rename(full_path_old, full_path_new)


files = os.listdir(path2data)

# Iterate through each file
for filename in files:
    print(filename[:3 ])
    # Check if the file has the specified pattern
    if filename[:2].isdigit() and filename[2] == '.' and filename[3].isdigit() and filename[4] == '.':
        # Replace the second dot with two zeros
        new_filename = filename[:4 ] + '00' + filename[5:]
        
        # Construct the full paths
        old_path = os.path.join(path2data, filename)
        new_path = os.path.join(path2data, new_filename)
        
        # Rename the file
        os.rename(old_path, new_path)
        print(f'Renamed: {filename} to {new_filename}')


# Create a dictionary to store files based on acq_number
acq_number_dict = {}

# Iterate through files in the folder
for filename in os.listdir(path2data):
    if filename.endswith(".miniseed"):  # Adjust the file extension as needed
        # Extract acq_number from the filename
        parts = filename.split('.')
        if len(parts) >= 3:
            acq_number = parts[1]
            
            # Create a folder for the acq_number if it doesn't exist
            acq_folder = os.path.join(path2data, acq_number)
            os.makedirs(acq_folder, exist_ok=True)
            
            # Move the file to the corresponding folder
            source_path = os.path.join(path2data, filename)
            destination_path = os.path.join(acq_folder, filename)
            shutil.move(source_path, destination_path)




#%% 
###########################################################################
## Rename miniseed files using .LOG files and geophone correspondance table
###########################################################################

def find_serial_number(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if 'Serial Number' in line:
                return line.split()[-1]

def find_2_digit_number(geo_table_path, serial_number):
    with open(geo_table_path, 'r') as geo_table:
        for line in geo_table:
            parts = line.split()
            if len(parts) == 2 and parts[1] == serial_number:
                return parts[0]
            

def process_files(files_directory, geo_table_path):
    for filename in os.listdir(files_directory):
        if filename.endswith(".LOG"):
            file_path = os.path.join(files_directory, filename)
            serial_number = find_serial_number(file_path)
            if serial_number:
                two_digit_number = find_2_digit_number(geo_table_path, serial_number)
                if two_digit_number:
                    # Extract the 4-digit number from the filename
                    original_number = filename.split('_')[1].split('.')[0]

                    # Replace the 4-digit number with the two-digit number
                    new_filename = filename.replace(original_number, two_digit_number)
                    
                    # Rename the file with the new filename
                    new_file_path = os.path.join(files_directory, new_filename)
                    os.rename(file_path, new_file_path)

                    print(f"Renamed {filename} to {new_filename}")


process_files(path2data, geophones_table_path)
            
            