#%% import libraries
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:19:51 2024

@author: moreaul
"""

import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pickle
import csv
#%% mettre les figures dans des fenetres
#import matplotlib
#matplotlib.use('TkAgg')
#%% code

plt.close('all')


date = '0210'
year = '2025'

#path2logfile = f'/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence pour ordi ludo
#geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

# windows:
#path2logfile = f'D:/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur ssd vasco
#geophones_table_path = 'D:/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur ssd vasco

#linux:
ordi = 'dell_vasco'

def routine(inputs,outputs,date,ordi='dell_vasco'):

    outputs[date] = {}

    #donnees exemple:

    if ordi=='adour':
#        if int(date)<=218:
#            path2logfile = f'/media/turbots/Shack25/Data/{date}/Geophones' # arborescence sur storageshared
#        elif int(date)>218:
#            path2logfile = f'/media/turbots/Shack25/Data/{date}/Geophones' # arborescence sur storageshared
        path2logfile = f'/media/turbots/Backup25/Data/{date}/Geophones'
        geophones_table_path = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared

    elif ordi=='dell_vasco':
    #    path2logfile = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur storageshared
    #    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
    #    path2logfile = f'D:/BicWin2025/Data/{date}/Geophones' # arborescence sur storageshared
        path2logfile = f'D:/copie_BicWin25_geophones/Data/{date}/Geophones'
        #path2logfile = f'B:/Data/{date}/Geophones'
        geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
        

    acquisition_numbers_to_plot = inputs[date]

    # Define the find_coordinates function
    def find_coordinates(file_path, date, year):
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
                        lat_avg = sum(latitudes) / len(latitudes) if latitudes else None
                        long_avg = sum(longitudes) / len(longitudes) if longitudes else None
                        GPS_coordinates_matrices.append({
                            "acquisition": num_snippets,
                            "num_GPS_logs": num_GPS_logs,
                            "avg_latitude": lat_avg,
                            "avg_longitude": long_avg
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
                lat_avg = sum(latitudes) / len(latitudes) if latitudes else None
                long_avg = sum(longitudes) / len(longitudes) if longitudes else None
                GPS_coordinates_matrices.append({
                    "acquisition": num_snippets,
                    "num_GPS_logs": num_GPS_logs,
                    "avg_latitude": lat_avg,
                    "avg_longitude": long_avg
                })

        return num_snippets, GPS_coordinates_matrices

    # Main function
    def main():
        log_files = [file for file in os.listdir(path2logfile) if file.startswith("DigiSolo") and file.endswith(".txt")]

        # Sort the log files based on DigiSolo number
        log_files.sort(key=lambda x: int(re.search(r'\d+', x).group()))

        all_GPS_coordinates_matrices = {}

        for log_file in log_files:
            log_file_path = os.path.join(path2logfile, log_file)
            num_snippets, GPS_coordinates_matrices = find_coordinates(log_file_path, date, '2024')
            acquisitions_dict = {entry['acquisition']: {
                'num_GPS_logs': entry['num_GPS_logs'],
                'avg_latitude': entry['avg_latitude'],
                'avg_longitude': entry['avg_longitude']
            } for entry in GPS_coordinates_matrices}

            all_GPS_coordinates_matrices[log_file] = acquisitions_dict

        return all_GPS_coordinates_matrices

    if __name__ == "__main__":
        all_matrices = main()

        # Print or use the matrices as needed
        for file, matrices in all_matrices.items():
            print(f"\nFile: {file}")
            for acquisition, matrix in matrices.items():
                print(f"Acquisition {acquisition}:")
                print(f"   GPS Coordinates Matrix: {matrix}")

    def plot_average_coordinates(all_matrices, acquisition_numbers, geophones_table_path):
        # Load geophones table
        geophones_table = {}
        with open(geophones_table_path, 'r') as table_file:
            # Skip the header
            next(table_file)
            for line in table_file:
                num, geo_sn = line.strip().split('\t')
                geophones_table[geo_sn] = num

        plt.figure(figsize=(10, 6))
        plt.title('Average GPS Coordinates for Acquisitions')

        for acquisition_number in acquisition_numbers:
            for file, matrices in all_matrices.items():
                if acquisition_number in matrices:
                    coordinates = matrices[acquisition_number]
                    if coordinates['avg_latitude'] is not None and coordinates['avg_longitude'] is not None:
                        avg_latitude = coordinates['avg_latitude']
                        avg_longitude = coordinates['avg_longitude']

                        # Extract the four-digit number from the filename using a correct regex pattern
                        file_number_match = re.search(r'DigiSolo_(\d{4})', file)
                        if file_number_match:
                            file_number = file_number_match.group(1)
                            geo_sn = "45302" + file_number.zfill(4)  # Add leading zeros and the prefix

                            # Look up the corresponding two-digit number from the geophones table
                            if geo_sn in geophones_table:
                                two_digit_number = geophones_table[geo_sn]
                                label = f"GN_{two_digit_number}"
                            else:
                                label = f"Unknown"

                            # Plot the average point for each file with the new label
                            plt.scatter(avg_latitude, avg_longitude, label=f'Acquisition {acquisition_number} - {label}')

        plt.xlabel('Latitude')
        plt.ylabel('Longitude')
        #plt.legend()
        
        output_file_path = os.path.join(path2logfile, f'average_coordinates_acqs_{acquisition_numbers}.png')
        plt.savefig(output_file_path)
        plt.show()

    # Specify the acquisition numbers you want to plot

    plot_average_coordinates(all_matrices, acquisition_numbers_to_plot, geophones_table_path)

    def save_average_coordinates(all_matrices, acquisition_number, path2logfile):
        # Prepare the full path to the text file
        output_file_name = os.path.join(path2logfile, f"GPS_coordinates_acq{acquisition_number}.txt")

        # Open the text file for writing
        with open(output_file_name, 'w') as output_file:
            # Write the header
            output_file.write("GN\tLatitude\tLongitude\n")

            # Write the average coordinates for each file
            for file, matrices in all_matrices.items():
                if acquisition_number in matrices:
                    coordinates = matrices[acquisition_number]
                    if coordinates['avg_latitude'] is not None and coordinates['avg_longitude'] is not None:
                        avg_latitude = coordinates['avg_latitude']
                        avg_longitude = coordinates['avg_longitude']

                        # Extract the file number from the filename (assuming it follows the DigiSolo_XX.LOG pattern)
                        file_number = int(re.search(r'\d+', file).group())

                        # Create the GN identifier (e.g., IGU_01)
                        gn_identifier = f"IGU_{file_number:02d}"

                        # Write the row in the text file
                        output_file.write(f"{gn_identifier}\t{avg_latitude}\t{avg_longitude}\n")
                        outputs[date]['acq'+str(acquisition_number)] = coordinates
        # Save the figure for this acquisition
        plot_average_coordinates(all_matrices, [acquisition_number], geophones_table_path)
        print(coordinates)
    # Save average coordinates and corresponding figures
    for acquisition_number_to_save in acquisition_numbers_to_plot:
        save_average_coordinates(all_matrices, acquisition_number_to_save, path2logfile)


########################################################################################

# %%
save = True
savefolder = 'B:/General/Summary_geophone_lines/'

inputs = {}
# anse à Mercier
inputs['0131'] = [1]
inputs['0201'] = [1]
inputs['0203'] = [1]

# baie du Ha! Ha!
inputs['0204'] = [1,2,3,4,5]
inputs['0206'] = [1]
inputs['0210'] = [1,2,3]
inputs['0212'] = [1,2]
inputs['0221'] = [1,2,3]
inputs['0224'] = [1,2]
inputs['0304'] = [1,2]

# baie de Rimouski
inputs['0227'] = [1,2]


Dates = ['0131','0201','0203','0204','0206','0210','0212','0221','0224','0227','0304']
#Dates = ['0204','0210','0224']


outputs = {}

for date in Dates:
    routine(inputs,outputs,date)

print(outputs)

plt.figure(figsize=(20,15))
count=0
for date in outputs:
    for acq_key in outputs[date]:    
        acq_num = int(acq_key[-1])
        if count<10:
            plt.plot(outputs[date][acq_key]['avg_longitude'],outputs[date][acq_key]['avg_latitude'],'o',label=date+', '+acq_key)
        elif count>=10:
            plt.plot(outputs[date][acq_key]['avg_longitude'],outputs[date][acq_key]['avg_latitude'],'s',label=date+', '+acq_key)


        count+=1

plt.legend()

# %%
# fonction simple qui enregistre un fichier csv avec 
# une premiere ligne qui contient les noms des colonnes et ensuites les colonnes

csv_file_path = savefolder + 'all_geophone_lines_coordinates.csv'

header = ['date', 'acq', 'longitude', 'latitude']

data = []
for d in Dates:
    print(d,type(d))
    for a in inputs[d]:
        data.append([year + d, a, outputs[d]['acq'+str(a)]['avg_longitude'], outputs[d]['acq'+str(a)]['avg_latitude']])

        #path2logfile = f'B:/Data/{d}/Geophones/'
        #t1_file = path2logfile + 't1_to_time_' + d + '_' + year + '.pkl'

with open(csv_file_path, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f,delimiter=',')

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(data)


#%%

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/vasco/field/')

from field.read_UTC_time_geophones import *

csv_file_path_times_UTC = savefolder + 'all_geophone_lines_times_UTC.csv'
header_times_UTC = ['date', 'acq', 'tstart', 'tend']


times_UTC = []

for d in Dates:
    print(d,type(d))
    for a in inputs[d]:
        tstart,tend = read_tstart_tend(year='2025',date=d,acqu_numb=str(a).zfill(4),ordi='dell_vasco')
        times_UTC.append([year + d, a,tstart,tend])

with open(csv_file_path_times_UTC, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f,delimiter=',')

    # write the header
    writer.writerow(header_times_UTC)

    # write multiple rows
    writer.writerows(times_UTC)


# %%
