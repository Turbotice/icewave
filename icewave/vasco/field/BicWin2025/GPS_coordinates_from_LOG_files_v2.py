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
#donnees exemple:

if ordi=='adour':
    path2logfile = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur storageshared
    geophones_table_path = f'/media/turbots/DATA/thiou/storageshared/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared

elif ordi=='dell_vasco':
#    path2logfile = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/{year}_BICWIN/{date}/Geophones' # arborescence sur storageshared
#    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared
    path2logfile = f'D:/BicWin2025/Data/{date}/Geophones' # arborescence sur storageshared
    geophones_table_path = f'X:/Banquise/Vasco/Startup_kit_Stage_MSIM/data/geophones_table' # arborescence sur storageshared


acquisition_numbers_to_plot = [2]#[1,2,3,4]

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

    # Save the figure for this acquisition
    plot_average_coordinates(all_matrices, [acquisition_number], geophones_table_path)

# Save average coordinates and corresponding figures
for acquisition_number_to_save in acquisition_numbers_to_plot:
    save_average_coordinates(all_matrices, acquisition_number_to_save, path2logfile)

########################################################################################

# %%
