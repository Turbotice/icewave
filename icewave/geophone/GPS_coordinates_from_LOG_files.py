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

plt.close('all')

date = '0221'
year = '2024'
path2logfile = f'/YOUR_DATA_PATH_HERE/{year}_BICWIN/{date}/Geophones'
geophones_table_path = '/icewave/geophone/geophones_table'
acquisition_numbers_to_plot = [1,2,3]  # Example: plot acquisitions 1, 2, and 3

# Define the find_coordinates function
def find_coordinates(file_path, date, year):
    num_snippets = 0
    target_date_found = False
    target_date_pattern = year + '/' + date[:2] + '/' + date[2:]
    GPS_coordinates_matrices = []
    current_snippet = ""

    with open(file_path, 'r') as file:
        for line in file:
            # Detect start of acquisition snippet
            if line.startswith("Start Acquisition FileName = "):
                #print(line)
                
                if current_snippet:
                    # Process the previous snippet before starting a new one
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
                current_snippet += line
            # Accumulate snippet lines once target date found or within snippet
            elif current_snippet:
                current_snippet += line

            # Check for target date to start reading snippets
            elif line.startswith("UTC Time ="):
                match = re.search(r'UTC Time = "(\d{4}/\d{2}/\d{2})', line)
                if match and match.group(1) >= target_date_pattern:
                    target_date_found = True
                    # Start capturing from now
                    current_snippet = ""

        # Process last snippet if file ends while capturing
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
        num_snippets, GPS_coordinates_matrices = find_coordinates(log_file_path, date, year)
        
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
    # for file, matrices in all_matrices.items():
    #     print(f"\nFile: {file}")
    #     for acquisition, matrix in matrices.items():
    #         print(f"Acquisition {acquisition}:")
    #         print(f"   GPS Coordinates Matrix: {matrix}")

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

    # Define a colormap
    colors = plt.cm.get_cmap('tab10', len(acquisition_numbers))
    for idx, acquisition_number in enumerate(acquisition_numbers):
        color = colors(idx)
        for file, matrices in all_matrices.items():
            if acquisition_number in matrices:
                coordinates = matrices[acquisition_number]
                print(coordinates)
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

                        # Plot the average point for each file with the new label and color
                        plt.scatter(avg_longitude,avg_latitude, label=f'Acquisition {acquisition_number} - {label}', color=color)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    #plt.legend()
    
    output_file_path = os.path.join(path2logfile, f'average_coordinates_acqs_{acquisition_numbers}.png')
    plt.savefig(output_file_path, format='png')
    plt.show()

# Specify the acquisition numbers you want to plot
plot_average_coordinates(all_matrices, acquisition_numbers_to_plot, geophones_table_path)

# def save_average_coordinates(all_matrices, acquisition_number, path2logfile):
#     # Prepare the full path to the text file
#     output_file_name = os.path.join(path2logfile, f"GPS_coordinates_acq{acquisition_number}.txt")

#     # Open the text file for writing
#     with open(output_file_name, 'w') as output_file:
#         # Write the header
#         output_file.write("GN\tLatitude\tLongitude\n")

#         # Write the average coordinates for each file
#         for file, matrices in all_matrices.items():
#             if acquisition_number in matrices:
#                 coordinates = matrices[acquisition_number]
#                 if coordinates['avg_latitude'] is not None and coordinates['avg_longitude'] is not None:
#                     avg_latitude = coordinates['avg_latitude']
#                     avg_longitude = coordinates['avg_longitude']

#                     # Extract the file number from the filename (assuming it follows the DigiSolo_XX.LOG pattern)
#                     file_number = int(re.search(r'\d+', file).group())

#                     # Create the GN identifier (e.g., IGU_01)
#                     gn_identifier = f"IGU_{file_number:02d}"

#                     # Write the row in the text file
#                     output_file.write(f"{gn_identifier}\t{avg_latitude}\t{avg_longitude}\n")

#     # Save the figure for this acquisition
#     plot_average_coordinates(all_matrices, [acquisition_number], geophones_table_path)

# # Save average coordinates and corresponding figures
# for acquisition_number_to_save in acquisition_numbers_to_plot:
#     save_average_coordinates(all_matrices, acquisition_number_to_save, path2logfile)
