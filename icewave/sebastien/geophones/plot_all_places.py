#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 22:02:25 2024

@author: moreaul
"""
%matplotlib qt
import csv
import numpy as np
import matplotlib.pyplot as plt
import gmplot

# Define the base file path
path_data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/'

# List of dates to process
dates = ['0206','0210','0211','0221','0223','0226','0306']  # Add more dates as needed

# Acquisition number
acqu = '1'

# Initialize a dictionary to hold the mean latitudes and longitudes for each date
mean_coordinates = {}

for date in dates:
    GPS_file = f"{path_data}{date}/Geophones/GPS_coordinates_acq{acqu}.txt"

    try:
        # Initialize lists to hold latitudes and longitudes
        latitudes = []
        longitudes = []

        # Open the file and read its contents
        with open(GPS_file, 'r') as file:
            # Read the header line
            header = file.readline()

            # Process each line in the file
            for line in file:
                # Split the line by whitespace
                parts = line.split()

                # Extract latitude and longitude values
                latitude = float(parts[1])
                longitude = float(parts[2])

                # Append values to lists
                latitudes.append(latitude)
                longitudes.append(longitude)

        if latitudes and longitudes:
            # Calculate mean latitude and longitude
            lat_mean = np.mean(latitudes)
            lon_mean = np.mean(longitudes)

            # Store the means in the dictionary
            mean_coordinates[date] = {'lat_mean': lat_mean, 'lon_mean': lon_mean}
        else:
            print(f"No data found in {GPS_file}")

    except FileNotFoundError:
        print(f"File not found: {GPS_file}")
    except Exception as e:
        print(f"An error occurred while processing {GPS_file}: {e}")

# Check if we have data to plot
if mean_coordinates:
    # Extract the lat_mean and lon_mean for plotting
    lat_means = [coords['lat_mean'] for coords in mean_coordinates.values()]
    lon_means = [coords['lon_mean'] for coords in mean_coordinates.values()]

    # Plot lat_mean vs lon_mean
    plt.figure(figsize=(10, 6))
    plt.scatter(lon_means, lat_means, color='blue', marker='o')

    # Add labels and title
    for date, lat, lon in zip(dates, lat_means, lon_means):
        plt.text(lon, lat, date, fontsize=12, ha='right')
    plt.xlabel('Longitude Mean')
    plt.ylabel('Latitude Mean')
    plt.title('Mean Latitude vs Mean Longitude for Different Dates')
    plt.grid(True)
    plt.show()
else:
    print("No valid data to plot.")





# Check if we have data to save
if mean_coordinates:
    # Create a CSV file to save the coordinates
    csv_file = "gps_coordinates.csv"
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Date', 'Latitude', 'Longitude'])  # Write header
        for date, coords in mean_coordinates.items():
            writer.writerow([date, coords['lat_mean'], coords['lon_mean']])
    print(f"GPS coordinates saved to {csv_file}")
else:
    print("No valid data to save.")

