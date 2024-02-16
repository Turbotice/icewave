import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pickle

path2logfile = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0206/Geophones'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'


pickle_file_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0206/Geometry/0001_Y_m_vs_X_mgeophone_line_16.pkl'

# Load the pickled data
with open(pickle_file_path, 'rb') as file:
    loaded_data = pickle.load(file)

# Now, 'loaded_data' contains the data from the pickled file
print(loaded_data)




def find_coordinates(log_content):
    latitude_pattern = re.compile(r'Latitude\s+=\s+([-+]?\d*\.\d+|\d+)')
    longitude_pattern = re.compile(r'Longitude\s+=\s+([-+]?\d*\.\d+|\d+)')

    start_acquisition_pattern = re.compile(r'Start Acquisition FileName')

    start_positions = [match.end() for match in start_acquisition_pattern.finditer(log_content)]

    GPS_coordinates_matrices = {}

    for i, start_position in enumerate(start_positions):
        end_position = start_positions[i + 1] if i + 1 < len(start_positions) else None
        snippet = log_content[start_position:end_position]

        latitude_values = latitude_pattern.findall(snippet)
        longitude_values = longitude_pattern.findall(snippet)

        # Append coordinates to the matrix
        coordinates_acquisition = list(zip(latitude_values, longitude_values))
        GPS_coordinates_matrices[i + 1] = coordinates_acquisition

    return GPS_coordinates_matrices

def main():
    log_files = [file for file in os.listdir(path2logfile) if file.startswith("DigiSolo") and file.endswith(".txt")]

    # Sort the log files based on DigiSolo number
    log_files.sort(key=lambda x: int(re.search(r'\d+', x).group()))

    all_GPS_coordinates_matrices = {}

    for log_file in log_files:
        log_file_path = os.path.join(path2logfile, log_file)
        
        with open(log_file_path, 'r') as file:
            log_content = file.read()
            GPS_coordinates_matrices = find_coordinates(log_content)

            all_GPS_coordinates_matrices[log_file] = GPS_coordinates_matrices

    return all_GPS_coordinates_matrices

if __name__ == "__main__":
    all_matrices = main()

    # Print or use the matrices as needed
    for file, matrices in all_matrices.items():
        print(f"\nFile: {file}")
        for acquisition, matrix in matrices.items():
            print(f"Acquisition {acquisition}:")
            print(f"   GPS Coordinates Matrix: {matrix}")

#matrix_1_1 = all_matrices['DigiSolo_07.LOG'][1]
#print(f"GPS Coordinates Matrix for File DigiSolo_08.LOG, Acquisition 1: {matrix_1_1}")



# Assuming 'all_matrices' is the dictionary returned from your main function
# (the one you printed in the original script)
# It contains the GPS coordinates matrices for each file and acquisition

def plot_average_coordinates(all_matrices, acquisition_number, geophones_table_path):
    # Load geophones table
    geophones_table = {}
    with open(geophones_table_path, 'r') as table_file:
        # Skip the header
        next(table_file)
        for line in table_file:
            num, geo_sn = line.strip().split('\t')
            geophones_table[geo_sn] = num

    plt.figure(figsize=(10, 6))
    plt.title(f'Average GPS Coordinates for Acquisition {acquisition_number}')

    for file, matrices in all_matrices.items():
        if acquisition_number in matrices:
            coordinates = matrices[acquisition_number]
            coordinates_np = np.array(coordinates, dtype=float)

            if coordinates_np.ndim == 1:
                avg_latitude = coordinates_np[0]
                avg_longitude = coordinates_np[1]
            else:
                avg_latitude = np.mean(coordinates_np[:, 0])
                avg_longitude = np.mean(coordinates_np[:, 1])

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
                plt.scatter(avg_latitude, avg_longitude, label=label)

    plt.xlabel('Latitude')
    plt.ylabel('Longitude')
    plt.legend()
    
    output_file_path = os.path.join(path2logfile, f'average_coordinates_acq{acquisition_number}.png')
    plt.savefig(output_file_path)
    plt.show()

acquisition_number = 2
plot_average_coordinates(all_matrices, acquisition_number, geophones_table_path)


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
                coordinates_np = np.array(coordinates, dtype=float)

                if coordinates_np.ndim == 1:
                    avg_latitude = coordinates_np[0]
                    avg_longitude = coordinates_np[1]
                else:
                    avg_latitude = np.mean(coordinates_np[:, 0])
                    avg_longitude = np.mean(coordinates_np[:, 1])

                # Extract the file number from the filename (assuming it follows the DigiSolo_XX.LOG pattern)
                file_number = int(re.search(r'\d+', file).group())

                # Create the GN identifier (e.g., IGU_01)
                gn_identifier = f"IGU_{file_number:02d}"

                # Write the row in the text file
                output_file.write(f"{gn_identifier}\t{avg_latitude}\t{avg_longitude}\n")




save_average_coordinates(all_matrices, acquisition_number+1, path2logfile)




