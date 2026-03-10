#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.io import loadmat
import pickle
from scipy.optimize import curve_fit

# %%
# open 100 images, named from Basler_acA2040-90um__23029848__20260127_182642777_0420
# to Basler_acA2040-90um__23029848__20260127_182642777_0429, and store their pixel values in numpy arrays
# (use a 3d numpy array to store the pixel values of all 10 images, with shape (10, height, width))

disk = 'D:'
img_dir = f"{disk}/manips_BsAs/Basler_images/20260127/acq2_expo2000_facq50/cam2/"
images = []
for i in range(100):
    filename = f"{img_dir}Basler_acA2040-90um__23029848__20260127_182642777_{str(350+i).zfill(4)}.tiff"
    image = plt.imread(filename)
    images.append(image)
images_array = np.array(images)
# %%
%matplotlib qt
plt.figure()
plt.imshow(images_array[0],cmap='gray')
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(images_array[0,1035,:750])
plt.show()

plt.figure(figsize=(10, 10))
plt.plot(np.gradient(images_array[0,1035,:750]))
plt.show()
#%%

def find_peak_near_x0(signal, x0, window_size=10,subpix=True):
    """
    Find the peak in a 1D signal near a given x0, within a specified window size.

    Parameters:
    signal (numpy array): The 1D signal to search for peaks.
    x0 (int): The x-coordinate around which to search for peaks.
    window_size (int): The size of the window around x0 to search for peaks.

    Returns:
    int: The x-coordinate of the peak found near x0, or None if no peak is found.
    """
    # Define the search range
    start = max(0, x0 - window_size)
    end = min(len(signal), x0 + window_size + 1)

    # Extract the segment of the signal to search for peaks
    segment = signal[start:end]

    # Find the index of the maximum value in the segment
    if len(segment) == 0:
        return None  # No data to search

    peak_index = np.argmax(np.abs(segment))

    # Convert the local index back to the original signal index
    peak_x = start + peak_index

    # subpixel peak detection using quadratic interpolation
    if (peak_index > 0 and peak_index < len(segment) - 1)&(subpix==True):
        x1 = peak_x - 1
        y1 = signal[x1]
        x2 = peak_x
        y2 = signal[x2]
        x3 = peak_x + 1
        y3 = signal[x3]

        A = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
        B = np.array([y1, y2, y3])
        Coefficients = np.linalg.solve(A, B)

        x_max_subpix = -Coefficients[1] / (2 * Coefficients[0])
        return x_max_subpix
    elif (peak_index == 0)|(peak_index == len(segment) - 1):
        print("Warning: Peak is at the edge of the search window, subpixel detection may not be accurate.")
        return peak_index
    elif subpix==False:
        return peak_x

# test the function on a simple signal
# signal is taken from the first image, using the row at y=1035 and looking for a peak near x0=750
y = 850
x0 = 1972
signal = np.gradient(images_array[2, y, :])
peak_x = find_peak_near_x0(signal, x0,window_size=10)
# plot the signal and the peak found
plt.figure(figsize=(10, 10))
plt.plot(signal)
plt.axvline(x=peak_x, color='r', linestyle='--')
plt.title(f"Signal and peak found at x={peak_x}")
plt.show()
#%%
# test the function on the first image, using the row at y=1035 and looking for a peak near x0=750
y = 1035
x0 = 198
signal = np.gradient(images_array[0, y, :])
peak_x = find_peak_near_x0(signal, x0)

peaks_arr = []
for i in range(len(images_array[:,0,0])):
    signal = np.gradient(images_array[i, y, :])
    peak_x = find_peak_near_x0(signal, x0)
    peaks_arr.append(peak_x)
peaks_arr = np.array(peaks_arr) - peaks_arr[0]
plt.figure(figsize=(10, 10))
plt.plot(peaks_arr)
plt.show()

# now take a peak_arr_ref
y_ref = 27
x0_ref = 894

peaks_arr_ref = []

for i in range(len(images_array[:,0,0])):
    signal_ref = np.gradient(images_array[i, y_ref, :])
    peak_x_ref = find_peak_near_x0(signal_ref, x0_ref)
    peaks_arr_ref.append(peak_x_ref)
peaks_arr_ref = np.array(peaks_arr_ref) - peaks_arr_ref[0]
plt.figure(figsize=(10, 10))
plt.plot(peaks_arr)
plt.plot(peaks_arr_ref )
plt.plot(peaks_arr - peaks_arr_ref, label='relative motion')
plt.legend()
plt.show()
# %%
# now import piv data for the reference

disk = 'D:'
csv_file_path = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/params_acquisitions/params_allacq.csv"

with open(csv_file_path, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=';')
    count=0
    data = []
    for lines in csvFile:
        if count==0:
            header = lines
        else:
            data.append(lines)
        count+=1
print(header)

data_csv = np.array(data)

# --- Convert all string entries to float ---
def str_to_float(value):
    try:
        # Replace comma with dot (European decimal format)
        value = value.replace(',', '.')
        # Remove spaces just in case
        value = value.strip()
        # Convert to float
        return float(value)
    except ValueError:
        if (value=='nan')|(value=='None')|(value=='#VALEUR !'):
            return np.nan
        else:
            return value  # Return the original value if it cannot be converted
# Apply the conversion to every element
#data_csv = np.vectorize(str_to_float)(data_csv)

data_csv_converted = np.zeros(data_csv.shape, dtype=object)  # Create an empty array with the same shape
for i in range(data_csv.shape[0]):
    for j in range(data_csv.shape[1]):
        data_csv_converted[i, j] = str_to_float(data_csv[i, j])
data_csv = data_csv_converted
#print(data_csv_converted)

# displacement vs frames with piv
idx = 40

use_ref = int(data_csv[idx, 20])

frame_frac = data_csv[idx, 12]
frame_force_application = data_csv[idx,11] #début de l'application de la force

frame_ref_cam1 = data_csv[idx, 16]
frame_ref_cam2 = data_csv[idx, 17]

date = str(int(data_csv[idx, 0])).zfill(4)
acq = int(data_csv[idx, 1])
serie = data_csv[idx, 2] # si pas de serie, mettre serie = None

if np.isnan(serie):
    serie = None
else:
    serie = int(serie)

if np.isnan(frame_ref_cam2):
    pass
else:
    frameref_diff = frame_ref_cam2 - frame_ref_cam1
    frame_frac += frameref_diff

i0 = int(data_csv[idx, 7])
refimg = int(data_csv[idx, 8])
W = int(data_csv[idx, 9])
n_passages = int(data_csv[idx, 10])


if type(serie)==int:
    path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}/serie{serie}'
else:
    path2dailydir = f'{disk}/manips_BsAs/Basler_images/{date}'

lisdir = os.listdir(path2dailydir)
for acqname in lisdir:
    if f'acq{acq}_' in acqname:
        acqname_found = acqname

matfile = f'{path2dailydir}/{acqname_found}/PIV_processed_{n_passages}passages_i0{i0}_refimg{refimg}_W{W}.mat'

mat_dict = loadmat(matfile)

# deux facons d'importer les données en fonction de la structure du fichier .mat
def reshape_array(arr):
    array_new = np.zeros((len(arr),arr[0].shape[0],arr[0].shape[1]))
    for i in range(len(arr)):
        array_new[i,:,:] = arr[i]
    return array_new

def clean_array(arr):
    # On va parcourir le tableau et garder seulement les sous-tableaux qui ont au moins un élément non nul ou non NaN
    cleaned = []

    for sublist in arr:
        for subarr in sublist:
            # Vérifie si le tableau contient au moins un élément significatif
            if np.any(~np.isnan(subarr) & (subarr != 0)):
                cleaned.append(subarr)
    # Convertir en tableau numpy propre
    cleaned_array = np.array(cleaned, dtype=object)
    return cleaned_array

if 'u_original' in mat_dict:
    u_original = mat_dict['u_original'][:,0]
    v_original = mat_dict['v_original'][:,0]
    u = reshape_array(u_original)
    v = reshape_array(v_original)
    xpix = mat_dict['x'][0][0][0]
    ypix = mat_dict['y'][0][0][:,0]

    concat_with_zero = True
else:
    #u = np.moveaxis(clean_array(mat_dict['u']), 2, 1)
    #v = np.moveaxis(clean_array(mat_dict['v']), 2, 1)
    u = clean_array(mat_dict['u'])
    v = clean_array(mat_dict['v'])

    xpix = mat_dict['xpix'].flatten()
    ypix = mat_dict['ypix'].flatten()
    concat_with_zero = False


# pour ajouter un zero au début (dans le cas imgref et piv faite à la mains, car  sinon ca décale par rapport aux mesures de forces)
if concat_with_zero:
    u = np.concatenate((np.zeros((1,u.shape[1],u.shape[2])), u))
    v = np.concatenate((np.zeros((1,v.shape[1],v.shape[2])), v))


# tracking using coordonnées des vis

xpx_vis = int(data_csv[idx, 3])
ypx_vis = int(data_csv[idx, 4])
xpx_vis_ref = int(data_csv[idx, 5])
ypx_vis_ref = int(data_csv[idx, 6])

xind = (np.abs(xpix - xpx_vis)).argmin()
yind = (np.abs(ypix - ypx_vis)).argmin()

xind_ref = (np.abs(xpix - xpx_vis_ref)).argmin()
yind_ref = (np.abs(ypix - ypx_vis_ref)).argmin()


dmm_sur_dpx = float(data_csv[idx, 14])
dmm = dmm_sur_dpx
dpx = 1

nb_neighbors2avg = 0

# average neighbours and smooth vs time :
#u_field = gaussian_filter1d(np.mean(u[:,yind-2:yind+2,xind-2:xind+2],axis=(1,2)),5)
if nb_neighbors2avg!=0:
    u_field = np.nanmean(u[:,yind-nb_neighbors2avg:yind+nb_neighbors2avg,xind-nb_neighbors2avg:xind+nb_neighbors2avg],axis=(1,2))
    u_field_ref = np.nanmean(u[:,yind_ref-nb_neighbors2avg:yind_ref+nb_neighbors2avg,xind_ref-nb_neighbors2avg:xind_ref+nb_neighbors2avg],axis=(1,2))
else:
    u_field = u[:,yind,xind]
    u_field_ref = u[:,yind_ref,xind_ref]

if use_ref:
    u_field_relative = u_field - u_field_ref
else:
    u_field_relative = u_field

plt.figure()
plt.plot(np.arange(u.shape[0])+i0,u_field, label='moving point')
plt.plot(np.arange(u.shape[0])+i0,u_field_ref, label='ref (point immobile)')
plt.plot(np.arange(u.shape[0])+i0,u_field_relative, label='relative motion')
# add the points from peaks_arr, but only for the frames where we have piv data (between i0 and i0+u.shape[0])
frames_with_piv = np.arange(i0, i0+u.shape[0])
plt.plot(np.arange(len(peaks_arr)) + i0, peaks_arr - peaks_arr[0], label='peaks')
plt.xlim(i0,frame_frac)
#plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
plt.legend()
plt.show()
# %%

def track_peak_and_plot(images_array, y, x0, y_ref=0, x0_ref=0,ref=False):
    peaks_arr = []
    for i in range(len(images_array[:,0,0])):
        signal = np.gradient(images_array[i, y, :])
        peak_x = find_peak_near_x0(signal, x0,window_size=30)
        peaks_arr.append(peak_x)
    peaks_arr = np.array(peaks_arr) - peaks_arr[0]
    if ref:
        peaks_arr_ref = []
        for i in range(len(images_array[:,0,0])):
            signal_ref = np.gradient(images_array[i, y_ref, :])
            peak_x_ref = find_peak_near_x0(signal_ref, x0_ref)
            peaks_arr_ref.append(peak_x_ref)
        peaks_arr_ref = np.array(peaks_arr_ref) - peaks_arr_ref[0]
        plt.figure(figsize=(10, 10))
        plt.plot(peaks_arr)
        plt.plot(peaks_arr_ref )
        plt.plot(peaks_arr - peaks_arr_ref, label='relative motion')
        plt.legend()
        plt.show()
        return peaks_arr, peaks_arr_ref
    else:
        plt.figure(figsize=(10, 10))
        plt.plot(peaks_arr)
        plt.show()
        return peaks_arr

track_peak_and_plot(images_array, y=825, x0=560, ref=False)
