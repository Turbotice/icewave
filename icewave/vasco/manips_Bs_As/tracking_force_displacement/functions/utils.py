import numpy as np
from scipy.io import loadmat
import csv


disk_example = 'D:'


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

def load_csv_params_allacq(csv_file_path = f"{disk_example}/manips_BsAs/Summary/tracking_force_displacement/params_acquisitions/params_allacq.csv"):
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
    return data_csv

def correspond_samplenum_acqnum(date='1111', acq=1,serie=1, disk='D:'):
    if serie==None:
        q = (acq - 1) //3
        r = (acq - 1) % 3
        sample_num_found = str(int(10 * (q + 1) + (r + 1)))
        fd = sample_num_found[0]
        sd = sample_num_found[1]        
    else:
        infos_series_path = f'{disk}/manips_BsAs/epaisseurs/{date}/infos_series.txt'
        data_infos_series = np.loadtxt(infos_series_path, skiprows=1)
        sample_nums = data_infos_series[:,0]
        series = data_infos_series[:,1]
        indices_this_serie = np.where(series==serie)[0]
        # acq sera le 3eme indicede ce tableau "indices_this_serie"
        # donc le sample num sera :
        sample_num_found = str(int(sample_nums[indices_this_serie[acq-1]]))
        fd = sample_num_found[0]
        sd = sample_num_found[1]
    return fd, sd

######################################
# fonctions utilisées pour charger les fichier .mat renvoyés par PIVlab
# iIl y en a 2 sortes : ceux faits lors du traitement à la main, et ceux 
# faits lors du traitement automatisé. Leurs structures sont un peu différentes

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

def load_displacement_data_PIV(matfile):


    mat_dict = loadmat(matfile)

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

    return u, v, xpix, ypix, concat_with_zero

def load_data_PIV_without_coord(matfile):

    mat_dict = loadmat(matfile)

    if 'u_original' in mat_dict:
        u_original = mat_dict['u_original'][:,0]
        v_original = mat_dict['v_original'][:,0]
        u = reshape_array(u_original)
        v = reshape_array(v_original)

        concat_with_zero = True
    else:
        u = clean_array(mat_dict['u'])
        v = clean_array(mat_dict['v'])

        concat_with_zero = False


    # pour ajouter un zero au début (dans le cas imgref et piv faite à la mains, car  sinon ca décale par rapport aux mesures de forces)
    if concat_with_zero:
        u = np.concatenate((np.zeros((1,u.shape[1],u.shape[2])), u))
        v = np.concatenate((np.zeros((1,v.shape[1],v.shape[2])), v))

    return u, v, concat_with_zero