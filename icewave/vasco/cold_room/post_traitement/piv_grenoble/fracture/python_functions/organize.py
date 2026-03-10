# import modules

import numpy as np
import h5py
import pickle
import csv

# functions

def mat_to_dict(mat_object,ref_matobj):
    """
    Recursively convert a MATLAB structure (HDF5 group or dataset) to a Python dictionary.
    

    INPUTS : - mat_object : matlab object extracted from a .mat file using h5py
             - ref_matobj : matlabo object of reference, main root of the matlab structure, used to dereference some values, 
                     needed for cell_arrays for instance 
                     
    OUTPUT : - whatever was in the .mat file : structure, substructures, cell_array, strings etc..
    
    """
    if isinstance(mat_object, h5py.Dataset):  # If it's a dataset, return its value
        data = mat_object[()]
                  
        # Handle MATLAB strings (stored as bytes)
        if data.dtype == 'uint16':  # Check if it's a string
            # Convert uint16 array to Python string (decode as Unicode characters)
            return ''.join(chr(code_point[0]) for code_point in data)

        # Handle case of cell array
        if data.dtype == 'O':
            new_data = np.empty(data.shape,dtype = object).ravel()
            for i,pointer in enumerate(data.flat):
                new_mat_object = ref_matobj[pointer]
                new_data[i] = mat_to_dict(new_mat_object,ref_matobj)
                
            new_data = new_data.reshape(data.shape)
            return new_data
    
    
        if isinstance(data, np.ndarray):
            data = np.squeeze(data) 
            if data.size == 1:  # If the array contains only one element, convert to float
                return float(data)
            else :
                return data
            
    
    elif isinstance(mat_object, h5py.Group):  # If it's a group (structure), create a dictionary
        result_dict = {}
        for key, item in mat_object.items():
            result_dict[key] = mat_to_dict(item,ref_matobj)  # Recursively call for each element
        return result_dict
    
    else:
        raise TypeError(f"Unsupported type {type(mat_object)}")

def matcell2dict_PIV(matcell,dim_keys = 0):
    """ Create a dictionnary for a 2xN matlab cell array 
    
    INPUT : - matcell, cell array converted as an array using h5py and function mat_to_dict, 
                row 0 -> keys and row 1 -> values 
                
    OUTPUT : - a python dictionnary whose keys correspond to keys stored in the first dimension 
    
    """
    my_dict = {}
    keys = matcell[dim_keys]
    for i,key in enumerate(keys) :
        my_dict[key] = matcell[1,i]
        
    return my_dict


def savedict(filename,namedict):
    with open(filename+'.pickle', 'wb') as handle:
        pickle.dump(namedict, handle, protocol=pickle.HIGHEST_PROTOCOL)

def opendict(filename):
    with open(filename+'.pickle', 'rb') as handle:
        dico = pickle.load(handle)
    return dico

def csv2dict(csv_file_path=str):
    # Initialize the dictionary
    input_params = {}
    # Read the CSV file
    with open(csv_file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)  # Assuming the file is tab-delimited
        for row in reader:
            if len(row) == 3:  # Ensure the row has two elements (key and value)
                key, typ, value = row
                print(row)
                # Try to convert the value to a float or int if possible
                try:
                    if ',' in value:  # Replace commas with dots for decimal values
                        value = value.replace(',', '.')
                    value = float(value) if '.' in value else int(value)
                except ValueError:
                    pass  # Keep the value as a string if conversion fails
                if typ == 'int':
                    input_params[key] = int(value)
                elif typ == 'float':
                    input_params[key] = float(value)
                elif typ == 'str':
                    input_params[key] = str(value)
                elif typ == 'bool':
                    input_params[key] = bool(value)
    return input_params            
