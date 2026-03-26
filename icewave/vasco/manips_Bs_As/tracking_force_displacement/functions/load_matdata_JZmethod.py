#%%
import numpy as np
import matplotlib.pyplot as plt
import h5py
#%%

def load_displacement_px(path2matfile = "D:/manips_BsAs/Summary/tracking_force_displacement/tracking_method_JZ/vasco_exemple_loc.mat"):
    with h5py.File(path2matfile, 'r') as f:
        mat_dict = {key: f[key][()] for key in f.keys()}

    loc_new = mat_dict['loc'].flatten()
    print(mat_dict['loc'].flatten())
    print(len(loc_new))
    return loc_new