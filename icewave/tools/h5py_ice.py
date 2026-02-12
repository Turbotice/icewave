




import h5py
import numpy as np

def write(filename,R):
#Writing data
    hf = h5py.File(filename, "w")
    for key in R.keys():
        tables={}
        headers={}
        for k in R[key]:
            if type(R[key][k]) is list:
                #go in a table
                tables[k] = np.asarray(R[key][k])
            else:
                headers[k]=R[key][k]
        grp = hf.create_group(str(key))

        for k in tables.keys():
            grp.create_dataset(str(k), data=tables[k])
        for k in headers.keys():
            grp.attrs[k] = headers[k] #all the same values of headers.
    hf.close()

#Reading data

def read(filename):
    hf1 = h5py.File(filename, "r")
    for name in hf1:
        dset = hf1[name]
        print(name)
        
    for key in hf1.attrs.keys():
        print(key)
        print(hf1.attrs[key])
    #hf1.close()
    return hf1
