
import numpy as np

"""def correspond_samplenum_acqnum(acq=1,serie=1):
    if serie!=None:
        fd = ((serie-1) * 3 + 1) + acq//3
    else:
        fd = 1+(acq)//3
    sd = acq%3
    if sd==0:
        sd=3
    if fd==0:
        fd = 1
    return fd, sd
"""
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