
import icewave.phone.pipeline as pl

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df

import icewave.field.multi_instruments as multi
import icewave.field.time as timest

import numpy as np

def load_data(date,phonelist,nums,tmin,tmax,orientation):
    tmin = multi.convert_time(tmin)
    tmax = multi.convert_time(tmax)

    data = {}
    names = ['Z','E','N']
    for phone in phonelist:
        print(f'Load data N1 for phone {phone}')
        data[phone]={}
        ts = []
        for c,name in zip(orientation,names):
            data[phone]['a_'+name]=[]
        
        for num in nums:
            hf = pl.load_lvl_1(date,phone,num)
            t = list(hf['t_sync'])
            t = timest.today_time(t)
            ts = ts+t
            for c,name in zip(orientation,names):
                data[phone]['a_'+name] = data[phone]['a_'+name]+list(hf['a'+c])

        ts = np.asarray(ts)
        indices = np.logical_and(ts>=tmin,ts<tmax)
        data[phone]['ts'] =  ts[indices]

        for c,name in zip(orientation,names):
            data[phone]['a_'+name] = np.asarray(data[phone]['a_'+name])[indices]
    return data

def save_S2(data,tmin,tmax,index):
    tmintag = tmin.replace(':','_')
    tmaxtag = tmax.replace(':','_')

    folder = summary_folder(date)
    filesave = folder+f'Sismo_{index}_{tmintag}_{tmaxtag}.h5'
    rw.write_h5(filesave,data)
    print('Writing S2 file : ')
    print(filesave)


