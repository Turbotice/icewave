import numpy as np
import pylab as plt
import glob
import os
from pprint import pprint
import scipy.interpolate as interp
import scipy.signal as sig

import icewave.gps.gps as gps
import icewave.display.graphes as graphes
import phonefleet.data as dataphone

import icewave.tools.rw_data as rw
import icewave.tools.datafolders as df

from matplotlib.colors import LinearSegmentedColormap as LSC
import icewave.field.time as timest
import icewave.phone.sismo as sismo
import icewave.phone.waves as waves


import argparse
import h5py


def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate multi instruments data")
    parser.add_argument('-date', dest='date', type=str,default='0211',help='select date to process data')
    parser.add_argument('-step', dest='step', type=int,default=1,help='select Step to be performed')
    parser.add_argument('-num', dest='num', type=int,default=None,help='Select recording number for display phone map') 

    args = parser.parse_args()
    print(args)
    return args

def get_folder(date):
    #base = f'/media/turbots/BlueDisk/Shack25_local/Data/'
    base = df.find_path(year='2025')#
    print(base)
    return base +f'{date}/Phone/'

def get_savefolder(date):
    folder = get_folder(date)+'Results/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def get_phonelist(date):
    folder = get_folder(date)

    folders = glob.glob(folder+'*')

    phonelist=[]
    for f in folders:
        try:
            num = int(os.path.basename(f))
        except:
            continue
        phonelist.append(num)
    print(f'Found {len(phonelist)} phone folders : ')
    print(phonelist)
    return phonelist

def parse_filename(filename):
    name = filename.split('/')[-1]
    c0 = name.split('-')[-1]
    if len(name.split('-'))>1:
        c1 = name.split('-')[-2]
        num = int(name.split('-')[-3])
        typ = name.split('-')[-4]
    else:
        return None,None,False
    if '.' in typ:
        typ = typ.split('.')[-1]

    control = not(c0=='1.csv')
    return num,typ,control

def get_filelist(date,keytest='accelerometer',display=False):
    phonelist = get_phonelist(date)
    folder = get_folder(date)

    folder 
    files = {}
    for phone in phonelist:
        files[phone]={}
        f = folder+str(phone)
        filelist = glob.glob(f+'/*')

        for filename in filelist:
            num,typ,control = parse_filename(filename)
            if not control:
                continue
            if not typ in files[phone].keys():
                files[phone][typ]={}
            files[phone][typ][num]=filename
        if display:
            if keytest in files[phone].keys():
                print(phone,len(files[phone][keytest].keys()))
            else:
                print(phone,f'{keytest} not found')
    filelist = glob.glob(folder+'Tsync/*.txt')
    files['Tsync']=filelist
    return files

def load_lvl_0(files,phone,num,keys=None,header_only=False):
    r={}
    r['num']= num
    r['phone']=phone
    if keys==None:
        keys = files[phone].keys()
    for k in keys:
        if k in files[phone].keys():
            #print(files[phone][k])
            if num in files[phone][k]:
                filename = files[phone][k][num]
                r['folder']=os.path.dirname(filename)
                r['filename']=os.path.basename(filename)

                r['date'] = '-'.join(r['filename'].split('T')[0].split('-')[1:4])
                if not header_only:
                    dic = dataphone.load_data(filename)
                    for key in dic.keys():
                        r[key] = dic[key]
                else:
                    print('Only loading header')
            else:
                print(f'No data for phone {phone}, num {num}')
                return None
        else:
            print(f'No GPS for phone {phone}')
            return None
    return r

def load_lvl_1(date,phone,num,year='2025'):
    folder = get_folder(date)
    filename = glob.glob(folder+str(phone)+f'/*_phone{phone}_num{num}*.h5')[0]
    hf = rw.read_h5(filename)
#    f'{year}_'+date[:2]+'_'+date[2:]+f'_L1_phone{phone}_num{num}.h5'
    return hf

def get_h5_filename_N1(folder,date,phone,num):
    filename = folder+f'/'+date.replace('-','_')+'_L1_phone'+str(phone)+'_num'+str(num)+'.h5'
    return filename

def save_to_h5(r):
    print(r['folder'])
    print(r['date'])
    filename = get_h5_filename_N1(r['folder'],r['date'],r['phone'],r['num'])
    print(filename)
    if not os.path.exists(filename):        
        hf = h5py.File(filename, 'w')
        for key in r.keys():
            #print(key)
            hf.create_dataset(key, data=r[key])
        hf.close()
    else:
        print(f'Filename {filename} already exists')
        
def h5_exist(r):
    if r is None:
        return False
    filename = get_h5_filename_N1(r['folder'],r['date'],r['phone'],r['num'])
    return os.path.exists(filename)

def N0_to_N1(files,synctime,phone,num):
    r = load_lvl_0(files,phone,num,header_only=True)
    #r['date'] = timest.today_date(r['t0'])
    if not h5_exist(r):
        r = load_lvl_0(files,phone,num,header_only=False)
    else:
        print('File already exist, skipping')
        return None
    print(phone,num)
    #print(r.keys())
    if r is None:
        print(f'No data found for {phone}, {num}')
        return None
    if phone in synctime.keys():
        t0 = synctime[phone]
        if 'ta' in r.keys():
            r['t_sync']=r['ta']+t0
            r['date']=timest.today_date(r['t_sync'][0])
        print(r.keys())
        print(r['date'])
        r = check_status(r)                    
        save_to_h5(r)
                
def from_N0_to_N1(date,key='accelerometer',imin=0,overwrite=False):
    files = get_filelist(date)
    phonelist = list(files.keys())

    synctime = find_timetable(date)

    keys = ['accelerometer', 'gyroscope', 'magnetic_field', 'gps']
    print(synctime.keys())
    for phone in phonelist[imin:]:
        if key in files[phone].keys():
            for num in files[phone][key].keys():
                N0_to_N1(files,synctime,phone,num)
        else:
            print(f'No acceleration data for phone {phone}')    

def from_N1_to_N2(date):
    folder = summary_folder(date)
    filename = folder+'Phone_Table.txt'
    params = parse_phone_table(filename)

    for i,param in enumerate(params):
        tag = param['tag']
        if tag=='waves':
            data = waves.load_data(param)#date,phonelist,nums,tmin,tmax,orientation)
            data = waves.smooth(data)
            waves.save_W2(data,param['tmin'],param['tmax'],i)
        elif tag=='sismo_active':
            data = sismo.load_data(param)
            sismo.save_S2(data,param['tmin'],param['tmax'],i)

        elif tag=='sismo_passive':
            data = sismo.load_data(param)
            sismo.save_S2(data,param['tmin'],param['tmax'],i)

        else:
            print('Tag key not recognized')

def summary_folder(date):
    #base = f'/media/turbots/BlueDisk/Shack25_local/Data/'
    base = df.find_path(year='2025')#
    print(base)
    return base +f'{date}/Phone/Summary/'


def parse_list(string):
    stringlist = string.split(',')
    numlist = []
    for elem in stringlist:
        if ':' in elem:
            imin,imax = elem.split(':')
            numlist = numlist+list(range(int(imin),int(imax)+1))
        else:
            numlist.append(int(elem))
    #print(numlist)
    return numlist

def parse_phone_table(filename):
    data = rw.read_csv(filename,delimiter='\t')
    print(data)
    dic = rw.csv2dict(data)

    params = []
    for i,tag in enumerate(dic['tag']):
        param={}
        param['tag']=tag
        for key in ['phonelist','nums']:
            print(key,dic[key][i])
            param[key] = parse_list(dic[key][i])
        param['orientation']=dic['orientation'][i].split(',')
        param['tmin'] = dic['tstart'][i]
        param['tmax'] = dic['tend'][i]
        params.append(param)
    return params

def select_data(param):
    #select the phonelist,
    date = param['date']
    #get_filelist(date,keytest='accelerometer',display=False)

def check_status(r):
    status = {}
    status['oural_time']=('t_sync' in r.keys())
    status['GPS_time']=('t_GPS' in r.keys())
    status['position_GPS']=('gpslon' in r.keys())  
    status['position_Garmin']=('gargps_lon' in r.keys())  
    status['complete'] = np.sum([key in r.keys() for key in ['ta','tm','tg','tgps']])==4
    for key in status.keys():
        r[key]=status[key]
    return r

def on_field(files,phone,num,synctime,tstart,tend):
    tmin = int(files[phone]['gps'][num].split('.')[-2].split('-')[-2])/10**6+synctime[phone]
    tmax = int(files[phone]['gps'][num].split('.')[-2].split('-')[-1])/10**6+synctime[phone]
    print(timest.display_time(timest.today_time([tmin,tmax])))
    print(tmin,tmax)

def attribute_tag(r):
    pass


def find_timetable(date):
    phone_to_sync = get_phonelist(date)
    folder = get_folder(date)
    filelist = glob.glob(folder+'Tsync/*.txt')

    timetable={}
    synctime={}
    for filename in filelist:
        phonelist = []
        synctable = rw.read_csv(filename,delimiter=',')
        synctable = rw.csv2dict(synctable)

        for key in synctable.keys():
            phone,i = key.split('_')
            phone = int(phone)
            if (not phone in phonelist) and (phone in phone_to_sync):
                phonelist.append(phone)
                phone_to_sync.remove(phone)
                timetable[phone]=[]
            if phone in timetable.keys():
                timetable[phone].append(synctable[key]['tlag'])

        for phone in phonelist:
            synctime[phone] = np.mean(timetable[phone])

    print(f'Phones not sync : {phone_to_sync}')
    return synctime

def display_times(files,synctime,phone):
    if phone in synctime.keys():
        t0 = synctime[phone]
    else:
        print(f'No time sync for phone {phone}')
    if 'gps' in files[phone].keys():
        for num in files[phone]['gps'].keys():
            r = load_lvl_0(files,phone,num,keys=['gps'])
            tmin = r['tgps'][0]+t0
            tmax = r['tgps'][-1]+t0

            #times = timest.to_UTC([tmin,tmax])
            times = timest.today_time(times)
            times = timest.display_time(times)
            print(num,times)
    else:
        print(f'No gps signal for phone {phone}!')

def is_moving(r,maxmotion=10,meanmotion=3):
    t = r['tgps']
    X,Y = gps.project(r['gpslon'],r['gpslat'],meter=True)
    d = np.sqrt(X**2+Y**2)    
    #plt.plot(t,d,'k.')
    
    dmax = np.round(np.max(d),decimals=1)
    dmean = np.round(np.mean(d),decimals=1)

    if dmax<maxmotion and dmean<meanmotion:
        move=False
    else:
        move=True
    return move

def position(r):
    X,Y = gps.project(r['gpslon'],r['gpslat'])
    Xmean = np.mean(X)
    Ymean = np.mean(Y)
    return Xmean,Ymean

def display_position(r,ax=None,eps = 10**(-6)):
    X,Y = position(r)
    phone = r['phone']
    ax.plot(X,Y,'rs',markersize=10)
    plt.plot(X,Y,'rs',markersize=10)
    
    sign = int((np.mod(phone,2)-0.5)*2)
    ax.text(X-eps/4,Y-eps*sign,str(phone),color='r',fontsize=20)

def find_relevant_map(r):
    sites = ['haha','capelans']
    Long = r['gpslon'][0]
    Lat = r['gpslat'][0]
    for name in sites:
        b = gps.check_box(Long,Lat,gps.boxes(name))
        print(b)
        if b:
            return name
    return None

def situation_map(files,num,date):
    phonelist = files.keys()
    fig,ax = plt.subplots(figsize=(20,10))

    phone = list(phonelist)[0]
    r = load_lvl_0(files,phone,num,keys=['gps'])
    name = find_relevant_map(r)
    print(name)
    ax,figs = gps.display_standard_map(ax,name,title=date)

    for phone in phonelist:
        r = load_lvl_0(files,phone,num,keys=['gps'])
        if r is not None:
            move = is_moving(r)
            #print(phone,move)
            if move==False:
                display_position(r,ax=ax)
    return fig,ax,figs

def moving_table(files,imin=1,imax=None):
    phonelist = files.keys()
    nphone = len(phonelist)
    imax=0
    
    length=[]
    key = 'accelerometer'
    for phone in phonelist:
        if key in files[phone].keys():
            n = len(files[phone][key].keys())
        else:
            print(f'no accelerometer data for phone {phone}')
            n=0
        length.append(n)
    imax = int(np.nanmedian(length))-3
    print(f'Median number of files : {imax}')
    
    indices = range(imin,imax)
    nt = len(indices)
    moving = np.zeros((nphone,nt))
    for i,phone in enumerate(phonelist):
        for j,num in enumerate(indices):
            r = load_lvl_0(files,phone,num,keys=['gps'])
            if r is not None:
                move = is_moving(r)
                #print(phone,move)       
                moving[i,j]=int(move)
            else:
                moving[i,j]=np.nan
    return moving,indices

def moving_map(moving,indices,phonelist):
    Indices,Phones = np.meshgrid(indices,phonelist)
    fig,ax = plt.subplots(figsize=(15,10))

    colors = [(0, 1, 0), (1, 0, 0)]# R -> G -> B
    cmap_name = 'onoff'
    n_bin = 2
    cmap = LSC.from_list(cmap_name, colors, N=n_bin)

    sc = ax.pcolormesh(Indices,Phones,moving,cmap=cmap)
    #plt.colorbar()
    cbar = fig.colorbar(sc, ax=ax,ticks=[0,1])
    cbar.ax.set_yticklabels(['Immobile','Moving'])

    figs = graphes.legende('number of recording','Phone number','')
    ax.set_xticks(indices)
    ax.set_yticks(phonelist)

    return fig,ax,figs

def situation(date,imax=None,num=None,save=True):
    files = get_filelist(date)
    phonelist = list(files.keys())
    savefolder = get_savefolder(date)

    moving,indices = moving_table(files,imax=imax)
    fig,ax,figs = moving_map(moving,indices,phonelist)
    if save:
        graphes.save_figs(figs,savedir=savefolder,prefix='Moving_map_',suffix='_'+date,overwrite=True)

    if num is None:
        s = input('Enter number to select :')
        num = int(s)

    fig,ax,figs = situation_map(files,num,date)
    if save:
        graphes.save_figs(figs,savedir=savefolder,prefix='Situation_map',suffix='_'+date,overwrite=True)

def generate_N1_selective():
    base = df.find_path(year='2025')#
    filename = base + 'Summary/Phone/Phone_Indexes_to_keep.pkl'
    tokeep = rw.load_pkl(filename)

    for date in tokeep.keys():
        files = get_filelist(date,keytest='accelerometer',display=False)
        synctime = find_timetable(date)
        print(date)
        for k in tokeep[date]:
            for phone in k['files'].keys():
                for num in k['files'][phone]:
                    N0_to_N1(files,synctime,phone,num)

if __name__=='__main__':
    args = gen_parser()
    generate_N1_selective()
    #situation(args.date,num=args.num,save=True)
    #from_N0_to_N1(args.date)
