


import os
import time
import numpy as np
from datetime import datetime
import tsync_modif as tsync
import run_gobannos as gob
import connect
import rw_data

import socket
import platform

def get_folder():
    ostype = platform.platform().split('-')[0]
    osname = socket.gethostname()
    if osname=='localhost':
        folder = '/storage/emulated/0/Documents/'
    elif ostype=='macOS':
        folder = 'Phonefleet/'
    elif ostype=='linux':
        folder= 'Phonefleet/'
    else:
        folder= 'Phonefleet/'
    return folder

def set_network():
    s = input('Type sub network number')
    network = int(s)
    return network

def set_phonelist():
    s = input('Type phonelist to adress :')
    lis = s.split(':')
    phonelist = []
    for i,l in enumerate(lis[:-1]):
        phonelist = phonelist+list(range(int(lis[i]),int(lis[i+1])+1))
    phonelist = list(set(phonelist))
    return phonelist

def get_filelist(network,phonelist,s,display=True):
    phone = int(s.split(' ')[1])
    if phone in phonelist:
        ip = connect.get_adress(phone,network=network)
        filelist = gob.get_file_list(ip)

        if display:
            for i,filename in enumerate(filelist):
                print(i,filename)
        return filelist
    else:
        print(f'Phone {phone} not found !')
        return []

def get_file(network,phonelist,s):
    try:
        phone = int(s.split(' ')[1])
        num = int(s.split(' ')[2])
    except:
        print('argument not valid, specify phone and file number,')
        print('exemple : pull 5 16')
        return None
    s = f'ls {phone}'
    filelist = get_filelist(network,phonelist,s)
    if len(filelist)==0:
        return None
    if num<len(filelist):
        filename = filelist[num]

        ip = connect.get_adress(phone,network=network)
        a = gob.get_file(ip,filename)
        data = a.decode('utf-8')
        
        save_file(filename,data,phone)
    else:
        print(f'Num {num} greater than the number of files')
        return None

def save_file(filename,data,phone):
    name = os.path.basename(filename)
    savefile = get_folder()+f'P{phone}_D{name}'
    print(savefile)
    print(data.split('\n')[0])
#    print(data)
    with open(savefile,'w') as f:
        f.write(data)
    
def check_status(network,phonelist):
    for phone in phonelist:
        get_status(phone,network=network)

def get_status(phone,network=1):
    ip = connect.get_adress(phone,network=network)
    status = gob.get_status(ip)
    print(status)
    
def stop_phonelist(network,phonelist):
    for phone in phonelist:
        ip = connect.get_adress(phone,network=network)
        status = gob.get_status(ip)
        print(status)
        if not status==b'STOPPED':
            gob.individual_stop(ip)

def start_phonelist(network,phonelist):
    for phone in phonelist:
        ip = connect.get_adress(phone,network=network)
        status = gob.get_status(ip)
        print(status)
        gob.individual_start(ip)

def time_sync(network,phonelist,iter=5):
    savefolder = get_folder()+'Gobannos/'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

    for i in range(iter):
        results={}
        for phone in phonelist:
            ip = connect.get_adress(phone,network=network)
            status = gob.get_status(ip)
            print(status)
            print(ip)
            Dt = tsync.time_sync_ip(ip,n=100,timeout=0.1)
            print(Dt)
            if Dt is not None:
                result = get_lag(Dt)
                results[phone]=result
        maintenant = datetime.now()
        filename = savefolder+f'tsync_'+str(maintenant.strftime("%Y-%m-%d %H:%M:%S"))+'.txt'
        print(filename)
        rw_data.writedict_csv(filename,results)

def get_lag(Dt):
    duration = Dt['duration']
    t0 = Dt['time']
    tmedian = np.median(duration)

    tmax = tmedian*1
    print(f'Median duration of UDP request : {np.round(tmedian*1000,decimals=3)} ms')

    indices = np.where(duration<=tmax)[0]
    tlag1 = np.asarray(Dt[1])[indices]
    tlag2 = np.asarray(Dt[2])[indices]
    tlag3 = np.asarray(Dt[3])[indices]
    tlag = (tlag1+tlag3)/2
    Dt = np.median(tlag)

    results={}
    results['tlag'] = Dt
    results['dtmedian'] = tmedian
    results['tmin'] = np.min(tlag)
    results['tmax'] = np.max(tlag)
    results['tstd'] = np.std(tlag)
    results['n'] = len(duration)
    results['t0'] = t0
    return results
                
def choose(network,phonelist):
    print('chose action among : ')
    actions = ['network','phones','status','time','start','stop','exit','ls options','pull options']
    descriptions = ['','','','','','','','','']
    for action,description in zip(actions,descriptions):
        print(action,description)
    s = input('')

    if s=='network':
        network = set_network()
        return s,network
    if s=='phones':
        phonelist = set_phonelist()
        return s,phonelist
    elif s=='status':
        print(f'Available phonelist : {phonelist}')
        check_status(network,phonelist)
    elif s=='time':
        time_sync(network,phonelist)
    elif s=='start':
        start_phonelist(network,phonelist)
    elif s=='stop':
        stop_phonelist(network,phonelist)
    elif s[:2]=='ls':
        get_filelist(network,phonelist,s)
    elif s[:4]=='pull':
        data = get_file(network,phonelist,s)
    elif s=='exit':
        print("exit")
    else:
        print('command not known, do nothing')
    return s,None

def defaults():
    ostype = platform.platform().split('-')[0]
    osname = socket.gethostname()
    if osname=='localhost':
        phonelist = [1,2,3,4,5]
        network = 189
    elif ostype=='macOS':
        phonelist = [8]
        network = 1
    else:
        phonelist = [5,6,7,8]
        network = 142
    return phonelist,network

def main():
    phonelist,network = defaults()

    s=''
    while not s=='exit':
        s,output = choose(network,phonelist)
        if s=='phones':
            phonelist=output
        elif s=='network':
            network = output
        print(phonelist)
        print(network)
    #action to code :
    #define phonelist (phones)
    #check connection (status)
    #timesync (time)
    #start acquisition (start)

if __name__=='__main__':
    main()
