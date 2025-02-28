# 2020-05-04 
# https://phyphox.org/wiki/index.php/Remote-interface_communication
# http://10.10.10.21:8080/control?cmd=start
# http://10.10.10.21:8080/control?cmd=stop
# http://10.10.10.21:8080/control?cmd=clear
# http://10.10.10.21:8080/get?pressure
# http://10.10.10.21:8080/get?pressure=201.69097|p_time&p_time=201.69097
# http://10.10.10.21:8080/get?pressure=262.40482|p_time&p_time=262.40482
 # http://10.10.10.21:8080/get?pressure=full&p_time=full

# Alle Messwerte ab secunde 35.88 (aus letzten Daten)
#  http://10.10.10.21:8080/get?pressure=full&p_time=35.88
# http://10.10.10.21:8080/get?illum
# http://10.10.10.21:8080/get?illum=57.182803|illum_time&illum_time=57.182803
import requests
import time
import zipfile
import numpy as np
import glob
import os


import phonefleet.connect as connect

#import aiohttp
import asyncio
import re

import json
import urllib

from pprint import pprint

global port
port = 8080

# accessible functions
commands = {'start':'/start?name',
            'stop':'/stop',
            'status':'/status',
            'sync':'/kick-sync?threshold'}
          
#import stephane.display.graphes as graphes
def get_base_url(ip):
    return "http://"+ip+f":{port}"

def get_phone(ip):
    return int(ip.split('.')[-1])-100

def get_status(ip):
#    ip = get_adress(phone)
    address = f"http://{ip}:{port}/status"
    print(address)
    a = urllib.request.urlopen(address).read()
    print('Status :'+str(a))
    return a

def get_file_list(ip):
    a = urllib.request.urlopen(f"http://{ip}:{port}/list-files").read()
    s = a.decode('utf-8')
    filelist = s[1:-1].split(', ')
    return filelist

def get_today_file_list(ip):
    data = '2025-01-21'
    filelist = gob.get_file_list(ip)
    
    rx = re.compile(r'experiment-{date}*')
    filelist = list(filter(rx.search, filelist))
    pprint('Number of today files :'+str(len(filelist)))
    return filelist

def load_data():
    pprint(filelist)

def individual_stop(ip):
    a = urllib.request.urlopen(f"http://{ip}:{port}"+commands['stop']).read()

def time_sync(phone,n=200,timeout=0.02):
    ip = connect.get_adress(phone)
    
    do_nothing = 0
    do_nothing_command = do_nothing.to_bytes(4, 'big', signed=False)
    respond = 1
    respond_command = respond.to_bytes(4, 'big', signed=False)
    stop = 2
    stop_command = stop.to_bytes(4, 'big', signed=False)

    socket.socket(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock_send = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock_receive = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock_receive.bind(("", 5001))
    sock_receive.settimeout(timeout)

    sock_send.sendto(do_nothing_command, (ip, 5000))

	# first query to warm up
	#except:
	#    print('initialisation fail')
	# activate udp sync on phone

    urllib.request.urlopen("http://" + ip + f":{port}/udp-sync").read()

    Dt = {}
    duration = []
    for i in range(1,4):
        Dt[i]=[]

    t0 = time.time()
    c=0
    for i in range(n):
		#print(i)
        if np.mod(i,100)==0:
            print(i)
        t1 = time.time_ns()
        sock_send.sendto(respond_command, (ip, 5000))
        t2 = time.time_ns()
        try:
            t_phone_bytes = sock_receive.recv(8)
        except:
            c=c+1
            #print(f'time longer than {timeout*1000} ms')
            continue
        t3 = time.time_ns()
		#print(t3)
        t_phone = int.from_bytes(t_phone_bytes, byteorder='big')
		#print(str((t_phone-t1)/1000000) + "    " + str((t3-t_phone)/1000000) + "    " + str((t3-t1)/1000000))
        duration.append((t3-t1)*10**(-9))
        Dt[1].append((t1-t_phone)*10**(-9))
        Dt[2].append((t2-t_phone)*10**(-9))
        Dt[3].append((t3-t_phone)*10**(-9))
		# stop sync
    sock_send.sendto(stop_command, (ip, 5000))
    tend = time.time()
    print('')
    print(f'Duration : {np.round((tend-t0)*1000,decimals=3)} ms \n Number of packets lost : {c}/{n}')

    sock_receive.close()
	
    duration = np.asarray(duration)
    for key in Dt.keys():
        Dt[key] = np.asarray(Dt[key])
    Dt['duration']=duration
    return Dt

def get_lag(Dt):
    duration = Dt['duration']
    tmedian = np.median(duration)
    tmax = tmedian*1.0

    indices = np.where(duration<tmax)[0]
	#print(Dt[2])
    tlag = np.mean(np.asarray(Dt[2])[indices])
    return tlag


def run_serie(iplist,T=10,folder=''):
    basefolder = connect.basefolder()+'Data/'
    folder = basefolder+folder

    if not os.path.exists(folder):
    	os.makedirs(folder)
   
    loop = asyncio.get_event_loop()
    funlist = [clear, run, stop, save]
    #funlist = [run]#, stop, save]

    coroutine =  execute(run,iplist)
    loop.run_until_complete(coroutine)

    t1=time.time()
    t2=t1
    print('Acquisition in progress ...')
    time.sleep(T)

    coroutine =  execute(stop,iplist)
    loop.run_until_complete(coroutine)

    print('Acquisition done')

def run_save(iplist,T=10,folder=''):
    basefolder = connect.basefolder()+'Phyphox/'
    folder = basefolder+folder

    #PP_CHANNELS = ["accX","accY","accZ"]
    #PP_CHANNELS_COUNT = len(PP_CHANNELS)
    print(folder)

    if not os.path.exists(folder):
    	os.makedirs(folder)
   
    loop = asyncio.get_event_loop()
    print('Save data ...')

    coroutine = execute(lambda x,y:save(x,y,folder=folder),iplist)
    loop.run_until_complete(coroutine)
#        print(fun.__name__,'Over')

def run_config(iplist):
    loop = asyncio.get_event_loop()
    coroutine = execute(lambda x,y:get_config(x,y),iplist)
    R = loop.run_until_complete(coroutine)
    return R

def run_fun(fun,iplist,**kwargs):
    loop = asyncio.get_event_loop()
    coroutine =  execute(fun,iplist)
    R = loop.run_until_complete(coroutine)
    return R

async def run_try(address,fun,session):
    #r = await fun(address,session)
    try:
        r = await fun(address,session)
        print(r)
    #fun(address,folder=folder)
    except:
        print('Cannot connect to '+address)
        r =  {'result':False}
    return r

async def execute(fun,iplist,**kwargs):
    async with aiohttp.ClientSession() as session:    
        R={}
        urls = [get_base_url(ip) for ip in iplist]
        R = await asyncio.gather(*map(lambda url:run_try(url,fun,session), urls))
#            R[get_phone(ip_s)] = r#["result"]
        print('Results :'+str(R))

        out = {}
        for i,ip in enumerate(iplist):
            phone = get_phone(ip)
            print(R[i])
            for r in R[i]:
                out[phone]=R[i]
    return out

async def get_config(address,session):
#    adress=get_base_url(ip)
    url = address+"/config"
    return await request_json(url,session)

async def request_json(url,session):
    print(url)
    async with session.get(url,timeout=5) as resp:
        r = await resp.json()#.text()
    return r  

def run_get(iplist,name):
    loop = asyncio.get_event_loop()
    #coroutine = execute(get,iplist,name=name)
    coroutine = execute(lambda x,y:get(x,y,name=name),iplist)

    R = loop.run_until_complete(coroutine)
    return R

async def get(address,session,name=''):
    #name=kwargs.get('name')
#    adress=get_base_url(ip)
    print(name)
#    name='all'
    if name=='location':
        s = "locLat&&locLon"
        url = address+"/get?"+s
        return await request_json(url,session)
    if name=='all':
        #rlist = await get_config(address,session)
        s= "accX&&accY&&accZ&&locLat&&locLon"
        #for r in rlist['buffers'][1:4]:
        #s=s+r['name']+"&&"
        #s=s[:-2]
        print(s)
        url = address+"/get?"+s
        return await request_json(url,session)
    else:
        print('not defined')
        s=""
        return None

async def request(url,session):
    print(url)
    async with session.get(url,timeout=600) as resp:
        r = await resp.read()#.text()
    return r

async def clear(address,session,folder=''):
    url = address + commands['clear']#"/control?cmd=clear"
    return await request_json(url,session)
    
async def run(address,session,folder=''):
    url = address + commands['start']#
    return await request(url,session)

async def stop(address,session,folder=''):
    url = address + commands['stop']#
    return await request(url,session)

async def status(address,session,folder=''):
    url = address + commands['status']#
    return await request(url,session)

async def save(address,session,folder,name='test'):
    ip_s = address.split(':')[1][2:]
    print('Save data for '+ip_s)
    url = address + '/export?format=1'

    out = await request(url,session)#requests.get(url=url)    
    print(ip_s)
    if not os.path.isdir(folder):
        os.makedirs(folder)
    name = folder+'/'+name+'_'+ip_s.replace('.','_')
    with open(name+".zip", "wb") as file:
        file.write(out)
    print('Done')

    return {'result':True}


async def get_buffer(address,session,buffername):
    url = adress+"/get?"+buffername
    return await request(url,session)

def get_buffer():
    pass
    
def load(folder):
    filelist = glob.glob(folder+'*.zip')
    print(filelist)

    data = {}
    for filename in filelist:
        with zipfile.ZipFile(filename,"r") as zip_ref:
            foldersave = filename.split('.')[0]
            zip_ref.extractall(foldersave)
        datafile = foldersave+'/Raw Data.csv'
        d = np.loadtxt(datafile, delimiter=',',skiprows=1)
        data[filename] = d
    return data

def showdata(data):
    import pylab as plt
    n = len(data.keys())

    colors = ['g','b','y']
    coords = ['x','y','z']

    fig,axs = plt.subplots(nrows=3,ncols=n,figsize=(10,8))
    for j,key in enumerate(data.keys()):
        d = data[key]
        for i,(ax,color,coord) in enumerate(zip(axs[:,j],colors,coords)):
            ax.plot(d[:,0],d[:,i+1],color)
        #    if j==0:
        #        graphes.legende('','$a_'+coord+'$ (m/s$^2$)','',ax=ax)
        #graphes.legende('Times (s)','$a_'+coord+'$ (m/s$^2$)','',ax=ax)

    plt.show()


if __name__=='__main__':
    basefolder = connect.basefolder()
    date = '121223_'+str(time.time())
    folder = basefolder+'Phyphox/'+'date'
    # [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,22,23,24,25]#
    tel_numbers =np.arange(1,64,dtype=int)
    print(tel_numbers)
    base = connect_phone.ipbase()
    print(base)
    iplist = [base + str(100+tel) for tel in tel_numbers]
    run_serie(iplist,T=600,folder=date)
    #data = load(date+'/')
    #showdata(data)

    
