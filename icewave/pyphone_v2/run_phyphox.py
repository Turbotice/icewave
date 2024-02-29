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

import icewave.pyphone_v2.connect as connect

import aiohttp
import asyncio

import json

import urllib

from pprint import pprint


#import stephane.display.graphes as graphes
def get_base_url(ip):
    return "http://"+ip+":8080"

def get_phone(ip):
    return int(ip.split('.')[-1])-100

def run_serie(iplist,T=10,folder=''):
    basefolder = connect.basefolder()+'Phyphox/'
    folder = basefolder+folder

    #PP_CHANNELS = ["accX","accY","accZ"]
    #PP_CHANNELS_COUNT = len(PP_CHANNELS)
    print(folder)

    if not os.path.exists(folder):
    	os.makedirs(folder)
   
    loop = asyncio.get_event_loop()
    funlist = [clear, run, stop, save]
    #funlist = [run]#, stop, save]

    coroutine =  execute(clear,iplist)
    loop.run_until_complete(coroutine)

    coroutine =  execute(run,iplist)
    loop.run_until_complete(coroutine)

    t1=time.time()
    t2=t1
    print('Acquisition in progress ...')
    """
    while (t2-t1)<T:
        print('get data')
        coroutine =  execute(lambda x,y:save(x,y,folder=folder),iplist)
        loop.run_until_complete(coroutine)
        t2 = time.time()
    """
    time.sleep(T)
    print('Acquisition done')

    coroutine =  execute(stop,iplist)
    loop.run_until_complete(coroutine)

#    time.sleep(10)
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
    coroutine =  execute(fun,iplist,**kwargs)
    R = loop.run_until_complete(coroutine)
    return R

def run_save(fun,iplist,folder,**kwargs):
    loop = asyncio.get_event_loop()
    coroutine =  execute_save(fun,iplist,folder,**kwargs)
    R = loop.run_until_complete(coroutine)
    return R

async def execute(fun,iplist,**kwargs):
    async with aiohttp.ClientSession() as session:    
        R={}
        for ip_s in iplist:
            address = get_base_url(ip_s)
#            try:
            r = await fun(address,session)
            R[get_phone(ip_s)] = r#["result"]
            print(r)
        #fun(address,folder=folder)
        #   except:
                #raise
         #       R[get_phone(ip_s)] = {'result':False}
         #       print('Cannot connect to '+ip_s)
    #for r in R:
    #    print(r)
    return R

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
        s= "accX"#&&accY&&accZ&&locLat&&locLon"
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
    async with session.get(url) as resp:
        r = await resp.read()#read()#.text()
    return 

async def clear(address,session,folder=''):
    url = address + "/control?cmd=clear"
    return await request_json(url,session)
    
async def run(address,session,folder=''):
    url = address + "/control?cmd=start"
    return await request_json(url,session)

async def stop(address,session,folder=''):
    url = address + "/control?cmd=stop"
    return await request_json(url,session)

async def save(address,session,folder,name='test'):
    ip_s = address.split(':')[1][2:]
    print('Save data for '+ip_s)
    url = str(address) + '/export?format=1'

    print(url)
#    response = await request_save(url)
#    out = response.content
    out = await request(url,session)#requests.get(url=url)    
    print(out)
    if not os.path.isdir(folder):
        os.makedirs(folder)
    name = folder+'/'+name+'_'+ip_s.replace('.','_')
    with open(name+".zip", "wb") as file:
        file.write(out)
    print('Done')

async def request_save(url):
    response = requests.get(url)
    return await response

def save_single(address,folder,name='test'):
    ip_s = address.split(':')[1][2:]
    print('Save data for '+ip_s)
    url = str(address) + '/export?format=1'
    
    out = request_save(url) 

    basefolder = connect.basefolder()
    folder = basefolder+folder

    if not os.path.isdir(folder):
        os.makedirs(folder)

    name = folder+'/'+name+'_'+ip_s.replace('.','_')
    print(name)
    with open(name+".zip", "wb") as file:
        file.write(out)
    print('Done')


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

    
