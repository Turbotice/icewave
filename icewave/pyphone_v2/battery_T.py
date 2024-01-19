import os
import subprocess
import shlex
import time
import pylab as plt
import pickle

import numpy as np
import connect    
    
def make_phone_dict(phonelist):
    phone_dict={}
    for phone in phonelist:
        phone_dict[phone]={}
        ip = connect.ipbase()+str(100+phone)
        port = str(5500+phone)
        phone_dict[phone]['ipfull']=ip+':'+port
	    
    return phone_dict
	
def smart_print(s,window=None):
    if window is not None:
        print('toto')
        print(s)
        window.pt(s)
    else:
        print(s)

def acquire(phone_dict,Tmax=10000,window=None,folder='Temperature/'):#default length is 3hours
    print(Tmax)

    basefolder = connect.basefolder()
    folder = basefolder+folder


    #print(window.console.text)
    data = {}


    for phone in phone_dict.keys():
	    data[phone]={}
	    data[phone]['time'] = []
	    data[phone]['T']=[]
	    data[phone]['level']=[]
	
    t0 = time.time()
    t1 = 0

    keys = ['time','T','level']
    
    date = str(t0).split('.')[0]

    filename = folder+'battery_data_'+date+'.txt'
    
    with open(filename, 'a') as f:
        line = "phone\t"
        for key in keys:
            line = line +key+"\t"
        line = line[:-1]+"\n"
        f.write(line)
                
    while t1<Tmax:
        for phone in phone_dict.keys():
            t1 = time.time()-t0
            data[phone]['time'].append(np.round(t1,decimals=1))

            ipfull = phone_dict[phone]['ipfull']
            output = subprocess.run(['adb','-s',ipfull,'shell','dumpsys','battery'],stdout=subprocess.PIPE)
            try:
	            batt = int(str(output.stdout).split('level')[1].split('scale')[0].split("\\n")[0][2:])
	            T = float(str(output.stdout).split('temperature')[1].split('technology')[0][2:5])/10
            except:
	            print('Battery access failed')	
            print(int(t1),phone,T,batt)
#            smart_print([int(t1),phone,T,batt],window=window)
            data[phone]['T'].append(T)
            data[phone]['level'].append(batt)
            
        time.sleep(1)

        
        #format of temperature data storage to be changed !
        with open(filename, 'a') as f:
            for phone in phone_dict.keys():
                line = str(phone)+" "
                for key in keys:
                    line = line + str(data[phone][key][-1])+"\t"
                line = line[:-1]+"\n"
                f.write(line)
            
    
if __name__=='__main__':
    adrlist,phonelist = connect.connect()

    phone_dict = make_phone_dict(phonelist)
    acquire(phone_dict,Tmax=100)
 
#for phone in phonelist.keys():
#	plt.plot(data['time'],data[phone]['T'],'ko')
	
#plt.show()
