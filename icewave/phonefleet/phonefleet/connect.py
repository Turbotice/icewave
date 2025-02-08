

import subprocess
import time
from pprint import pprint
import urllib

global local
local = '192.168'
global network
network = 2
global port
port = 8080

def ipbase():
    return f'{local}.{network}.'

def basefolder():
    return 'Bic25/'

def get_adress(phone):
    return ipbase()+str(100+phone)

def get_adresslist(phonelist):
    adresslist = []
    for phone in phonelist:
        adresslist.append(get_adress(phone))
    return adresslist
    
def scan():
    a = subprocess.run(['nmap','-sn',f'{local}.{network}.0/24'],text=True,capture_output=True)
    pprint(a.stdout.split('\n'))

    lines = a.stdout.split('\n')
    phonelist = []
    for i,line in enumerate(lines):
        if 'Host is up' in line:
            print('Host is up, '+lines[i-1])
            num = lines[i-1].split(f'{local}.{network}.')[-1]
            print(num)
            try:
                phonelist.append(int(num)-100)
            except:
                print('format does not correspond to a phone')
    return phonelist

def connect():
    localfolder = '/Users/stephane/Documents/PMMH/Telephones/test/'
    folder = '/storage/self/primary/'

    phonelist = scan()
    adrlist = get_adresslist(phonelist)
    print(phonelist)
    print("Number of phone connected :"+str(len(phonelist)))
    return adrlist,phonelist
    

def test_connect():
    c=subprocess.run(['adb','devices'],text=True,capture_output=True)
    log = c.stdout#.split('\n')
        #self.pt('toto')
    print(log)

if __name__=='__main__':
    launch_phyphox(0)
#    test_connect()
