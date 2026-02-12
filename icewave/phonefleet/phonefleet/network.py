import subprocess
from pprint import pprint
import urllib

global local
local = '192.168'
global network
network = 2
global port
port = 8080

def scan():
    a = subprocess.run(['nmap','-sn',f'{local}.{network}.0/24'],text=True,capture_output=True)
    pprint(a.stdout.split('\n'))

    lines = a.stdout.split('\n')
    phonelist = []
    for i,line in enumerate(lines):
        if 'Host is up' in line:
            print('Host is up, '+lines[i-1])
            num = lines[i-1].split(f'{local}.{network}')[-1]
            try:
                int(num)
                phonelist.append(num[-2:])
            except:
                print('format does not correspond to a phone')
    return phonelist
    
def get_ip(phone,base='1'):
    return f'{local}.{network}.{base}{phone}'
    
def test_network(phonelist):
    c = 2
    found = []
    for phone in phonelist:
        print(phone)
        ip = get_ip(phone)
        a = subprocess.run(['ping','-c',str(c),ip],text=True,capture_output=True)
        #print(a.stdout)
        lines = a.stdout.split('\n')
        line = lines[-3]
        pprint(line)
        
        value = int(line.split(' received,')[0][-1])
        print(value)
        if value==c:
            found.append(phone)
            
    return found
    
def get_status(phone):
    ip = get_ip(phone,base='1')
    a = urllib.request.urlopen(f"http://{ip}:{port}/status").read()
    print(a)
    return a

def test_gobannos_link(phonelist):
    for phone in phonelist:
        a=get_status(phone)
        stop = (a==b'STOPPED')
        print(phone,a,stop)
        
        
    
