import subprocess
import os


date = '20241127'

global network
network = 214

def ls(date,phone,network=network):
    ip = f"192.168.{network}.1{phone}"
    path = 'storage/self/primary/DCIM/Camera/'
    s = subprocess.run(['adb','-s',ip,'shell','ls',path+'VID_'+date+'*'],text=True,capture_output=True)

    return s.stdout.split('\n')[:-1]
    

def download(filename,savefolder,phone,network=network):
    ip = f"192.168.{network}.1{phone}"

    name = filename.split('/')[-1]
    print(name)
    s = subprocess.run(['adb','-s',ip,'pull',filename,savefolder+name],text=True,capture_output=True)
#    print(s)


phone = '00'
filelist = ls(date,phone)

base = f'/Volumes/Fabien_2024/Grenoble/{date}/Phone_{phone}/'


for i,filename in enumerate(filelist):
    savefolder = base + f'{i+1}/'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    print(filename)
    download(filename,savefolder,phone)
    
print(filelist)
