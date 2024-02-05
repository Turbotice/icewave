

import icewave.tools.datafolders as df
import icewave.tools.browse as browse
import icewave.pyphone_v2.connect as connect

import subprocess

def download(date='',folder=''):
    phonedict= get_datalist(date=date,folder=folder)
    
    base = df.instrument_folder(key='T')
    for phone in phonedict.keys():
        #the numbers could be mixed up if some phones have not recorded data !!
        id = phonedict[phone]['id']
        base_adb = phonedict[phone]['folder']
        folderlist = phonedict[phone]['foldersave']

        for i,folder in enumerate(folderlist):
            dest_folder = base+df.ndigit(i,n=4)+'/T'+str(df.ndigit(phone,n=2))+''#+folder.split('/')[-1]
            browse.create_folder(dest_folder)
            print(base_adb+folder,dest_folder)
            exemple = ['adb','pull',base_adb+folder+'',dest_folder]
            output = subprocess.run(exemple)#,stdout=subprocess.PIPE)


def get_datalist(date='',folder=''):
    if folder=='':
        phonedict = ls_all()
    else:
        phonedict = ls_all(folder)
    if date=='':
        date = df.get_current_date()
    date = df.date_phox(date)
    

    for phone in phonedict.keys():
        phonedict[phone]['foldersave']=[]
        for folder in phonedict[phone]['folderlist']:
            if date in folder:
                #folder = folder.replace(' ','\ ')
                phonedict[phone]['foldersave'].append(folder)
    return phonedict

def ls_all(folder='storage/self/primary/phyphox/'):
    phonedict = connect.get_connected()    

    for phone in phonedict.keys():
        id = phonedict[phone]['id']
        out = ls(id,folder=folder)
        phonedict[phone]['folder']=folder
        phonedict[phone]['folderlist']=out
    return phonedict


def ls(id,folder=''):
    print(id)
    exemple = ['adb','-s',id,'shell','ls',folder]
    output = subprocess.run(exemple,stdout=subprocess.PIPE)
    return str(output.stdout)[2:].split(('\\n'))
        #output#out = output.split('/n')

if __name__=='__main__':
    ls()
