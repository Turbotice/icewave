

import subprocess
import time
import numpy as np
import pylab as plt
import os
import pickle
import timeit
import glob	

import connect

def run():
    adrlist,phonelist = connect.connect()
    return compute_timetable(adrlist,phonelist)

def gen_autosyn(adrlist,phonelist):
    base = '/home/turbots/Documents/Codes/Python/pyphone_v2/'
    filename = base + 'autosyn.sh'
    with open(filename,'w') as f:
        incipit = "#!/usr/bin/env bash\n"
        f.write(incipit)
        print(filename)
        fold = 'timelog_2/'
        for adr,phone in zip(adrlist,phonelist):
    #        line = 'adb -s '+str(adr)+' shell "echo $EPOCHREALTIME &&date +%s.%N && echo $EPOCHREALTIME">'+fold+'t'+str(phone)+'.txt\n'
            line = 'time (echo $EPOCHREALTIME && adb -s '+str(adr)+' shell /system/bin/time -p echo \$EPOCHREALTIME && echo $EPOCHREALTIME)>>'+fold+'t'+str(phone)+'.txt 2>&1\n'
            print('get time phone'+str(phone))
            f.write(line)
    return filename

def read_timetable(folder):
    filelist = glob.glob(folder+'*.txt')
    savefile = folder+'timetable.txt'
    
    tablist = ""
    with open(savefile,'w') as fi:
        for filename in filelist:
            if 'timetable' in filename:
                continue
            with open(filename,'r') as f:
                num = filename.split('/')[-1].split('.txt')[0][1:]
                try:
        #            t1 = float(f.readline().split('\n')[0].replace(',','.'))
        #            t2 = float(f.readline().split('\n')[0])
        #            t3 = float(f.readline().split('\n')[0].replace(',','.'))
                    t1 = float(f.readline().split('\n')[0].replace(',','.'))
                    t2 = float(f.readline().split('\n')[0])
                    dt3 = float(f.readline().split('\n')[0].split('real ')[1])
                    t3 = t2+dt3
                    f.readline()
                    f.readline()
                    t4 = float(f.readline().split('\n')[0].replace(',','.'))
                    f.readline()
                    dt4check =  float(f.readline().split('\n')[0].replace(',','.').split('m')[1][:-1])

                    Dt = (t2+t3-t1-t4)/2
        #            print(t4-t1-dt4check)
        #            print(Dt)
            #        Dt = t2-t1
                except:
                    #raise
                    print('Cannot read phone #'+str(num)) 
                    Dt=0.0
            line = str(num)+'\t'+str(Dt)+'\n'
            tablist = tablist+line#append(line.split('\t'))
            fi.write(line)
        
    fi.close()
    return tablist

def compute_timetable(adrlist,phonelist):
    #write bash script 
    filename = gen_autosyn(adrlist,phonelist)

    print('run autosyn.sh')
    subprocess.call(filename)
    print('done')
    
    folder = os.path.dirname(filename)
    tablist = read_timetable(folder)
    return tablist
    
if __name__=='__main__':
    t0 = time.time()
    t1 = t0
    T = 1800
    c=0

    adrlist,phonelist = connect.connect()
    filename = gen_autosyn(adrlist,phonelist)

    while (t1-t0)<T:
        t1 = time.time()
        print(str(c)+', run autosyn.sh')
        subprocess.call(filename)
        c=c+1

