

import phonefleet.run_gobannos as gob
import phonefleet.connect as connect

import socket, time,  urllib.request
import matplotlib.pyplot as plt
import numpy as np

import icewave.tools.rw_data as rw


import time

import threading


global results
results = {}

def run(phonelist,n=10000,iter=10):
    results = {}
    for i in range(iter):
        for phone in phonelist:
            print('')
            print(f"Phone {phone}, iteration {i}")
            #protocol command
            Dt = time_sync(phone,n=n)
#            thread = threading.Thread(target=tsync,args=(phone,results),kwargs={'n':n})
#            thread.start()
#            thread.join()
#            if thread.is_alive():
                # operation not completed
                # you must do cleanup here and stop the thread
#                print('Shit happen')
#            else:
#            Dt,duration = gob.time_sync(phone,n=50)
            if Dt is not None:
                ip = connect.get_adress(phone)
                gob.individual_stop(ip)
                
                result = {}
                result['phone']=phone
                result['iter']=i
                r = get_lag(Dt)
                result.update(r)
                results[f"{phone}_{i}"]=result
                time.sleep(0.1)
            else:
                print(f'Connexion failed with {phone}, iteration {i}')
        #print(results)
    #print(results)
    return results

def time_sync(phone,n=1000,timeout=0.1):
    ip = connect.get_adress(phone)
    port = 8080
    
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

    try:
        urllib.request.urlopen("http://" + ip + f":{port}/udp-sync").read()
    except:
        print(f'Enable to connect to {ip}')
        sock_receive.close()
        return None

    Dt = {}
    duration = []
    for i in range(1,4):
        Dt[i]=[]

    t0 = time.time()
    c=0
    lost=[]
    for i in range(n):
		#print(i)
        #if np.mod(i,100)==0:
        #    print(i)
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
    print(f'Duration : {np.round((tend-t0)*1000,decimals=3)} ms \n Number of packets lost : {c}/{n}')
    Dt['time']= t0
    sock_receive.close()
	
    duration = np.asarray(duration)
    for key in Dt.keys():
        Dt[key] = np.asarray(Dt[key])
    Dt['duration']=duration

    return Dt

def get_lag(Dt):
    duration = Dt['duration']
    t0 = Dt['time']
    tmedian = np.median(duration)

    tmax = tmedian*1
    print(f'Median duration of UDP request : {np.round(tmedian*1000,decimals=3)} ms')

    #plt.hist(duration,100)
    #plt.yscale('log')

    #plt.show()

    indices = np.where(duration<tmax)[0]
    #plt.plot(np.asarray(Dt[2])[indices])
    #plt.show()
	#print(Dt[2])
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


def save_table(Dts):
    pass

    # save table of Dt, with:
    # Each line correspond to a phone, 
    # Phone     Dtmean      Dtstd   Dtmax   Dtmin   Number of samples

def main():
#    phonelist = connect.scan()
    phonelist = [0,1,3,4,5]
    print(phonelist)

    for i in range(50):
        t1 = time.time()
        results = run(phonelist,n=50,iter=100)
        #function to write the Dts in a file, with following lines:
        t2 = time.time()
        print(f'elapsed time :{np.round(t2-t1,decimals=2)} s')
        print(results)

        savefolder = '/home/turbots/Documents/Bicwin2024/git/phonefleet/phonefleet/Bic25/Tsync/' 
        filename = savefolder+'tsync_'+str(int(np.round(time.time())))
        rw.writedict_csv(filename,results)

if __name__=='__main__':    
    main()