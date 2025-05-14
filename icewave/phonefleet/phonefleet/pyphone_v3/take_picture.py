import os
import subprocess
import shlex
import time
import pylab as plt
import pickle
import connect_phone
output = subprocess.run(['adb','devices'],stdout=subprocess.PIPE)

ip = connect_phone.get_adress(7)#:'192.168.0.107:5507'

t0 = time.time()
t1 = 0
Tmax = 3600*12

data={}
data['time']=[]
	
while t1<Tmax:
	t1 = time.time()-t0
	print(t1)
	data['time'].append(t1)
	
	subprocess.run(['adb','-s',ip,'shell','input','tap','550','2000'],stdout=subprocess.PIPE)
	time.sleep(60)
	
   
#for phone in phonelist.keys():
#	plt.plot(data['time'],data[phone]['T'],'ko')
	
#plt.show()
