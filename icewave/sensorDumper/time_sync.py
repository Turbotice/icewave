import socket, time,  urllib.request
import numpy as np

def get_adress(network,phone):
	return f"192.168.{network}.{phone}" # phone address

def time_sync(network,phone,n=200,timeout=0.02):
	address = f"192.168.{network}.{phone}" # phone address

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

	sock_send.sendto(do_nothing_command, (address, 5000))

	# first query to warm up
	#except:
	#    print('initialisation fail')
	# activate udp sync on phone

	urllib.request.urlopen("http://" + address + ":8080/udp-sync").read()

	Dt = {}
	duration = []
	for i in range(1,4):
		Dt[i]=[]

	t0 = time.time()
	for i in range(n):
		if np.mod(i,100)==0:
			print(i)
		t1 = time.time_ns()
		sock_send.sendto(respond_command, (address, 5000))
		t2 = time.time_ns()
		t_phone_bytes = sock_receive.recv(8)
		t3 = time.time_ns()
		#print(t3)
		t_phone = int.from_bytes(t_phone_bytes, byteorder='big')
		#print(str((t_phone-t1)/1000000) + "    " + str((t3-t_phone)/1000000) + "    " + str((t3-t1)/1000000))
		duration.append((t3-t1)*10**(-6))
		Dt[1].append((t1-t_phone)*10**(-6))
		Dt[2].append((t2-t_phone)*10**(-6))
		Dt[3].append((t3-t_phone)*10**(-6))
		# stop sync
	sock_send.sendto(stop_command, (address, 5000))
	tend = time.time()
	print(tend-t0)
	sock_receive.close()
	
	duration = np.asarray(duration)
	for key in Dt.keys():
		Dt[key] = np.asarray(Dt[key])
	return Dt,duration
	
def get_lag(Dt,duration):
	
	tmedian = np.median(duration)
	tmax = tmedian*1.0

	indices = np.where(duration<tmax)[0]
	#print(Dt[2])
	print(indices)
	Dt = np.mean(np.asarray(Dt[2])[indices])
	
	return Dt

	
