import numpy as np
import subprocess


def run():
	filename = '/home/turbots/Documents/Codes/Python/pyphone_v2/adb_usb.sh'
	#filename = 'adb_usb.sh'
	print('run adb connect')
	c = subprocess.run(['bash',filename],text=True,capture_output = True,shell=False)
	log = c.stdout.split('\n')
	for l in log:
	    print(l)
	print('done')
