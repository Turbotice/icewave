


import subprocess


def get_phone():
    pass

def adb_command(phone):
    exemple = ['adb','-s','192.168.0.100','shell','dumpsys','battery']
    output = subprocess.run(exemple,stdout=subprocess.PIPE)
