import tsync


global phonelist
phonelist = []

global network


folder = '/storage/emulated/0/Documents/Gobannos/'

def set_network():
    s = input('Type sub network number')
    network = int(s)

def get_phonelist():
    s = input('Type phonelist to adress :')
    lis = s.split(':')

    for i,l in enumerate(lis[:-1]):
        phonelist = phonelist+list(range(int(lis[i]),int(lis[i+1])+1))
    phonelist = list(set(phonelist))

def check_status(phonelist):
    for phone in phonelist:
        get_status(phone)

def get_status(phone):
    connect.ipbase


def choose():
    s = input('Chose action among :')
    actions = ['network','phones','status','time','start','stop','exit']
    descriptions = ['','','','','','','']
    for action,description in zip(actions,descriptions):
        print(action,description)

    if s=='network':
        set_network()
    if s=='phones':
        get_phonelist()
    elif s=='status':
        print(f'Available phonelist : {phonelist}')
        get_status(phonelist)
    elif s=='time':
        time_sync(phonelist)
    elif s=='start':
        start(phonelist)
    elif s=='stop':
        stop(phonelist)
    elif s=='exit':
        print("exit")
    else:
        print('command not known, do nothing')
    return s

def main():

    s=''
    while not s=='exit':
        s = choose()
    #action to code :
    #define phonelist (phones)
    #check connection (status)
    #timesync (time)
    #start acquisition (start)
    #

if __name__=='__main__':
    main()
