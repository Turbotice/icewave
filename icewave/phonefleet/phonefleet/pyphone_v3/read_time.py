import glob	


folder = '/run/user/1000/gvfs/afp-volume:host=thiou.local,user=sperrard,volume=labshared2/Telephones/Time/timelog_181223/'

filelist = glob.glob(folder+'*.txt')
print(filelist)

tref = 0.04898786544799805
for filename in filelist:
    f=open(filename,'r')
    num = filename.split('/')[-1].split('.txt')[0][1:]
    t1 = float(f.readline().split('\n')[0].replace(',','.'))
    t2 = float(f.readline().split('\n')[0])
    t3 = float(f.readline().split('\n')[0].replace(',','.'))
    
    print(num,t2-t1-tref)
