
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

global base
base = df.find_path(disk='Hublot24')

def get_records(date):
    srtfiles = get_srtfiles(date)
    records = {}
    for key in srtfiles.keys():
        records[key]={}
        for i,srtfile in enumerate(srtfiles[key]):
            print(i,srtfile)
            name = srtfile.split('/')[-2]
            record = get_flighrecord(srtfile,drone=key)
            records[key][name]=record
    return records
    
def get_srtfiles(date):
    srtfiles = {}
    print(base)
    drones = ['mesange','Bernache','Fulmar']
    for key in drones:
        srt = glob.glob(base+date+'/Drones/'+key+'/*/*.SRT')#/*/*.srt')
        pprint(srt)
        if len(srt)>0:
            srtfiles[key] = srt
        else:
            print(f"No data for {key} on {date}")
    return srtfiles

def get_flighrecord(srtfile,step=100,drone='mesange'):
    #convert all times to UTC
    if drone=='mesange':
        h0 = -1
    elif drone=='Bernache':
        h0 = 5
    elif drone=='Fulmar':
        print('Time to be checked')
        h0 = 0
    else:
        h0=0
        print(drone)
        print('Drone unknown')
        
    data = rw_data.read_csv(srtfile)
    n = int(len(data)/6)
    print('number of records : '+str(n))
    record = {}
    for i in range(0,n-1,step):
        event = data[i*6:(i+1)*6]
        if int(event[0][0])==i+1:  
            record[i]={}
            record[i]['record_time']=event[1]
            record[i]['date']=event[3][0].split(' ')[0]

            time = event[3][0].split(' ')[1][:-4]
            h,m,s = time.split(':')
            hnew = str(int(h)+h0)
            t = f"{hnew}:{m}:{s}"            
            record[i]['time']=t
            params = event[4][0]
            latitude = float(params.split('latitude: ')[1].split(']')[0])
            longitude = float(params.split('longitude: ')[1].split(']')[0])

            #print(event[3],latitude,longitude)
            record[i]['latitude']=latitude
            record[i]['longitude']=longitude
            record[i]['params'] = params
#pprint(d[6:12])
    return record
    
