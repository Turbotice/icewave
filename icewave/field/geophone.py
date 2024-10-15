
import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data
import glob
from pprint import pprint

global base
base = df.find_path(disk='Hublot24')

def get_records(date,year='2024'):
    files = glob.glob(base+date+'/Geophones/DigiSolo*.txt')#/*/*.srt')

    records = {}
    records['geophones'] = {}
    for filename in files:
        record = read_digiSolo(filename)
        num = filename.split('/DigiSolo_')[1].split('.txt')[0]
        #find only the data from this day
        select={}
        for key in record.keys():
            elem=record[key]
            if elem['date'].split('/')==[year,date[:2],date[2:]]:
                print(f'date matching for {num}')
                select[key]=elem
        #.keys()#.keys()#['04-waves_001'][0].keys()
        if len(select)>0:         
            records['geophones'][num] = select#record    
    return records

def read_digiSolo(filename):
    print(filename)
    geo = rw_data.read_csv(filename)
    record = {}
    for i,line in enumerate(geo[:-1]):
        l = geo[i+1]
        if len(l)>0:
            if 'Serial Number' in l[0]:
                serial = l[0].split(' = ')[1]
                #print(serial)
            if 'Start Acquisition FileName' in l[0]:
                rec = l[0].split('Seis')[1].split('.DLD')[0]
                rec = 'Seis'+rec #keep the full name as a key for the record dictionnary
                record[rec]={}
                record[rec]['date']=line[0].split('= "')[1]
                record[rec]['time']=[line[1][:-1]]
                record[rec]['latitude']=[]
                record[rec]['longitude']=[]
            elif len(line)>0:
                if 'RTC Time' in line[0] or 'UTC Time' in line[0]:
                    record[rec]['time'].append(line[1][:-1])
                elif 'Latitude' in line[0]:
                    lat = line[0].split(' = ')[1]
                    record[rec]['latitude'].append(float(lat))
                elif 'Longitude' in line[0]:
                    lon = line[0].split(' = ')[1]
                    record[rec]['longitude'].append(float(lon))
    
    return record    
