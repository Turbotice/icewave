import stephane.display.graphes as graphes
import stephane.tools.Smath as smath

import icewave.tools.datafolders as df
import icewave.tools.rw_data as rw_data

#import sympy #symoblic python
#import mpmath as math
#import cv2
import glob
import csv
import os
import numpy as np

#import icewave.phone.rw_pyphone as rw
import icewave.tools.rw_data as rw_data

import icewave.field.drone as drone
import icewave.field.phone as phone
import icewave.field.geophone as geophone
import icewave.field.buoys as buoys
import icewave.field.gps as gps


import argparse

def gen_parser():    
    parser = argparse.ArgumentParser(description="Manipulate multi instruments data")
    parser.add_argument('-date', dest='date', type=str,default='0226',help='select date to process data')
    parser.add_argument('-year', dest='year', type=str,default='2025',help='select year to process data')

    #parser.add_argument('-step', dest='step', type=int,default=3,help='select Step to be performed')

#    print(parser)   
    args = parser.parse_args()
    print(args)
    return args

def get_records(date,year='2025'):
    records = gps.get_records(date,year=year)
    records.update(drone.get_records(date,year=year))
    records.update(phone.get_records(date,year=year))
    records.update(geophone.get_records(date,year=year))
    #records.update(buoys.get_records(date,year))
    return records

def save_records(date,year):
    records = get_records(date)
    base = df.find_path(year)
    filename = base+date+'/Summary/records_'+date+'.pkl'

    folder = os.path.dirname(filename)
    print(folder)
    if not os.path.exists(folder):
        os.makedirs(folder)
    rw_data.write_pkl(filename,records)
    return filename

def compute_timeline(date,date_format='2024/02/03',year='2025'):
    filename = save_records(date,year)
    records = rw_data.load_pkl(filename)

    savefolder= os.path.dirname(filename)
    print(records.keys())
    fig,ax,figs = display_timeline(records,date,date_format=date_format)
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    #ax.set_xlim(['11:30:00','15:00:00'])
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%h-%m-%s'))
    #fig.autofmt_xdate()
    graphes.save_figs(figs,savedir=savefolder,overwrite=True,prefix=date+'_')
    #graphes.save_figs(figs,savedir=savefolder_local,overwrite=True,prefix=date+'_')

def display_timeline(records,date,date_format='2024/02/23',year='2025'):
    fig = plt.figure(figsize=(15,10))
    ax = fig.add_subplot(111)

    positions = {'Fulmar':1,'mesange':1,'Bernache':2,'geophones':2.6,'buoys':3.8,'phones':4.5}
        
    import matplotlib.dates as mdates    
    bs = {'Fulmar':0,'mesange':1,'Bernache':1}
    
    record = records['drones']
    for (drone,name) in selection:
    #record = records['drones'][k]
        pos = positions[drone]
        try:
            title = name.split('-')[1]#.split('_')[0]
        except:
            title = name
        if name in record[drone].keys():
            times = record[drone][name][0]['time']# for key in record[name][0].keys()]
            x = [times[0],times[-1]]
            x = np.asarray([multi.convert_time(x[0]),multi.convert_time(x[1])])

            r=np.random.random()*bs[drone]-0.5
            y = pos+r
            ax.plot(x/3600, [y,y],colors[drone]+'-')
            ax.text(x[0]/3600,y+0.03,title)

    key = 'geophones'
    record = records[key]

    pos = positions[key]
    X={}
    first={}
    for i,num in enumerate(record.keys()):
        first[num]=True
        for k in record[num].keys():
#            date_format = str("2024/"+date[:2]+"/"+date[2:])#'2024/02/23'
            if record[num][k]['date']==date_format:# 2024/02/23date_format:
                t0 = record[num][k]['time'][0]
                t1 = record[num][k]['time'][-1]
                x = [t0,t1]
                x = np.asarray([multi.convert_time(x[0]),multi.convert_time(x[1])])

                y = pos+i*0.05#np.random.random()*b
                ax.plot(x/3600, [y,y],colors[key]+'-')
                if first[num]:
                    ax.text(x[0]/3600-0.25,y,num)
                    first[num]=False

    x0 = 0.1
    for key in ['buoys','phones']:
        record = records[key]
        pos = positions[key]
        X={}
        first={}
        for i,num in enumerate(record.keys()):
            first[num]=True
            for k in record[num].keys():
                elem = record[num][k]
                t0 = elem['time'][0]
                t1 = elem['time'][-1]
                x = [t0,t1]
                X[i] = np.asarray([multi.convert_time(x[0]),multi.convert_time(x[1])])
                #np.asarray([convert(t0)-3600*6,convert(t1)-3600*6])

                print(i,X[i][0]/3600)
                if (X[i][0]/3600)<16:
                    if i>0:
                        X[i]=X[i-1]
                    else:
                        continue
                y = pos+i*0.1#np.random.random()*b
                ax.plot(X[i]/3600, [y,y],colors[key]+'-')
                if first[num]:
                    ax.text(X[i][0]/3600-x0,y,num)
                    first[num]=False
                    
    ax.set_ylim([0,6])

    hours = np.asarray([11,12,13,14,15])+6
    ax.set_xticks(hours,[str(h)+':00' for h in hours])

    # access legend objects automatically created from data
    handles, labels = plt.gca().get_legend_handles_labels()

    legends = timeline_legend()
    # add manual symbols to auto legend
    handles.extend(legends)
    plt.legend(handles=handles,fontsize=20)
    figs = graphes.legende('UTC time','Instruments',date)
    #plt.show()
    return fig,ax,figs

def timeline_legend():
    from matplotlib.lines import Line2D

    names = {'Smartphones':('r','s'),'Buoys':('m','o'),'Geophones':('g','>'),'Bernache':('b','o'),'Mesange':('k','o')}
    legends = []
    for key in names.keys():
        point = Line2D([0], [0], label=key, marker=names[key][1], markersize=20, 
             markeredgecolor=names[key][0], markerfacecolor=names[key][0], linestyle='')
        legends.append(point)
    return legends
#base = '/media/turbots/Hublot24/Share_hublot/Data/'

def convert_time(t):
    h,m,s = t.split(':')
    return int(h)*3600+int(m)*60+int(s)

def get_time_interval(record):
    times = [record[key]['time'] for key in record.keys()]
    t0 = times[0]
    t1 = times[-1]
    return t0,t1

def get_avg_position(record):
    latitudes = [record[key]['latitude'] for key in record.keys()]
    longitudes = [record[key]['longitude'] for key in record.keys()]
    latitude = np.mean(latitudes)
    longitude = np.mean(longitudes)
    return latitude,longitude

def read_BA_timeline(date):
    base = df.find_path()
    filename = base + f'Summary/{date}_path_drone.txt'
    
    filelist = rw_data.read_csv(filename)
    keydrones = [f[0].split('/')[-1] for f in filelist]
    drones = [f[0].split('/')[-3] for f in filelist]

    out = [(drone,keydrone) for (drone,keydrone) in zip(drones,keydrones)]
    return out

def main(args):
    if args.date=='all':
        dates = ['0210','0211','0220','0221','0223','0226','0306']
        for date in dates:
            save_records(date,year=args.year)
    else:
        save_records(args.date,year=args.year)
    
if __name__ =='__main__':
    args = gen_parser()
    main(args)
