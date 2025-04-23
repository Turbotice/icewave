import pylab as plt
import numpy as np

import icewave.field.multi_instruments as multi
import icewave.geometry.display as display
import icewave.display.graphes as graphes

global colors
colors = display.colors

def display_records(records,date,date_format='2024/02/26',ax=None):
    if ax==None:
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)

    positions = {'Fulmar':1,'mesange':1,'Bernache':2,'geophones':2.6,'buoys':3.8,'phones':4.5}
        
    import matplotlib.dates as mdates
    h0 = 0
    
    bs = {'Fulmar':0,'mesange':1,'Bernache':1}
    
    record = records['drones']
    for drone in records['drones'].keys():
        for name in records['drones'][drone].keys():#(drone,name) in selection:
            record = records['drones']#[drone]
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
                ax.plot(x/3600+h0, [y,y],colors[drone]+'-')
                ax.text(x[0]/3600+h0,y+0.03,title)

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
                ax.plot(x/3600+h0, [y,y],colors[key]+'-')
                if first[num]:
                    ax.text(x[0]/3600-0.25+h0,y,num)
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

                #print(i,X[i][0]/3600)
                if (X[i][0]/3600)<16:
                    if i>0:
                        X[i]=X[i-1]
                    else:
                        continue
                y = pos+i*0.1#np.random.random()*b
                ax.plot(X[i]/3600+h0, [y,y],colors[key]+'-')
                if first[num]:
                    ax.text(X[i][0]/3600-x0+h0,y,num)
                    first[num]=False
                    
    ax.set_ylim([0,6])

    hours = np.asarray([11,12,13,14,15])+6
    ax.set_xticks(hours,[str(h)+':00' for h in hours])

    from matplotlib.lines import Line2D
    # access legend objects automatically created from data
    handles, labels = plt.gca().get_legend_handles_labels()
    names = {'Smartphones':('r','s'),'Buoys':('m','o'),'Geophones':('g','>'),'Bernache':('b','o'),'Mesange':('k','o')}
    legends = []
    for key in names.keys():
        point = Line2D([0], [0], label=key, marker=names[key][1], markersize=20, 
             markeredgecolor=names[key][0], markerfacecolor=names[key][0], linestyle='')
        legends.append(point)

    # add manual symbols to auto legend
    handles.extend(legends)

    #ax.legend(handles=handles,fontsize=20)
    figs = graphes.legende('UTC time','Instruments',date,ax=ax)
    #plt.show()
    return figs
