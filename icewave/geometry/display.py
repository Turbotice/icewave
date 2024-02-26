

import pylab as plt
import numpy as np
import glob
from pprint import pprint

import stephane.display.graphes as graphes

global labels

labels = {  'G':'g>',\
            'T':'rs',\
            'B':'mo',\
            'D':'kp',\
            'C':'gs',\
            'S':'bo',\
            'X':'m+',\
            'P':'r*',\
            'CTD':'b>',\
            'H':'bv'}

def show(table,dim=2,eps=1,sx=6,sy=6,display=True):

    norme = read_norme()
    n = len(table)

    fig,ax = plt.subplots(figsize=(sx,sy))
    notes = {}
    notes['pos']=[]
    notes['label']=[]
    c=1
    for i in range(1,n):
        tag = table[i][0]
        typ,num = tag.split('_')
        label = norme[typ]['label']
        x = table[i][1]
        y = table[i][2]

        b=0
        #print(norme[typ]['size'])
#        notes.append([x,y+eps])
        if i>2:
            d = np.linalg.norm(np.asarray(notes['pos'][:-1])-np.asarray([x,y+eps]),axis=1)
            if np.min(d)==0:
                j = np.argmin(d)
                c=4
                print(label,notes['label'][j])
                if label==notes['label'][j]:
                    b=0
                else:
                    b=1
            else:
                c=1
                b=0
                #print(np.min(d),c,label)

        ax.plot(x,y+b*eps,label,markersize=norme[typ]['size'])
        tagd = '$'+typ+'^{'+num[0]+'}_{'+str(int(num[1:]))+'}$'#.replace('_','')
        ax.annotate(tagd,(x,y+c*eps),fontsize=20)
        notes['pos'].append([x,y+c*eps])
        notes['label'].append(label)
        
    title = ''
    figs = graphes.legende('$X$ (m)','$Y$ (m)',title)
    plt.axis('equal')
    if display:
        plt.show()
    
    return ax,figs


def read_norme():
#    print(glob.glob('*'))
    filename = "/Users/stephane/Documents/git/icewave/icewave/display/Nomenclature_plots.txt"
    with open(filename,'r') as f:
        out = f.read()
    lines = out.split('\n')
    header = lines[0].split('\t')
    dtable = {}
    for line in lines[1:]:
        line = line.split('\t')
        key = line[0]
        dtable[key]={}
        for head,item in zip(header,line):
            dtable[key][head]=item
#    table = np.asarray([line.split('\t') for line in lines])
#    dtable = {tab[0]:tab[1] for tab in table}
#    pprint(dtable)
    return dtable
