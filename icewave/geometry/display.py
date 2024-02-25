

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

    for i in range(1,n):
        tag = table[i][0]
        typ,num = tag.split('_')
        label = norme[typ]['label']
        x = table[i][1]
        y = table[i][2]

        ax.plot(x,y,label,markersize=norme[typ]['size'])
        tagd = '$'+typ+'_{'+num+'}$'#.replace('_','')
        ax.annotate(tagd,(x,y+eps))

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
