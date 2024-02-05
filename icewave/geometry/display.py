

import pylab as plt
import numpy as numpy
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

def show(table,dim=2,eps=1):
    n = len(table)

    fig,ax = plt.subplots(figsize=(6,6))

    for i in range(1,n):
        tag = table[i][0]
        label,num = tag.split('_')
        x = table[i][1]
        y = table[i][2]

        ax.plot(x,y,labels[label])
        ax.annotate(tag,(x,y+eps))

    title = ''
    graphes.legende('$X$ (m)','$Y$ (m)',title)
    plt.axis('equal')
    plt.show()
    