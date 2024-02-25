import numpy as np
from pprint import pprint

import icewave.geometry.display as display
from icewave.geometry.define import *

import stephane.display.graphes as graphes

def form(table):
    for i,line in enumerate(table):
        if i>0:#do not modify for the header
            tag,num = line[0].split('_')
            if len(num)>2:
                num = int(num)
                num = f"{num:04d}"
                line[0] = tag+'_'+num
                table[i] = line
    return table


def save(table,filename):
    table = form(table)
    with open(filename,'w') as f:
        for line in table:
            print(table[0])
            #add missing zeros
            
            line_s = map(str,line)
            s = '\t'.join(line_s)
            print(s)
            s = s+'\n'
            f.write(s)
            
if __name__ == '__main__':
    pass
#    save()    


    
