import numpy as np
from pprint import pprint

import icewave.geometry.display as display
from icewave.geometry.define import *

import stephane.display.graphes as graphes


def save(table,filename):
    with open(filename,'w') as f:
        for line in table:
            line_s = map(str,line)
            #print(line_s)
            s = '\t'.join(line_s)
            s = s+'\n'
            f.write(s)
            
if __name__ == '__main__':
    pass
#    save()    


    
