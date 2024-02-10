import numpy as np
from pprint import pprint

import icewave.geometry.display as display
from icewave.geometry.define import *


def exemple1():
    n=16
    L=5
    D=10
    table = gen_line(n,L,instrument='G')
    lines = gen_sources(D,L,n)
    table = add_lines(table,lines)
    return table

def exemple2():
    n=4
    ni = n-1
    L=20
    D=20
    table = gen_rectangle(n,n,L,L,instrument='G')

    table = add_lines(table,gen_sources(D,L,n))
    table = add_lines(table,gen_sources(D,L,n,n0=3,y0=ni*L))
    table = add_lines(table,gen_sources(D,L,n,n0=5,axis=1))
    table = add_lines(table,gen_sources(D,L,n,n0=7,axis=1,x0=ni*L))
    return table

def sag24_0206_Geophone1():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=1,instrument='G')
    table.append(gen_point(1,-1,0,0,instrument='S'))
    table.append(gen_point(2,15,0,0,instrument='S'))
    table.append(gen_point(3,33,0,0,instrument='S'))
    table.append(gen_point(4,46,0,0,instrument='S'))
    return table

def sag24_0206_Geophone2():
    n=8
    L=6
    D = 1
    table = gen_line(n,L,n0=1,dn=2,instrument='G')
    table = add_lines(table,gen_line(n,L,n0=2,dn=2,x0=3,y0=-5.2,instrument='G')[1:])

#    table.append(gen_point(1,-1,0,0,instrument='S'))
#    table.append(gen_point(2,15,0,0,instrument='S'))
#    table.append(gen_point(3,33,0,0,instrument='S'))
#    table.append(gen_point(4,46,0,0,instrument='S'))
    return table    
    
    
def sag24_0210_G_Line1():
    n=16
    L=3
    D=1
    table = gen_line(n,L,n0=1,instrument='G')


    table.append(gen_point(1,-1,0,0,instrument='S'))
    table.append(gen_point(2,15,0,0,instrument='S'))
    table.append(gen_point(3,33,0,0,instrument='S'))
    table.append(gen_point(4,46,0,0,instrument='S'))
    return table


if __name__ == '__main__':
    table = exemple2()
    display.show(table)
