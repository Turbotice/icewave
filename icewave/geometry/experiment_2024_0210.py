import numpy as np
from pprint import pprint

import icewave.geometry.display as display
import icewave.geometry.save as save
from icewave.geometry.define import *

import stephane.display.graphes as graphes

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
    

def sag24_0210_Geophone_line1():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=101,instrument='G')
    table = add_lines(table,gen_S123(101,x0=0,y0=0,z0=0,axis=0))
    print(table)
    return table

def sag24_0210_Geophone_line2():
    n=16
    L=3
    D=1
    table = gen_line(n,L,n0=201,instrument='G',axis=1,direction=-1)
    table = add_lines(table,gen_S123(201,x0=0,y0=0,z0=0,axis=1,direction=-1))
    print(table)
    return table

def sag24_0210_Geophone_line3():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=301,y0=-45,instrument='G')
    table = add_lines(table,gen_S123(301,x0=0,y0=-45,z0=0,axis=0))
    print(table)
    return table

def sag24_0210_Telephones():
    table = []
    header = ['#','X','Y','Z']
    table.append(header)
    
    phonelist = [1,6,7,13,16,17,18]
    for phone in phonelist:
        table.append(gen_point(100+phone,0,0,0,instrument='T'))
    table.append(gen_point(100+12,0,-22.5,0,instrument='T'))
    table.append(gen_point(100+4,22.5,-22.5,0,instrument='T'))
    table.append(gen_point(100+8,0,-45,0,instrument='T'))
    table.append(gen_point(100+9,22.5,-45,0,instrument='T'))
    table.append(gen_point(100+11,45,-45,0,instrument='T'))
    return table    

def sag24_0210_Buoys():
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    table.append(gen_point(100+2,0,-22.5,0,instrument='B'))
    table.append(gen_point(100+5,22.5,-22.5,0,instrument='B'))
    return table

def sag24_0210_Tomo():
    L = 30
    table = gen_rectangle(4,4,L,L,n0=401,x0=0,y0=0,z0=0,directionx=1,directiony=-1,instrument='G')
    return table


def sag24_0210():
    table1 = sag24_0210_Geophone_line1()
    figs={}
    #ax,figs = display.show(table1,sx=10,sy=2,display=False)
 #   figs.update(fig)
    
    table2 = sag24_0210_Geophone_line2()
    #ax,figs = display.show(table2,sx=2,sy=10,display=False)
  #  figs.update(fig)
    
    table3 = sag24_0210_Geophone_line3()
    #ax,figs = display.show(table3,sx=10,sy=2,display=False)
   # figs.update(fig)

    table4 = sag24_0210_Telephones()

    table5 = sag24_0210_Buoys()

    table6 = sag24_0210_Tomo()
    
    tables = table1
#    tables = add_lines(tables,table1,n0=1)
    tables = add_lines(tables,table2,n0=1)
    tables = add_lines(tables,table3,n0=1)
    tables = add_lines(tables,table4,n0=1)
    tables = add_lines(tables,table5,n0=1)
    tables = add_lines(tables,table6,n0=1)

    ax,figs = display.show(tables,sx=10,sy=10,display=False)
    #figs.update(fig)

    return figs,tables
    
if __name__ == '__main__':
    figs,tables = sag24_0210()

    filename= 
    save.save(table,filename)

    
