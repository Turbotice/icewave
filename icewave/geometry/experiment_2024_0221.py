import numpy as np
from pprint import pprint

import icewave.geometry.display as display
from icewave.geometry.define import *
import pylab as plt
import stephane.display.graphes as graphes

def Geophones_line1():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=201,instrument='G')
    table = add_lines(table,gen_S123(201,x0=0,y0=0,z0=0,axis=0))
    print(table)
    return table

def Geophones_line2():
    n=16
    L=3
    D=1
    table = gen_line(n,L,n0=101,instrument='G',axis=1,direction=-1)
    table = add_lines(table,gen_S123(101,x0=0,y0=0,z0=0,axis=1,direction=-1))
    print(table)
    return table

def Geophones_line3():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=301,y0=-45,instrument='G')
    table = add_lines(table,gen_S123(301,x0=0,y0=-45,z0=0,axis=0))
    print(table)
    return table

def Geophones_Tomo1():
    Lx = 22.5
    Ly = 22.5
    table = gen_rectangle(4,4,Lx,Ly,n0=401,x0=-22.5,y0=22.5,z0=0,directionx=1,directiony=-1,instrument='G')
    return table

def Geophones_Quinconce1():
    n=6
    L=30
    table = gen_line(6,L,n0=501,x0=110,y0=-100,instrument='G')
    table = add_lines(gen_line(7,L,n0=507,x0=110+26,y0=-100,instrument='G'))
#    table = add_lines(table,gen_S123(301,x0=0,y0=-45,z0=0,axis=0))
    return table

def Telephones_1():
    table = []
    n = 46
    L=1

    table = gen_line(n,L,n0=111,instrument='T',axis=1,direction=-1)
    table.append(gen_point(100+11,0,0,0,instrument='T'))


#    phonelist = [1,6,7,13,16,17,18]
#    for phone in phonelist:
#        table.append(gen_point(100+phone,0,0,0,instrument='T'))
    #table.append(gen_point(100+11,0,-22.5,0,instrument='T'))
    #table.append(gen_point(100+12,22.5,-22.5,0,instrument='T'))
    #table.append(gen_point(100+13,45,-22.5,0,instrument='T'))
    #table.append(gen_point(100+16,22.5,-45,0,instrument='T'))
    #table.append(gen_point(100+21,22.5,0,0,instrument='T'))
    return table

def Buoys_1():
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    table.append(gen_point(100+2,0,-22.5,0,instrument='B'))
    table.append(gen_point(100+5,22.5,-22.5,0,instrument='B'))
    return table

def Buoys_2():
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    table.append(gen_point(100+2,0,-22.5,0,instrument='B'))
    table.append(gen_point(100+5,22.5,-22.5,0,instrument='B'))
    return table


def Sag24_0221():
    table1 = Geophones_line1()
    figs={}
    #ax,figs = display.show(table1,sx=10,sy=2,display=False)
 #   figs.update(fig)
    
    table2 = Geophones_line2()
    #ax,figs = display.show(table2,sx=2,sy=10,display=False)
  #  figs.update(fig)
    
    table3 = Geophones_line3()
    #ax,figs = display.show(table3,sx=10,sy=2,display=False)
   # figs.update(fig)

    table4 = Telephones_1()

    table5 = Geophones_Tomo1()

    #table8 = Geophones_Quinconce1()
    
    tables = table1
#    tables = add_lines(tables,table1,n0=1)
    tables = add_lines(tables,table2,n0=1)
    tables = add_lines(tables,table3,n0=1)
    tables = add_lines(tables,table4,n0=1)
    tables = add_lines(tables,table5,n0=1)
    #tables = add_lines(tables,table6,n0=1)
    #tables = add_lines(tables,table7,n0=1)
    #tables = add_lines(tables,table8,n0=1)

    ax,figs = display.show(tables,sx=20,sy=20,display=False)
    #figs.update(fig)

    return figs,tables
    
if __name__ == '__main__':
    Sag24_0221()    
    plt.show()
    
