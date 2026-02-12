import numpy as np
from pprint import pprint


import icewave.geometry.save as save
import icewave.geometry.display as disp
from icewave.geometry.define import *

import pylab as plt
import stephane.display.graphes as graphes

def Geophones_line1():
    n=16
    L=3
    D = 1
    table = gen_line(n,L,n0=101,instrument='G')
    table = add_lines(table,gen_S123(101,x0=0,y0=0,z0=0,axis=0),n0=0)
    table = add_lines(table,gen_S123(107,x0=0,y0=0,z0=0,axis=0),n0=0)
    return table

def Geophones_pyramide1():
    L = 15

    table = gen_line(1,0,n0=201,instrument='G')
    for i,num in enumerate([5,8,10]):
        table.append(gen_point(200+num,(i+1)*L,0,0,instrument='G'))

    for i,num in enumerate([2,6,9]):
        table.append(gen_point(200+num,7.5+i*L,-13,0,instrument='G'))

    for i,num in enumerate([11,13,16]):
        table.append(gen_point(200+num,7.5+i*L,13,0,instrument='G'))

    for i,num in enumerate([3,7]):
        table.append(gen_point(200+num,15+i*L,-13*2,0,instrument='G'))

    for i,num in enumerate([12,15]):
        table.append(gen_point(200+num,15+i*L,13*2,0,instrument='G'))

    table.append(gen_point(204,22.5,-13*3,0,instrument='G'))
    table.append(gen_point(214,22.5,13*3,0,instrument='G'))

#    print(table)

    return table


def Telephones_1():
    table = []
    x0 = 45+20
    table = gen_line(1,0,x0=x0,n0=119,instrument='T',axis=0,direction=1)
    table.append(gen_point(100+18,10+x0,0,0,instrument='T'))
    table.append(gen_point(100+17,20+x0,0,0,instrument='T'))
    table.append(gen_point(100+16,30+x0,0,0,instrument='T'))
    table.append(gen_point(100+13,40+x0,0,0,instrument='T'))
    table.append(gen_point(100+11,50+x0,0,0,instrument='T'))
    table.append(gen_point(100+9,60+x0,0,0,instrument='T'))
    table.append(gen_point(100+6,70+x0,0,0,instrument='T'))
    table.append(gen_point(100+4,80+x0,0,0,instrument='T'))
    table.append(gen_point(100,90+x0,0,0,instrument='T'))
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

    table.append(gen_point(100+1,65,0,0,instrument='B'))
    table.append(gen_point(100+2,90,0,0,instrument='B'))
    table.append(gen_point(100+3,115,0,0,instrument='B'))
    table.append(gen_point(100+5,140,0,0,instrument='B'))

    return table

def Buoys_2():
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    table.append(gen_point(100+2,0,-22.5,0,instrument='B'))
    table.append(gen_point(100+5,22.5,-22.5,0,instrument='B'))
    return table


def Sag24(display=True):
    table1 = Geophones_line1()
    figs={}
    #ax,figs = display.show(table1,sx=10,sy=2,display=False)
 #   figs.update(fig)
    
    table2 = Telephones_1()

    table3 = Buoys_1()

    table4 = Geophones_pyramide1()

    #table8 = Geophones_Quinconce1()
    
    tables = table1
#    tables = add_lines(tables,table1,n0=1)
    tables = add_lines(tables,table2)
    tables = add_lines(tables,table3)
    tables = add_lines(tables,table4)
    #tables = add_lines(tables,table6,n0=1)
    #tables = add_lines(tables,table7,n0=1)
    #tables = add_lines(tables,table8,n0=1)

    tables = save.form(tables)

    if display:
        ax,figs = disp.show(tables,sx=20,sy=10,display=False)
    else:
        figs=None
    #figs.update(fig)

    
    return figs,tables
    
if __name__ == '__main__':
    Sag24()    
    #plt.show()
    
