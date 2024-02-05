import numpy as np
from pprint import pprint

import icewave.geometry.display as display

instruments = {'G':'Geophones',\
               'T':'Telephones',\
               'B':'Buoys',\
               'D':'Drone',\
               'C':'Camera',\
               'S':'Source',\
               'X':'GPS',\
               'CTD':'Conductivity, Temperature, Depth',\
               'H':'Ice Thickness'}
                
def gen_point(num,x0,y0,z0,instrument='G'):
    return [instrument+'_'+ndigit(num),x0,y0,z0]

def gen_sources(D,L,n,n0=1,x0=0,y0=0,z0=0,axis=0):
    lines = []
    if axis==0:
        lines.append(gen_point(n0,x0-D,y0,z0,instrument='S'))
        lines.append(gen_point(n0+1,x0+(n-1)*L+D,y0,z0,instrument='S'))
    if axis==1:
        lines.append(gen_point(n0,x0,y0-D,z0,instrument='S'))
        lines.append(gen_point(n0+1,x0,y0+(n-1)*L+D,z0,instrument='S'))
    return lines

def gen_line(n,L,x0=0,y0=0,z0=0,instrument='G'):
    #create a line with n elements, spaced by L
    #by default, the line is oriented along x
    table = []
    header = ['#','X','Y','Z']
    table.append(header)
    for i in range(n):
        line = [instrument+'_'+ndigit(i),i*L+x0,y0,z0]
        table.append(line)
    return table

def gen_rectangle(nx,ny,Lx,Ly,x0=0,y0=0,z0=0,instrument='G'):
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    num=0
    for j in range(ny):
        for i in range(nx):
            num=num+1
            line = [instrument+'_'+ndigit(num),x0+Lx*i,y0+Ly*j,z0]
            table.append(line)
    return table   

def ndigit(i,n=2):
    n0 = len(str(i))
    if n0<n:#only works for 2 digits
        return '0'+str(i)
    else:
        return str(i)
    
def add_lines(table,lines):
    #sort the lines by x position ??
    for line in lines:
        table.append(line)
    return table

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


if __name__ == '__main__':
    table = exemple2()
    display.show(table)