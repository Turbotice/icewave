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

def add_sources(Dlist,n0,x0,y0,z0,sign=1,axis=0):
    lines = []
    for i,D in enumerate(Dlist):
        if axis==0:
            lines.append(gen_point(n0+i,x0-D*sign,y0,z0,instrument='S'))
        if axis==1:
            lines.append(gen_point(n0+i,x0,y0-D*sign,z0,instrument='S'))
    return lines

def gen_sources(Dlist,L,n,n0=1,x0=0,y0=0,z0=0,axis=0,direction=1):
    lines = []

    lines1 = add_sources(Dlist,n0,x0,y0,z0,sign=1*direction,axis=axis)

    if axis==0:
        lines2 = add_sources(Dlist,n0+len(Dlist),x0+(n-1)*L,y0,z0,sign=-1*direction,axis=axis)
    if axis==1:
        lines2 = add_sources(Dlist,n0+len(Dlist),x0,y0+(n-1)*L*direction,z0,sign=-1*direction,axis=axis)
    lines = lines1+lines2

    return lines

def gen_S123(num,x0=0,y0=0,z0=0,axis=0,direction=1):
    L = 3
    n = 16
    lines = gen_sources([5,8,11],L,n,n0=num,x0=x0,y0=y0,z0=z0,axis=axis,direction=direction)
    return lines
    
def gen_line(n,L,n0=1,dn=1,x0=0,y0=0,z0=0,instrument='G',axis=0,direction=1):
    #create a line with n elements, spaced by L
    #by default, the line is oriented along x
    table = []
    header = ['#','X','Y','Z']
    table.append(header)
    for i,num in enumerate(range(n0,n*dn+n0,dn)):
        if axis==0:
            line = [instrument+'_'+ndigit(num),direction*i*L+x0,y0,z0]
        if axis==1:
            line = [instrument+'_'+ndigit(num),x0,direction*i*L+y0,z0]

        table.append(line)
    return table

def gen_arccircle(n,R,dtheta,xc,yc,x0,y0):

    X,Y = Smath.pol2cart(R,theta)

def gen_rectangle(nx,ny,Lx,Ly,x0=0,y0=0,z0=0,n0=1,directionx=1,directiony=1,instrument='G'):
    table = []
    header = ['#','X','Y','Z']
    table.append(header)

    num = n0
    for j in range(ny):
        for i in range(nx):
            line = [instrument+'_'+ndigit(num),x0+Lx*i*directionx,y0+Ly*j*directiony,z0]
            num=num+1
            table.append(line)
    return table   

def ndigit(i,n=2):
    n0 = len(str(i))
    if n0<n:#only works for 2 digits
        return '0'+str(i)
    else:
        return str(i)
    
def add_lines(table,lines,n0=1):
    #sort the lines by x position ??
    for line in lines[n0:]:
        table.append(line)
    return table

def exemple1():
    n=16
    L=5
    D=10
    table = gen_line(n,L,instrument='G')
    lines = gen_sources([D],L,n)
    table = add_lines(table,lines)
    return table

def exemple2():
    n=4
    ni = n-1
    L=20
    D=20
    table = gen_rectangle(n,n,L,L,instrument='G')

    table = add_lines(table,gen_sources([D],L,n))
    table = add_lines(table,gen_sources([D],L,n,n0=3,y0=ni*L))
    table = add_lines(table,gen_sources([D],L,n,n0=5,axis=1))
    table = add_lines(table,gen_sources([D],L,n,n0=7,axis=1,x0=ni*L))
    return table

if __name__ == '__main__':
    table = sag24_0210_Geophone_line1()
    display.show(table)
