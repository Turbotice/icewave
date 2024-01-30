import numpy as np

instruments = {'G':'Geophones',\\
               'T':'Telephones',\\
               'B':'Buoys',\\
               'D':'Drone',\\
               'C':'Camera',\\
               'S':'Source',\\
               'CTD':'Conductivity, Temperature, Depth',\\
               'H':'Ice Thickness'}
                

def line(n,L,x0=0,y0=0,z0=0,instrument='G'):
    #create a line with n elements, spaced by L
    #by default, the line is oriented along x

    buffer = []
    header = ['#','X','Y','Z']
    for i in range(n):
        line = [instrument+2digit(i),i*L+x0,y0,z0]
        buffer.append(line)
    return line

def 2digit(i,n=2):
    n0 = len(str(i))
    if n0<n:#only works for 2 digits
        return '0'+str(i)
    else:
        return str(i)
    
def add_point(num,x0=0,y0=0,z0=0,instrument='G'):
    return [instrument+2digit(i),x0,y0,z0]
