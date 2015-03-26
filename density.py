#/usr/bin/env python
#-*- encoding: utf-8 -*-
#möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import sys
from misc import vol_unitcell
import math
import cmath
import numbers
from math import sqrt


'''
Created on 26.03.2015

@author: Daniel Kratzert
'''


def cell_from_fcf(filename):
    '''
    reads unit cell from fcf
    '''
    a, b, c, al, be, ga, list_code = None, None, None, None, None, None, None  
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('_cell_length_a'):
                    a = float(line.split()[1])
                if line.startswith('_cell_length_b'):
                    b = float(line.split()[1])
                if line.startswith('_cell_length_c'):
                    c = float(line.split()[1])
                if line.startswith('_cell_angle_alpha'):
                    al = float(line.split()[1])
                if line.startswith('_cell_angle_beta'):
                    be = float(line.split()[1])
                if line.startswith('_cell_angle_gamma'):
                    ga = float(line.split()[1])
                if line.startswith('_shelx_refln_list_code'):
                    list_code = line.split()[1]
                if all([a, b, c, al, be, ga, list_code]):
                    break
    except(IOError):
        print('Unable to read {} file.'.format(filename))
        sys.exit()
    return (a, b, c, al, be, ga, list_code)
                

def read_data_from_fcf(filename, type):
    '''
    read Fo and Fc data from fcf
    '''
    data = []
    read_data = False
    count = 0
    try:
        with open(filename, 'r') as f:
            for line in f:
                if read_data:
                    line = line.split()
                    try:
                        line[1]
                    except:
                        continue
                    if type == '6':
                        # h k l Fc^2 Fo^2
                        line = [int(line[0]), int(line[1]), int(line[2]), 
                                float(line[5]), sqrt(abs(float(line[3])))]
                        data.append(line)
                    if type == '4':
                        data.append([int(line[0]), int(line[1]), int(line[2]), 
                                    float(line[3]), float(line[4])])
                    continue
                if line.startswith(' _refln'):
                    count = count+1
                    if count == 7:
                        read_data = True
    except(IOError):
        print('Unable to read {} file.'.format(filename))
        sys.exit()
    return data
    
def fourier(hklfcfo, xyz, c, diffden=False):
    '''
    rho = 1/V*SUM([Fo-Fc]*exp(-i2pi(hx+ky+lz)))
    '''
    vc = 1/vol_unitcell(c[0], c[1], c[2], c[3], c[4], c[5])
    expo = [ n[0]*xyz[0] + n[1]*xyz[1] + n[2]*xyz[2] for n in hklfcfo ]
    if diffden:
        Foc = [sqrt(x[4])-sqrt(x[3]) for x in hklfcfo]
    else:
        Foc = [sqrt(abs(x[4])) for x in hklfcfo]
    #print(Foc)
    summe = 0.0
    for F, e in zip(Foc, expo):
        summe = summe+F*cmath.exp(-2j*math.pi*(e))
    #print(vc*summe)
    rho = sqrt(vc*summe.real)+sqrt(vc*summe.imag)
    print(rho)


if __name__ == '__main__':
    filename = 'D:\Programme\DSR\example\p21c.fcf'
    c = cell_from_fcf(filename)
    print(c)
    da = read_data_from_fcf(filename, c[-1])
    
    #print('h     k     l      Fc^2       Fo^2')
    #for i in da[:10]:
    #    print(i)
    # 0.6918  0.6019  0.2177
    print('four:')
    for r in  range(0, 20):
        num = float('0.69{}18'.format(r))
        fourier(da, [num,  0.6019,  0.2177], c)
        #fourier(da, [0.2763,  0.3712,  0.4868], c)
        #fourier(da, [0.6918,  0.6019,  0.2177], c)
        #fourier(da, [0.6806,  0.5286,  0.2073], c)
    
    
    