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
from math import sqrt, cos, sin, exp, pi, radians


'''
Created on 26.03.2015

@author: Daniel Kratzert
'''


def cell_from_fcf(filename):
    '''
    reads unit cell from fcf
    '''
    a, b, c, al, be, ga, list_code, F000 = (None, )*8  
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
                if line.startswith('_exptl_crystal_F_000'):
                    F000 = float(line.split()[1])
                if all([a, b, c, al, be, ga, list_code, F000]):
                    break
    except(IOError) as e:
        print(e)
        print('Unable to read {} file.'.format(filename))
        sys.exit()
    return (a, b, c, al, be, ga, list_code, F000)
                

def read_data_from_fcf(filename, type):
    '''
    read Fo and Fc data from fcf
    '''
    data = []
    read_data = False
    count = 0
    if type == '4':
        print('You need a LIST 6 fcf to make the fourier synthesis.')
        sys.exit()
    if not type == '6':
        print('Wrong fcf type {}'.format(type))
        sys.exit()
    try:
        with open(filename, 'r') as f:
            for line in f:
                if read_data:
                    line = line.split()
                    try:
                        line[1]
                    except:
                        continue
                    # h k l Fo Fc phi
                    # 0 1 2 3  4  5
                    line = [int(line[0]), int(line[1]), int(line[2]), 
                            round(sqrt(float(line[3])), 4), float(line[5]), radians(float(line[6]))]
                    data.append(line)
                    #data.append([int(line[0]), int(line[1]), int(line[2]), 
                    #            float(line[3]), float(line[4])])
                    continue
                if line.startswith(' _refln'):
                    count = count+1
                    if count == 7:
                        read_data = True
    except(IOError):
        print('Unable to read {} file.'.format(filename))
        sys.exit()
    return data
    
def fourier(hklfcfo, xyz, c, diffden=True):
    '''
    rho = 1/V*SUM([Fo-Fc]*exp(i*phi)*exp(-i2pi(hx+ky+lz)))
    '''
    v = vol_unitcell(c[0], c[1], c[2], c[3], c[4], c[5])
    expo = [ (n[0]*xyz[0] + n[1]*xyz[1] + n[2]*xyz[2]) for n in hklfcfo ]
    #Fo -> x[3]
    #Fc -> x[4]
    if diffden:
        Foc = [abs(x[3])-abs(x[4]) for x in hklfcfo]
    else:
        Foc = [x[4] for x in hklfcfo]
    philist = [p[5] for p in hklfcfo]
    F000 = c[7]
    summe = F000
    for F, hklxyz, phi in zip(Foc, expo, philist):
        summe = summe+abs(F)*cmath.exp(-1j*2*pi*(hklxyz-phi))
    print
    print((1/v)*summe)
    #summe = (1/v)*(sqrt(summe.real)+sqrt(summe.imag))
    #summe = (1/v)*summe.real
    summe = (1/v)*abs(summe)
    rho = abs(summe)
    print(rho)


if __name__ == '__main__':
    #filename = '/Users/daniel/Documents/DSR/p21c.fcf'
    #filename = 'D:\Programme\DSR\example\p21c.fcf'
    filename = r'D:\Programme\DSR\example\nte03a.fcf'
    
    c = cell_from_fcf(filename)
    print(c)
    da = read_data_from_fcf(filename, type=c[6])
    #print('h     k     l      Fc^2       Fo^2')
    #for i in da[:10]:
    #    print(i)
    # 0.6918  0.6019  0.2177
    print('four:')
    for r in  range(0, 1):
        #num = float('0.69{}18'.format(r))
        #fourier(da, [num,  0.6019,  0.2177], c)
        #fourier(da, [0.406763,   0.436327,   0.352405], c)
        fourier(da, [0.6015, -0.0444,  0.9054], c)
        #fourier(da, [0.7427,  0.5689,  0.2124], c)
        #fourier(da, [0.177383,    0.298455,   0.290893], c)
    
    
    