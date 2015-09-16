#!/usr/bin/python
import symmetrize_sympy
import read
import sympy
import numpy as np
import sys,os
import argparse

from funcs import *
from rename import rename
from symmetrize_tensor import symmetr_ten, print_tensor

#this finds the location of the main.py file and ads this location to the path where modules are searched
#this way the modules have to be present only in the install directory and not in the run directory
dirname, filename = os.path.split(os.path.abspath(__file__))
sys.path.append(str(dirname))

#parses the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('op1', help='Type of the first operator.')
parser.add_argument('op2', help='Type of the second operator.')
parser.add_argument('order', help='Order of the tensor.')
parser.add_argument('-p','--projection',help='Sets a projection on an atom.',default=-1)
parser.add_argument('-b','--basis',help='Sets a coordinate basis: abc for conventional crystallographic basis, i for the one used in input. \
        ',default='i')
parser.add_argument('--debug',help='Controls if debug output is printed. debug for standard debug output, Y for also printing the matrix Y\
       time for printing information about how much time it takes.',default='')
args = parser.parse_args()

op1=args.op1 #type of the first operator
op2=args.op2 #type of the second operator
order = int(args.order)

atom=int(args.projection)

basis=args.basis.split(',') #a list specifying the basis to be used
debug = args.debug.split(',')

debug_normal = False
debug_Y = False
debug_time = False

if 'debug' in debug:
    debug_normal = True
if 'Y' in debug:
    debug_normal = True
    debug_Y = True
if 'time' in debug:
    debug_time = True

lines = sys.stdin.readlines()
[vec_a,vec_b,vec_c] = read.r_basis(lines)
syms = read.r_sym(lines)

if 'abc' in basis:

    X = symmetr_ten(syms,'s','v',order,proj=atom,T=None,debug=debug_normal,debug_Y=debug_Y,debug_time=debug_time)

    print 'using basis:'
    print 'a = ', vec_a[0], '* x + ', vec_a[1], '* y + ', vec_a[2], '* z'
    print 'b = ', vec_b[0], '* x + ', vec_b[1], '* y + ', vec_b[2], '* z'
    print 'c = ', vec_c[0], '* x + ', vec_c[1], '* y + ', vec_c[2], '* z'
    print ''

    print_tensor(X)

if 'i' in basis:


    T = sympy.zeros(3)

    T[0,0] = vec_a[0]
    T[1,0] = vec_a[1]
    T[2,0] = vec_a[2]

    T[0,1] = vec_b[0]
    T[1,1] = vec_b[1]
    T[2,1] = vec_b[2]

    T[0,2] = vec_c[0]
    T[1,2] = vec_c[1]
    T[2,2] = vec_c[2]

    T = make_rational(T)

    X = symmetr_ten(syms,'s','v',order,proj=atom,T=T,debug=debug_normal,debug_Y=debug_Y,debug_time=debug_time)

    print 'In the input basis:'
    print ''

    print_tensor(X)

    #X = convert_tensor(X,T)

