#!/usr/bin/python
import re
import sys
import os

import symmetrize_sympy
import read
from tensors import matrix, mat2ten

import sympy
import numpy as np

from mpmath import cos as mcos
from mpmath import sin as msin
from mpmath import acos as macos
from mpmath import radians as mradians

import argparse

#this finds the location of the main.py file and ads this location to the path where modules are searched
#this way the modules have to be present only in the install directory and not in the run directory
dirname, filename = os.path.split(os.path.abspath(__file__))
sys.path.append(str(dirname))


#parses the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('op1', help='Type of the first operator')
parser.add_argument('op2', help='Type of the second operator')
parser.add_argument('-p','--projection',help='Sets a projection on an atom. Atom numbers start from 0.',default=-1)
parser.add_argument('-b','--basis',help='Sets a coordinate basis: abc for conventional crystallographic basis, i for the one used in input. \
        abc_o for orthogonalized crystalographic basis. all for outputing all.',default='i')
parser.add_argument('-e','--equivalent',help='finds response matrices for equivalent magnetic configurations. Needs output of finddsym with\
        zero moments as an input.')
parser.add_argument('--debug',help='Controls if debug output is printed. all means all debug output is printed, symmetrize means debug\
        output for symmetrizing, rename for renaming, equiv for finding the equivalent configurations',default='')
args = parser.parse_args()

#opens the output from findsym
#data = open('sym.out')
lines = sys.stdin.readlines()

op1=args.op1 #type of the first operator
op2=args.op2 #type of the second operator
atom=int(args.projection) #number of atom on which projection is done, -1 means no projection
basis=args.basis.split(',') #a list specifying the basis to be used
equiv = args.equivalent #contains name of the nonmagnetic findsym output
#debug controls whether debug output is printed
debug = args.debug.split(',')

debug_sym = False
debug_rename = False
debug_equiv = False

if 'symmetrize' in debug or 'all' in debug:
    debug_sym = True
if 'rename' in debug or 'all' in debug:
    debug_rename = True
if 'equiv' in debug or 'all' in debug:
    debug_equiv = True


#reads the input
#vec_a,b,c are needed to know the basis transformation
#syms contain the symmetries in the form that is needed by symmetr
[vec_a,vec_b,vec_c] = read.r_basis(lines)
syms = read.r_sym(lines)


#this returns the symmetrical form of spin-response tensor for atom with index atom
#if atom is -1 no projections are done
#operator types are given by op1 and op2
X = symmetrize_sympy.symmetr(syms,op1,op2,atom,debug_sym)

#outputs the tensor in the crystallographic basis used in the findsym output
#this is the one returned by symmetr
if 'abc' in basis or 'all' in basis:
    print 'using basis:'
    print 'a = ', vec_a[0], '* x + ', vec_a[1], '* y + ', vec_a[2], '* z'
    print 'b = ', vec_b[0], '* x + ', vec_b[1], '* y + ', vec_b[2], '* z'
    print 'c = ', vec_c[0], '* x + ', vec_c[1], '* y + ', vec_c[2], '* z'
    print ''

    print 'Symmetrized matrix in the abc basis even part:'
    X[0].pprint()
    print 'Symmetrized matrix in the abc basis odd part:'
    X[1].pprint()


#outputs the tensor in the input basis
if 'i' in basis or 'all' in basis:

    #transformation matrix from the original basis to the abc basis:
    T = matrix(0,3)

    T[0,0] = vec_a[0]
    T[1,0] = vec_a[1]
    T[2,0] = vec_a[2]

    T[0,1] = vec_b[0]
    T[1,1] = vec_b[1]
    T[2,1] = vec_b[2]

    T[0,2] = vec_c[0]
    T[1,2] = vec_c[1]
    T[2,2] = vec_c[2]
    #transforms the response matrix back to the input basis
    X_I = symmetrize_sympy.convert_X(X,T,debug=debug_rename)

    print ''
    print 'Symmetrized matrix in the input basis even part:'
    X_I[0].pprint()
    print 'Symmetrized matrix in the input basis odd part:'
    X_I[1].pprint()

if 'abc_o' in basis or 'all' in basis:
    #transform to cubic
    #this transforms the results to an orthogonal basis defined in the following way:
    #vector z=c
    #vector y lies in the bc plane and is ortogonal to z
    #vector x is ortogonal to y and z and they form a right-handed system
    #length of x is a, of y is b and of z is c

    T2 = matrix(0,3)

    abc = read.r_abc(lines)

    al = mradians(float(abc[3]))
    bet = mradians(float(abc[4]))
    gam = mradians(float(abc[5]))
    gam2 = macos((mcos(gam)-mcos(bet)*mcos(al))/(msin(al)*msin(bet)))

    T2[0,0] = msin(gam2)*msin(bet)
    T2[1,0] = mcos(gam2)*msin(bet)
    T2[2,0] = mcos(bet)

    T2[0,1] = 0
    T2[1,1] = msin(al)
    T2[2,1] = mcos(al)

    T2[0,2] = 0
    T2[1,2] = 0
    T2[2,2] = 1

    X_O = symmetrize_sympy.convert_X(X,T2,debug=debug_rename)

    print ''
    print 'EXPERIMENTAL: Symmetrized matrix in the orthogonalized basis even part:'
    X_O[0].pprint()
    print 'EXPERIMENTAL: Symmetrized matrix in the orthogonalized basis odd part:'
    X_O[1].pprint()

if equiv:
    #outputs also the form of the tensor for all equivalent magnetic configurations
    #equiv is a name of the nonamgnetic findsym output
    #always outputs in the input basis

    #reads data from the nonmagnetic output file
    with open(equiv) as f:
        lines_nm = f.readlines()
    syms_nm = read.r_sym(lines_nm) 
    [vec_a_nm,vec_b_nm,vec_c_nm] = read.r_basis(lines_nm)

    #transformation matrix from the magnetic basis to the input one

    T_m = np.zeros((3,3))

    T_m[0,0] = vec_a[0]
    T_m[1,0] = vec_a[1]
    T_m[2,0] = vec_a[2]

    T_m[0,1] = vec_b[0]
    T_m[1,1] = vec_b[1]
    T_m[2,1] = vec_b[2]

    T_m[0,2] = vec_c[0]
    T_m[1,2] = vec_c[1]
    T_m[2,2] = vec_c[2]

    o_m = np.array(read.r_origin(lines))

    try:
        X_I
    except:
        X_I = symmetrize_sympy.convert_X(X,T_m)

    #transformation matrix from the nonmagnetic basis to the input one
    T_nm = np.zeros((3,3))

    T_nm[0,0] = vec_a_nm[0]
    T_nm[1,0] = vec_a_nm[1]
    T_nm[2,0] = vec_a_nm[2]

    T_nm[0,1] = vec_b_nm[0]
    T_nm[1,1] = vec_b_nm[1]
    T_nm[2,1] = vec_b_nm[2]

    T_nm[0,2] = vec_c_nm[0]
    T_nm[1,2] = vec_c_nm[1]
    T_nm[2,2] = vec_c_nm[2]

    o_nm = np.array(read.r_origin(lines_nm))

    #the shift from the magnetic to the nonmagnetic
    shift = np.dot(np.linalg.inv(T_nm),o_nm)

    #atomic positions including magnetic moments
    pos = read.r_pos(lines)
    #converted positions to the input basis
    pos_t = symmetrize_sympy.convert_pos(pos,T_m,o_m)

    #this outputs all the equaivalen configurations
    #C is a conf class, it contains both the configurations and the transformed tensors
    C = symmetrize_sympy.find_equiv(X_I,op1,op2,atom,syms_nm,pos_t,T_nm,shift,debug_equiv)
    C.pprint()



