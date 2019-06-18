from __future__ import print_function
from __future__ import absolute_import
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from builtins import str
from builtins import range
import sys
import os
import re
import numpy as np
import sympy 
from .tensors import matrix
from .fslib import r_sym
from .symmetry import Symmetry,findsym2sym

dirname, filename = os.path.split(os.path.abspath(__file__))
sys.path.append(str(dirname))

def read_all_syms(hexag):
    """
    Generates a list of all crystallographic symmetry operations.

    Args:
        hexag(boolean): true if the conventional crystallographic group of the nonmagnetic crystal
            is hexagonal or trigonal, false otherwise.
    """

    if not hexag:
        syms_file = '/findsym/syms_table.dat'
    else:
        syms_file = '/findsym/syms_table_hex.dat'

    with open(dirname+syms_file) as f:
        lines = f.readlines()

    syms_list = []
    syms_list.append("_space_group_symop.magn_operation_mxmymz")
    for i,line in enumerate(lines):
        split = re.findall('[a-zA-Z0-9\-,\+]+',line)
        sym = str(i) + '  ' + split[1]+ ',+1  ' + split[2] 
        if not hexag:
            syms_list.append(sym)
        else:
            if 'T' not in line:
                syms_list.append(sym)


    syms_list.append(" ")
    syms = r_sym(syms_list,syms_only=True)
    mats = []
    for sym in syms:
        sym = findsym2sym(sym)
        mats.append(sym.get_R('x'))

    return mats

def noso_syms(syms,mag_conf,hexag,prec=1e-5,debug=False):
    """
    Takes symmetry operations in a nonmagnetic crystal and returns symmetry operations of a
        magnetic crystal without spin-orbit coupling.

    Args:
        syms: list of symmetry operations as returned by r_sym
        mag_conf([sympy(3)]): list of magnetic moments of every atom 
        hexag(boolean): true if the conventional crystallographic group of the nonmagnetic crystal
            is hexagonal or trigonal, false otherwise.

    Returns:
        syms_noso: list of symmetry operations without spin-orbit coupling. Format is different than 
            normal, the symmetry operations are in matrix form. Some functions can take argument
            sym_format='mat' to switch to this format.

    The algorithm:
        * Without SOC an arbitrary spin-rotation is a symmetry in a non-
            magnetic system. In a magnetic system,only spin-rotations that respect the magnetic order 
            are symmetries.
        * Thus we take every symmetry operation of the nonmagnetic system and try to see if we find some
            spin-rotations such that the symmetry+spin-rotation is a symmetry of the magnetic system.
        * The problem is that there is usually infinitely many such spin-rotations. It's not trivial to find them 
            in general and also to choose some finite number out of these.
        * Thus we use a trick. We do everything in the conventional coordinate basis because there the symmetry 
            operations have a simple form. We look for a combination of the symmetry operation+spin-rotation, which
            keeps the magnetic order invariant. We look at every crystallographic symmetry operation and see if it
            keeps the magnetic order invariant.
        * Rather than looking for the spin-rotatin we directly look for the combination of the symmetry +
            spin-rotation. Since the spin-rotation is arbitrary this is fine.
        * Since the spin-rotation is a proper rotation it must have determinant 1. Thus if there is time reversal,
            the total matrix (spin-rotation + the symmetry operation of the nonmagnetic crystal) must then have
            determinant -1. If there is no time-reversal in the symmetry operation it must have det +1. This is
            the only requirement.
        * In this way we get usually more than 1 symmetry operation for every nonmagnetic symmetry operation.
        * However, we do not get all the symmetry operations, in fact usually there is infinitely many of them.
        * This approach will probably work fine if the magnetic moments lie some high-symmetry directions. If they
            are not, then this approach will not find the correct spin-rotations and the symmetry will thus not
            be complete.
        * The symmetry should not be wrong, but can be incomplete - some elements can be shown as nonzero even
            though they should vanish or some relations missing.
    """

    if debug:
        if hexag:
            print('The conventional coordinate system is hexagonal.')
        else:
            print('The conventional coordinate system is not hexagonal.')
        print('')

    if debug:
        print('Magnetic moments are:')
        for i,mag in enumerate(mag_conf):
            print(i+1,' ', end=' ')
            sympy.pprint(mag)

    mats = read_all_syms(hexag)
    if debug:
        print('')
        print('list of all crystallographic symmetry operations:')
        for i,mat in enumerate(mats):
            print('')
            print(i+1)
            sympy.pprint(mat)
            print('')

    syms_noso = []
    start_new = 0
    for nsym,sym in enumerate(syms):
        start_new = len(syms_noso)

        if debug:
            print('')
            print('taking symmetry operation ', nsym+1)
            print('space part:')
            sympy.pprint(sym[0])
            print('magnetic part:')
            sympy.pprint(sym[2])
            print('time reversal:')
            sympy.pprint(sym[3])
            print('permutations:')
            print(sym[4])

        nmag = 0
        mag_is = []
        for i,mag in enumerate(mag_conf):
            if sympy.sqrt(mag[0]**2+mag[1]**2+mag[2]**2) > prec:
                nmag += 1
                mag_is.append(i)
        
        Y = sympy.zeros(nmag*3,10)

        for i in range(nmag):
            for j in range(3):
                for k in range(3):
                    Y[i*3+j,j*3+k] = mag_conf[mag_is[i]][k]
                #Y[i*3+j,9] = mag_conf[sym_type(mag_is[i]+1,sym)-1][j]
                Y[i*3+j,9] = mag_conf[sym.permutations[mag_is[i]+1]-1][j]

        if debug: 
            print('')
            print('Matrix of the system to be solved.')
            sympy.pprint(Y)
        
        Rs = matrix('s',3,name='R') 
        sol = sympy.solve_linear_system(Y,Rs.x[0,0],Rs.x[0,1],Rs.x[0,2],Rs.x[1,0],Rs.x[1,1],Rs.x[1,2],Rs.x[2,0],\
                                         Rs.x[2,1],Rs.x[2,2])
        if debug:
            print('solution of the system')
            print(sol)

        for s in sol:
            Rs.subs(s,sol[s])

        if debug:
            print('')
            print('general form of transformation matrix that keeps the magnetic order invariant:')
            sympy.pprint(Rs.mat())
            print('')

        for mat in mats:
            if debug:
                print('Considering the matrix:')
                sympy.pprint(mat)
            fits = True
            for i in range(3):
                for j in range(3):
                    if Rs.x[i,j] in sol:
                        if mat[i,j] != sol[Rs.x[i,j]]:
                            fits = False

            if debug:
                if fits:
                    print('The matrix is compatible with the magnetic order.')
                else:
                    print('The matrix is not compatible with the magnetic order.')

            if fits:
                if mat.det() == 1 and sym.has_T:
                    fits = False
                    if debug:
                        print('Improper spin rotation, thus not taking this matrix.')
            if fits:
                sym_new = sym.copy()
                sym_new.Rs = mat
                syms_noso.append(sym_new)

        if debug:
            print('original symmetry operation:')
            print('space part:')
            sympy.pprint(sym[0])
            print('magnetic part:')
            sympy.pprint(sym[2])
            print('')

            print('new symmetry operations (magnetic part only):')
            for sm in syms_noso[start_new:]:
                sympy.pprint(sm[2])
                print('')

    return syms_noso

