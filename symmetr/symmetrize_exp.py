from __future__ import print_function
from __future__ import absolute_import
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from builtins import str
from builtins import range
import re
import copy
import sys
import math
import time

import sympy
import numpy as np

from .tensors import matrix, mat2ten, tensor
from .fslib import transform_position
from .symmetrize import symmetr,params_trans
from .symT import convert_vec

def create_rank2(ten,n_op=2,xyz=False):
    """
    Creates a rank 2 tensor that includes magnetic moment explicitely.
    """

    X = matrix(0,3)
    X.x = ten.x
    X.v = ten.v

    if ten.dim2 > 2:

        m = {}
        for i in range(3):
            if xyz:
                names = ['m_x','m_y','m_z']
                name = names[i]
            else:
                name = 'm%s' % i
            m[i] = sympy.symbols(name)


        for ind in ten:
            M = 1
            for i in range(2,len(ind)):
                M *= m[ind[i]]
            X[ind[0],ind[1]] += ten[ind]*M


    else:
        for ind in ten:
            X[ind] = ten[ind]
    return X

def sub_m(ten,n_op,xyz=False):
    
    X = tensor(0,3,n_op)

    if ten.dim2 > n_op:
        m = {}
        for i in range(3):
            if xyz:
                names = ['m_x','m_y','m_z']
                name = names[i]
            else:
                name = 'm%s' % i
            m[i] = sympy.symbols(name)
        for ind in ten:
            M = 1
            for i in range(n_op,len(ind)):
                M *= m[ind[i]]
            ind_m = []
            for i in range(n_op):
                ind_m.append(ind[i])
            ind_m = tuple(ind_m)
            X[ind_m] += ten[ind]*M

    else:
        for ind in ten:
            X[ind] = ten[ind]

    return X

def print_tensor(ten,n_op,latex=False,xyz=False,no_newline=False):
    """
    Prints the expansion tensor in a nice form.

    Not tested for higher order than 1!!!
    """
    
    X = sub_m(ten,n_op,xyz=xyz)

    if not latex:
        X.pprint()
    else:
        X.pprint(latex=True,no_newline=no_newline)

def simplify_tensor(ten,xyz=False,index_from_1=False):
    """
    Renames the variables of the tensor and simplifies it.
    """

    X = create_rank2(ten,xyz=xyz)
    xinds = list(set(re.findall(r'x[0-9]+',sympy.srepr(X))))

    #contains the new indices
    xn = {}
    if not index_from_1:
        for i in range(len(xinds)):
            xn[i] = sympy.symbols('x'+str(i))
    else:
        for i in range(len(xinds)):
            xn[i] = sympy.symbols('x'+str(i+1))

    for ind in X:
        for i in range(len(xinds)):
            X[ind] = X[ind].subs(xinds[i],xn[i])
        X[ind] = sympy.simplify(X[ind])

    return X

def index_from_1(X,rank=2):
    """
    Takes a rank 3 tensor and rename the indeces so that the numbering starts from 1 and not 0.
    """

    if rank == 2:
    
        xinds = list(set(re.findall(r'x[0-9]+',sympy.srepr(X))))
        xn = {}
        xnz = {}
        for i in range(3):
            for j in range(3):
                xnz[(i,j)] = sympy.symbols('z'+str(i+1)+str(j+1))

        for i in range(3):
            for j in range(3):
                xn[(i,j)] = sympy.symbols('x'+str(i+1)+str(j+1))

        for ind in X:
            for i in range(3):
                for j in range(3):
                    X[ind] = X[ind].subs(X.x[i,j],xnz[(i,j)])

        for ind in X:
            for i in range(3):
                for j in range(3):
                    X[ind] = X[ind].subs(xnz[(i,j)],xn[(i,j)])

        return X

    if rank == 3:

        xinds = list(set(re.findall(r'x[0-9]+',sympy.srepr(X))))
        xn = {}
        xnz = {}
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    xnz[(i,j,k)] = sympy.symbols('z'+str(i+1)+str(j+1)+str(k+1))

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    xn[(i,j,k)] = sympy.symbols('x'+str(i+1)+str(j+1)+str(k+1))

        for ind in X:
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        X[ind] = X[ind].subs(X.x[i,j,k],xnz[(i,j,k)])

        for ind in X:
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        X[ind] = X[ind].subs(xnz[(i,j,k)],xn[(i,j,k)])

        return X

def is_collinear(mags,prec=1e-5):
    res = True
    first = sorted(list(mags.keys()))[0]
    for i in mags:
        if abs(mags[first].dot(mags[i])) < 1-prec:
            res = False
            break

    return res

def convert_mags(mags,sym):
    mags_R = {}
    for atom in mags:
        mags_R[atom] = convert_vec(mags[atom],sym.get_R('s'))
    mags_Rp = {}
    perms = sym.permutations
    for atom in mags:
        mags_Rp[perms[atom]] = mags_R[atom]

    return mags_Rp

def get_L_trans(mags,sym,debug=False):

    initial_signs = []
    
    first = sorted(list(mags.keys()))[0]
    for atom in mags:
        if mags[first].dot(mags[atom]) > 0:
            initial_signs.append(1)
        else:
            initial_signs.append(-1)

    mags_R = convert_mags(mags,sym)
    signs = []
    for atom in mags_R:
        if mags_R[first].dot(mags_R[atom]) > 0:
            signs.append(1)
        else:
            signs.append(-1)

    if debug:
        print('')
        print('transformed magnetic moments:')
        for atom in  sorted(list(mags_R.keys())):
            print(atom, ":")
            sympy.pprint(mags_R[atom])

        print('initial signs of magnetic moments')
        print(initial_signs)

        print('signs of the trasnformed moments')
        print(signs)

    if initial_signs != signs:
        return None
    else:
        mag_0_R = convert_vec(mags[first],sym.get_R('s'))
        sign_0 = mag_0_R.dot(mags[first])
        if sign_0 > 0:
            sign = 1
        else:
            sign = -1
        return sign * sym.get_R('s')

def def_syms_L(mags,syms,prec=1e-5,debug=False):

    #select nonzero magnetic moments
    mags_dict = {}
    for i,mag in enumerate(mags):
        if mag.norm() > prec:
            mags_dict[i+1] = mag

    if debug:
        print('Starting defining L transformation')
        print('Nonzero magnetic moments')
        print(mags,mags_dict)
        for atom in  sorted(list(mags_dict.keys())):
            print(atom, ":")
            sympy.pprint(mags_dict[atom])

    if len(mags_dict) == 0:
        print("!!!Warning!!!: no magnetic moments defined in the input. Assuming a ferromagnetic system.")
        for sym in syms:
            sym.def_custom_R('L',sym.get_R('s'))
        syms_L = syms

    else:

        if not is_collinear(mags_dict,prec):
            raise Exception("Expansions only work for collinear magnetic systems")

        syms_L = []
        for sym in syms:

            if debug:
                print('')
                print('Taking symmetry:')
                print(sym)

            L_trans = get_L_trans(mags_dict,sym,debug=debug)
            if debug:
                print('L_trans:')
                print(L_trans)
            if L_trans is not None:
                sym.def_custom_R('L',L_trans)
                syms_L.append(sym)

    return syms_L







            





