import re
import copy
import sys
import math
import time

import sympy
import numpy as np

from tensors import matrix, mat2ten, tensor
from read import transform_position
from funcs import *
from symmetrize import symmetr,params_trans

def symmetrize_exp(symmetries,op1,op2,order,proj=-1,T=None,debug=False,debug_Y=False,debug_time=False):
    """
    Returns a symmetrical form of a tensor, which describes a term in an expansion of linear response tensor in magnetization.

    Args:
        symmetries: The symmetry operations.
        op1 (string): The type of the first operator.
            Can be set to either 's', which means spin or 'v' which means velocity.
        op2 (string): The type of the second operator.
        order: The order of the expansion term. Zeroth order corresponds to rank 2 tensor, first order to rank 3 etc.
        proj (Optional[int]): Determines projection on atom. Defaults to -1.
            If set to -1, there is no projection.
            If set to positive integer, it determines atom number.
        debug (Optional[boolean]): Defaults to false. If set to true, additional debug output is printed.
        debug_Y (Optional[boolean]): Defaults to false. If set to true, the Y matrix is printed. This can be very large.
        debug_time (Optional[boolean]): Defaults to false. If set to true, information about how long each step takes.

    Returns:
        X: The most general form of the tensor allowed by symmetry.
    """

    #This defines the starting tensor.
    #Uses a symbolic tensor class made using sympy.
    
    X = tensor('s',3,order+2)

    if debug:
        print ''
        print '======= Starting symmetrizing ======='

    #we do a loop over all symmetry operations, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
    #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
    syms_sel = []
    for sym in symmetries:
        
        if debug:
            print 'Symmetry:' 
            print sym
            print ''
            if proj != -1:
                print 'Symmetry transforms the atom ', proj, ' into atom ', sym_type(proj,sym)
                if sym_type(proj,sym) != proj:
                    print 'Skipping symmetry'
                    print ''

        #if there is a projection set up we only consider symmetries that keep the atom invariant
        if proj == -1 :
            take_sym = True
        elif sym_type(proj,sym) == proj:
            take_sym = True
        else:
            take_sym = False
        
        #we only take the symmetry operations that don't contain time-reversal
        #those with time-reversal give no new information
        if sym[3] == '-1':
            take_sym = False

        if take_sym:
            syms_sel.append(sym)

    params = params_trans(op1,op2,None,order % 2,T,'findsym')
    X = symmetr(syms_sel,X,transform_exptensor_params,debug=debug,debug_time=debug_time,debug_Y=debug_Y)

    return X


def create_rank2(ten,xyz=False):
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

def print_tensor(ten,latex=False,xyz=False,no_newline=False):
    """
    Prints the expansion tensor in a nice form.

    Not tested for higher order than 1!!!
    """
    
    X = create_rank2(ten,xyz=xyz)

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






