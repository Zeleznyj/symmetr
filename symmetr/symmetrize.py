from __future__ import print_function
from __future__ import absolute_import
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from builtins import range
from builtins import object
import re
import copy
import sys
import math
import time

import sympy
import numpy as np
import numpy.linalg
import mpmath   

from .tensors import matrix, mat2ten
from .fslib import transform_position
from .conv_index import *

class params_trans(object):
    def __init__(self,op1,op2,op3,l,T=None,sym_format='findsym'):
        self.op1 = op1
        self.op2 = op2
        self.op3 = op3
        self.l = l
        self.T = T
        self.sym_format = sym_format

class SymmetrOpt(object):
    def __init__(self,num_prec=None,debug=False,debug_time=False,debug_Y=False,round_prec=None):
        self.num_prec = num_prec
        self.debug = debug
        self.debug_time = debug_time
        self.debug_Y = debug_Y
        self.round_prec = round_prec

def symmetr(syms,X,trans_func,params,opt=None):
    """
    This symmetrizes a tensor X given a list of symmetries and a transformation function.

    This function should be quite general and is now used for all symmetrizing.

    Args:
        syms: list of symmetry operations
        X: tensor - must be a tensor class
        trans_func: function that transforms the tensor X using symmetry sym
            trans_func must work in the following way:
            X_trans = trans_func(X,sym,params)
            If trans_func returns None then the symmetry operation is ignored
        params: parameters to be sent to function trans_func

    Returns:
        X_trans: the symmetry restricted form of tensor X
    """
    if opt is None:
        opt = SymmetrOpt()
    num_prec = opt.num_prec
    debug = opt.debug
    debug_time = opt.debug_time
    debug_Y = opt.debug_Y

    if debug:
        print('')
        print('======= Starting symmetrizing =======')

    #we do a loop over all symmetry operations, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
    #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
    for sym in syms:
        
        if debug:
            print('Symmetry:') 
            print(sym)
            print('')

        X_trans = trans_func(X,sym,params)
        if X_trans is None:
            continue
            
        if debug:
            print('')
            print('Current form of the tensor:')
            print('')
            X.pprint()
            print('')
            print('Transformed tensor:')
            print('')
            X_trans.pprint()
        #X.pprint()
        #X_trans.pprint()
        #The tensor must be equal to the transformed tensor, this give us a system of linear equations.
        #matrix Y is a matrix that represents this system, ie the system X-X_trans = 0
        #we reverse the order of the rows
        # it doesn't really matter but the results are more natural this way

        if debug_time:
            t1 = time.clock()

        Y = sympy.zeros(X.dim1**X.dim2)

        rev_inds = list(reversed(X.inds))
        #we do a loop over all rows of the matrix Y - ie over all linear equations
        n = -1
        for ind1 in X:
            n += 1
            m = -1

            #if this is zero, then we do not have to do any substituting so this saves quite a lot of time
            if X[ind1]-X_trans[ind1] == 0:
                is_zero = True
            else:
                is_zero = False

            speed_test = True
            if speed_test:
                inds = re.findall(r'x[0-9]+',sympy.srepr(X[ind1]-X_trans[ind1]))
                for ind2 in reversed(inds):
                    m_index = (int(i) for i in re.findall(r'[0-9]',ind2))
                    m = rev_inds.index(tuple(m_index))
                    Y_p = X[ind1]-X_trans[ind1]
                    #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                    for ind3 in inds:
                        if ind2 == ind3:
                            Y_p = Y_p.subs(X[ind3],1)
                        else:
                            Y_p = Y_p.subs(X[ind3],0)

                    Y[n,m] = Y_p

            else:
                for ind2 in rev_inds:
                    m += 1
                    Y_p = X[ind1]-X_trans[ind1]
                    #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                    for ind3 in X.inds:
                        if ind2 == ind3:
                            Y_p = Y_p.subs(X.x[ind3],1)
                        else:
                            Y_p = Y_p.subs(X.x[ind3],0)

                    Y[n,m] = Y_p

        if debug_time:
            t2 = time.clock()
            print('Time for constructing Y: ', t2-t1)
            t1 = time.clock()

        #this transforms the matrix into the Reduced row echelon form
        #piv are the indeces o the pivot columns
        if num_prec is None:
            [rref,piv] = Y.rref()
        else:
            def zerofunc(x):
                try:
                    a = abs(x) < num_prec
                    if a:
                        result = True
                    else:
                        result = False
                except:
                    result = None
                #print(result)
                return result
            #sympy.pprint(Y)
            #[rref,piv] = Y.rref()        
            [rref,piv] = Y.rref(iszerofunc=lambda x:abs(x)<num_prec)        
            #[rref,piv] = Y.rref(iszerofunc=zerofunc)        

        if debug_time:
            t2 = time.clock()
            print('Time for reducing Y to reduced row echelon form: ', t2-t1)

        if debug_Y:
            print('')
            print('Matrix representing the linear equation system that has to be satisfied: (right hand side is zero)')
            sympy.pprint(Y)
            print('')
            print('Reduced row echelon form and indeces of the pivot columns:')
            sympy.pprint(rref)
            print(piv)
            print('')

        #a loop over all the pivots: it's the pivots that give interesting information
        for j in list(reversed(piv)):

            
            #find the row of pivot j
            found = False
            i = X.dim1**X.dim2-1
            while found == False:
                if rref[i,j] == 1:
                    found = True
                else:
                    i = i-1
            
            if debug:
                print('')
                print('considering pivot ', i,j)

            tmp = 0
            #now we just make use of the linear equation that holds for this pivot
            #keep in mind that the rows are in reversed order
            for ll in range(j+1,X.dim1**X.dim2):
                tmp = tmp - rref[i,ll]*X.x[rev_inds[ll]]
            X = X.subs(X.x[rev_inds[j]],tmp)

            if debug:
                print('substituting ', end=' ')
                sympy.pprint(X.x[rev_inds[j]])
                print(' for ', end=' ')
                sympy.pprint(tmp)
                print('')

        if debug:
            print('Current form of the tensor:')
            X.pprint()
            print('')

    if opt.round_prec is not None:
        X.round(opt.round_prec)
        
    if debug:
        print('Symmetrized tensor:')
        X.pprint()
        print('')
        print('======= End symmetrizing =======')

    return X

def even_odd(Xs):
    """Finds whether the first part of the response tensor is even or odd
    
    Args:
        op1,op2,op3: the operator types
        Returns: either ('even','odd') or ('odd','even')
    """
    if Xs[0].is_even() and not Xs[1].is_even():
        return('even','odd')
    elif not Xs[0].is_even() and Xs[1].is_even():
        return('odd','even')
    else:
        raise Exception('Wrong transformation under time-reversal')


def symmetrize_res(symmetries,X,proj=-1,s_opt=None):

    syms_sel = []
    for sym in symmetries:
        
        if s_opt.debug:
            print('Symmetry:') 
            print(sym)
            print('')
            if proj != -1:
                print('Symmetry transforms the atom ', proj, ' into atom ', sym.permutations[proj])
                if  sym.permutations[proj] != proj:
                    print('Skipping symmetry')
                    print('')

        #if there is a projection set up we only consider symmetries that keep the atom invariant
        if proj == -1 :
            take_sym = True
        elif sym.permutations[proj] == proj:
            take_sym = True
        else:
            take_sym = False

        if take_sym:
            syms_sel.append(sym)

    def trans_func(X,sym,params):
        return X.transform(sym)
    X = symmetr(syms_sel,X,trans_func,None,s_opt)

    return X

def symmetrize_same_op(X,s_opt=None):

    perm = {}
    perm[0] = 1
    perm[1] = 0
    perms=[perm] 

    def trans_func(X,perm,params):
        X_T = X.copy0()

        for ind in X:
            ind_T = [0]*len(ind)
            for i in range(len(ind)):
                ind_T[i] = ind[perm[i]]
            ind_T = tuple(ind_T)
            X_T[ind_T] = X.T_comp * X[ind]

        ind_types = [0] * X.dim2
        for i in range(X.dim2):
            ind_types[i] = X.ind_types[perm[i]]
        X_T.ind_types = tuple(ind_types)

        for i in range(X.dim2):
            if X_T.ind_types[i] != X.ind_types[i]:
                X_T.reverse_index(i)

        return X_T

    X = symmetr(perms,X,trans_func,None,s_opt)
        
    return X

def symmetr_AB(syms,X,atom1,atom2,round_prec=None):
    """
    Tries to transform the tensor projected on one atom to a different atom

    Args:
        syms: The symmmetry operations. Format as outputted by read.py
        X: The input tensor.
        op1: The first operator.
        op2: The second operator.
        atom1: The atom on which X is projected.
        atom2: The atom on which X is transformed.
        T (Optional[matrix]): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.

    Returns:
        X_trans: The transformed tensor.
    """

    X_trans = []

    found = False
    for sym in syms:
        #there will usually be more symmetries that transform from atom1 to atom2, we need only one, as they all
        #give the same results
        if sym.permutations[atom1] == atom2 and not found:
            found = True
            for l in range(2):
                X_trans.append(X[l].transform(sym))

    if found:
        if round_prec is not None:
            for X in X_trans:
                X.round(round_prec)
        return X_trans
    else:
        return None
