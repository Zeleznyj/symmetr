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
from scipy.linalg import lu,qr
from scipy.io import savemat
import numpy as np
import numpy.linalg
from numpy.linalg import svd
import mpmath   

from .tensors import matrix, mat2ten, tensor2Y
from .fslib import transform_position
from .conv_index import *

def num_rref(Y,prec=15):

    Y = np.matrix(Y)

    pl,U = lu(Y,permute_l=True)
    #P,L,U = mpmath.lu(Y)

    #print('L,U decomposition')
    #print(pl)
    #print(U)

    pivots = []
    for i in range(U.shape[0]):
        for j in range(U.shape[1]):
            pivot_found = False
            if U[i,j] != 0 and not pivot_found:
                pivots.append(j)
                U[i,:] /= U[i,j]
                break
    try:
        np.linalg.inv(pl)
    except:
        print('The PL matrix is singular. The output of the code may be wrong!!!!')
    return U,pivots

def QR_rref(V,num_prec=1e-10):
    Q,R = qr(V)
    return U_to_rref(R,num_prec=num_prec)

def U_to_rref(U,num_prec=1e-10):
    U_out = U.copy()
    pivots = []
    for row in range(U_out.shape[0]):
        found_pivot = False
        for col in range(U_out.shape[1]):
            if abs(U_out[row,col]) >= num_prec:
                if col not in pivots:
                    pivots.append(col)
                    U_out[row,:] = U_out[row,:] / U_out[row,col]
                    found_pivot = True
                    break
                else:
                    U_out[row,:] -= U_out[pivots.index(col),:] * U_out[row,col]
        if not found_pivot:
            raise Exception('Did not find pivot!')
            
    pivots_i = [(p,i) for i,p in enumerate(pivots)]
    pivots_i.sort()
    sorted_pivots,permutation = zip(*pivots_i)
    if permutation != tuple(range(len(pivots))):
        U_out_perm = U_out.copy()
        for row in range(U_out.shape[0]):
            U_out_perm[row,:] = U_out[permutation[row],:]
    
        pivots = sorted_pivots
        U_out = U_out_perm

    for piv_row,col in enumerate(pivots):
        for row in range(U_out.shape[0]):
            if row < piv_row:
                if abs(U_out[row,col]) >= num_prec:
                    U_out[row,:] -= U_out[row,col] * U_out[piv_row,:]
    return U_out,pivots

def get_rref(U,num_prec=1e-10):
    """
    This transforms the matrix in to a reduced row echelong form. The algorithm is taken
    from octave and should be exactly the same.

    Columns that are smaller than num_prec are set to zero.
    """
    U_o = U.copy()
    pivots = []
    row = 0
    rows,columns = U.shape
    for col in range(columns):
        m,pivot = [np.max(np.abs(U_o[row:,col])),np.argmax(np.abs(U_o[row:,col]))]
        pivot = pivot + row
        
        if m < num_prec:
            U_o[row:,col] = np.zeros((rows-row))
        else:
            pivots.append(col)
            U_o[[pivot,row],:] = U_o[[row,pivot],:]
            U_o[row,:] = U_o[row,:] / U_o[row,col]
            for r in range(rows):
                if r != row:
                    U_o[r,:] -= U_o[row,:] * U_o[r,col]
            if row == rows-1:
                break
            else:
                row += 1
    
    return U_o,pivots

class params_trans(object):
    def __init__(self,op1,op2,op3,l,T=None,sym_format='findsym'):
        self.op1 = op1
        self.op2 = op2
        self.op3 = op3
        self.l = l
        self.T = T
        self.sym_format = sym_format

class SymmetrOpt(object):
    def __init__(self,num_prec=None,debug=False,debug_time=False,debug_Y=False,round_prec=None,numX=False):
        self.num_prec = num_prec
        self.debug = debug
        self.debug_time = debug_time
        self.debug_Y = debug_Y
        self.round_prec = round_prec
        self.numX = numX

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
    if opt.numX and num_prec is None:
        raise Exception('In the numX mode, num_prec must be specified!')

    if debug_time:
        print('Symmetrize starting')

    if num_prec is None:
        algo = 1
        algo_solve = 0
    else:
        algo = 3
        algo_solve = 1

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

        if debug_time:
            t0 = time.perf_counter()

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
            t1 = time.perf_counter()
            print('Time for transforming the tensor: ',t1-t0) 
            t0 = time.perf_counter()

        #equal = (X_trans == X)
        if debug_time:
            t3= time.perf_counter()
            print('Time for checking equality: ',t3-t0)
        #if equal:
        #    continue

        if not opt.numX:
            Y = tensor2Y(X-X_trans,Y_format='sympy',algo=algo,reverse=True,td=False)
        else:
            Y = X.t - X_trans.t

        """
        algo_solve == 0: We use the sympy rref function. This seems to have issues though even if num_prec is set.
        algo_solve == 1: We svd decomposition to find the nullspace and the convert to rref using get_rref.
        """
        if algo_solve == 0:

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
                [rref,piv] = Y.rref(iszerofunc=lambda x:abs(x)<num_prec)        
                #[rref,piv ] = num_rref(Y)
                #[rref,piv] = Y.rref(iszerofunc=zerofunc)        

            if debug_time:
                t2 = time.perf_counter()
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

        elif algo_solve == 1:
            """
            We use the svd decomposition to find the nullspace, this is given by V2:
            rows of V2 are vectors which form the basis of the Y nullspace. This means that 
            any combination of these vectors is a solution of Yx = 0.

            Therefore any solution of Yx =0, must be written as:
            x = a_i * v_i,
            where v_i are the rows of V2. 
            By transforming V2 into the reduced row echelon form, we can then
            eliminate some variables from X. The information is extracted from columns of the V2_rref matrix
            (nor rows like in the case of Y_rref). The pivot rows simply tell us that a_i = x_i, the non-pivot
            columns then give the relation of the other x_i's in terms of the pivot x_i's.
            """

            if debug_time:
                ts0 = time.perf_counter()

            if opt.numX:
                U,S,V = svd(Y)
            else:
                U, S, V = svd(np.array(np.flip(Y, axis=1), dtype=np.float))
            if debug_time:
                ts1 = time.perf_counter()
                print('Time for svd: {}'.format(ts1-ts0))
            zero_singulars = [i for i,s in enumerate(S) if abs(s) < num_prec]

            #if there are no singular values, then the tesor has to be zero:
            #thus we don't have to do any substituting
            if len(zero_singulars) == 0:
                X = X.copy0()
            else:
                #V2 = sympy.zeros(len(zero_singulars),Y.shape[1])
                V2 = np.zeros((len(zero_singulars),Y.shape[1]))
                for i in range(len(zero_singulars)):
                    V2[i,:] = V[zero_singulars[i],:]
                if debug_Y:
                    print('Singular values: ',S)
                    print(V2)
                V2_rref,pivs = get_rref(V2,num_prec)
                if debug_Y:
                    print(V2_rref)
                if debug_time:
                    ts0 = time.perf_counter()
                    print('Time for rref: {}'.format(ts0-ts1))

                #This is a consistency check: the rows of V2 should be linearly indpenedent and so
                #there should be no zero rows in V2_rref
                if len(pivs) != len(zero_singulars):
                    Vn = np.array(V2.evalf(),dtype=np.float)
                    np.save('debug_V.npy',Vn)
                    print(pivs)
                    print(zero_singulars)
                    print([S[i] for i in zero_singulars])
                    for i in range(len(zero_singulars)):
                        for j in range(Y.shape[1]):
                            V2[i,j] = round(V2[i,j],3)
                            V2_rref[i,j] = round(V2_rref[i,j],3)

                    sympy.pprint(V2)
                    sympy.pprint(V2_rref)
                    raise Warning('Issue with rref of V2, results are likely to be wrong!!!')

                if not opt.numX:
                    ais = []
                    for i in range(len(zero_singulars)):
                        ais.append(X.x[X.inds[pivs[i]]])

                    for j in range(Y.shape[1]):
                        if j not in pivs:
                            tmp = 0
                            print_debug = False
                            for i in range(len(zero_singulars)):
                                coeff = round(V2_rref[i,j],opt.round_prec)
                                tmp += coeff * ais[i]
                                #tmp += V2_rref[i,j] * ais[i]
                                if abs(coeff) > 10 or (abs(coeff) < 0.1 and abs(coeff) > 1e-3):
                                    print('Warning: Large or small coefficient {}, this can signify error'.format(coeff))
                                    print_debug = True
                            if debug or print_debug:
                                print('Substituting {} -> {}'.format(X.x[X.inds[j]],tmp))
                            X = X.subs(X.x[X.inds[j]],tmp)
                else:
                    for j in range(Y.shape[1]):
                        if j not in pivs:
                            sub_vec = np.zeros(Y.shape[1])
                            for i in range(len(zero_singulars)):
                                sub_vec[pivs[i]] = V2_rref[i,j]
                            X.subs(j,sub_vec)

            if debug_time:
                ts1 = time.perf_counter()
                print('Time for subs: {}'.format(ts1-ts0))
        else:
            raise Exception('Unknown algo_solve')

        #if debug_time:
        #    t1 = time.perf_counter()
        #    print('Time for constructing tensor ', t1-t2)
        #    print('')


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
    for i in range(X.dim2):
        perm[i] = i
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

def symmetrize_sym_inds(X,sym_inds,asym_inds,s_opt=None):

    def trans_func(X,perm,sign):

        X_T = X.copy0()

        for ind in X:
            ind_T = [0]*len(ind)
            for i in range(len(ind)):
                ind_T[i] = ind[perm[i]]
            ind_T = tuple(ind_T)
            X_T[ind_T] = sign * X[ind]

        ind_types = [0] * X.dim2
        for i in range(X.dim2):
            ind_types[i] = X.ind_types[perm[i]]
        X_T.ind_types = tuple(ind_types)

        for i in range(X.dim2):
            if X_T.ind_types[i] != X.ind_types[i]:
                X_T.reverse_index(i)

        return X_T

    if sym_inds is not None:

        perms_sym = []
        for si in sym_inds:
            perm = list(range(X.dim2))
            perm[si[0]] = si[1]
            perm[si[1]] = si[0]
            perms_sym.append(perm)
        X = symmetr(perms_sym,X,trans_func,1,s_opt)

    if asym_inds is not None:

        perms_asym = []
        for si in asym_inds:
            perm = list(range(X.dim2))
            perm[si[0]] = si[1]
            perm[si[1]] = si[0]
            perms_asym.append(perm)

        X = symmetr(perms_asym,X,trans_func,-1,s_opt)

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
