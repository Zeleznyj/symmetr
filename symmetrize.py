import re
import copy
import sys
import math

import sympy
import numpy as np

from tensors import matrix, mat2ten
from read import transform_position
from funcs import *
from conv_index import *

class params_trans:
    def __init__(self,op1,op2,op3,l,T=None,sym_format='findsym'):
        self.op1 = op1
        self.op2 = op2
        self.op3 = op3
        self.l = l
        self.T = T
        self.sym_format = sym_format

def symmetr(syms,X,trans_func,params,debug=False,debug_time=False,debug_Y=False):
    """
    This symmetrizes a tensor X given a list of symmetries and a transformation function.

    This function should be quite general and is now used for all symmetrizing.

    Args:
        syms: list of symmetry operations
        X: tensor - must be a tensor class
        trans_func: function that transforms the tensor X using symmetry sym
            trans_func must work in the following way:
            X_trans = trans_func(sym,X,params)
        params: parameters to be sent to function trans_func

    Returns:
        X_trans: the symmetry restricted form of tensor X
    """

    if debug:
        print ''
        print '======= Starting symmetrizing ======='

    #we do a loop over all symmetry operations, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
    #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
    for sym in syms:
        
        if debug:
            print 'Symmetry:' 
            print sym
            print ''

        X_trans = trans_func(X,sym,params)
            
        if debug:
            print ''
            print 'Current form of the tensor:'
            print ''
            X.pprint()
            print ''
            print 'Transformed tensor:'
            print ''
            X_trans.pprint()

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
            print 'Time for constructing Y: ', t2-t1
            t1 = time.clock()

        #this transforms the matrix into the Reduced row echelon form
        #piv are the indeces o the pivot columns
        [rref,piv] = Y.rref()        

        if debug_time:
            t2 = time.clock()
            print 'Time for reducing Y to reduced row echelon form: ', t2-t1

        if debug_Y:
            print ''
            print 'Matrix representing the linear equation system that has to be satisfied: (right hand side is zero)'
            sympy.pprint(Y)
            print ''
            print 'Reduced row echelon form and indeces of the pivot columns:'
            sympy.pprint([rref,piv])
            print ''

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
                print ''
                print 'considering pivot ', i,j

            tmp = 0
            #now we just make use of the linear equation that holds for this pivot
            #keep in mind that the rows are in reversed order
            for ll in range(j+1,X.dim1**X.dim2):
                tmp = tmp - rref[i,ll]*X.x[rev_inds[ll]]
            X = X.subs(X.x[rev_inds[j]],tmp)

            if debug:
                print 'substituting ',
                sympy.pprint(X.x[rev_inds[j]])
                print ' for ',
                sympy.pprint(tmp)
                print ''

        if debug:
            print 'Current form of the tensor:'
            X.pprint()
            print ''

    if debug:
        print 'Symmetrized tensor:'
        X.pprint()
        print ''
        print '======= End symmetrizing ======='

    return X

def symmetrize_linres(symmetries,op1,op2,op3=None,proj=-1,debug=False,debug_time=False,debug_Y=False,\
        T=None,sym_format='findsym'):
    """
    Returns a symmetrical form of a response matrix for a given atom and given list of symmetries.

    Args:
        symmetries: A list of symmetry operations. In a format outputted by read.py.
        op1 (string): First operator type.
            's' for spin
            'v' for velocity operator
            'x' for position
        op2 (string): Second operator type. 
        proj (Optional[int]): Determines projection on atom. Defaults to -1.
            If set to -1, there is no projection.
            If set to positive integer, it determines atom number.
        debug (Optional[boolean]): Defaults to false. If set to true, additional debug output is printed.
        T (Default[sympy matrix): A linear response matrix. Defaults to None
            If it is set, then the linear response matrix is used to transform the symmetry operations.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.
            Primarily for debugging.

    Outputs:
        X ([X[0],X[1]]): A list which contains symmetrized form of the even and the odd part of the linear response tensor.
    """

    #this defines starting response matrix
    #we repeat it twice, once for the even part and once for the odd part
    #matrices are stored as tensor class, which is made using sympy
    #matrix is a sublcass of tensor
    
    if op3 == None:
        X1 = matrix('s',3)
        X2 = matrix('s',3)
    else:
        X1 = tensor('s',3,3)
        X2 = tensor('s',3,3)

    X = []
    X.append(X1)
    X.append(X2)

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

        if take_sym:
            syms_sel.append(sym)


    for l in range(2):
        params = params_trans(op1,op2,op3,l,T,sym_format)
        X[l] = symmetr(syms_sel,X[l],transform_tensor_params,params,debug=debug,\
                debug_time=debug_time,debug_Y=debug_Y)

    return X

def symmetr_AB(syms,X,op1,op2,atom1,atom2,T=None):
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
        if sym_type(atom1,sym) == atom2 and not found:
            found = True
            for l in range(2):
                X_trans.append(transform_matrix(X[l],sym,op1,op2,l,T=T))

    if found:
        return X_trans
    else:
        return None

