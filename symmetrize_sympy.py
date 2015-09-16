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


def symmetr(symmetries,op1,op2,proj=-1,debug=False,T=None):
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
    
    X1 = matrix('s',3)
    X2 = matrix('s',3)
    X = []
    X.append(X1)
    X.append(X2)

    if debug:
        print ''
        print '======= Starting symmetrizing ======='

    #we do a loop over all symmetry operations, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
    #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
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
            #we do everything separately for the even and odd parts
            #most things are the same, the only difference in the physics is that when time-reversal is present, odd part  transformation has
            #minus compared to the even part transformation
            for l in range(2):
                if debug:
                    if l == 0:
                        print ''
                        print 'Even part'
                        print ''
                    if l == 1:
                        print ''
                        print 'Odd part'
                        print ''

                #this transforms the matrix by the symmetry operation
                X_trans = transform_matrix(X[l],sym,op1,op2,l,debug=debug,T=T)
                
                if debug:
                    print ''
                    print 'Current form of the matrix:'
                    print ''
                    sympy.pprint(X[l].mat())
                    print ''
                    print 'Transformed matrix:'
                    print ''
                    sympy.pprint(X_trans.mat())

                #the matrix must be equal to the transformed matrix, this give us a system of 9 linear equations
                #matrix Y is a matrix that represents this system, ie the system X-X_trans = 0
                #we reverse the order of the rows - ie the first row corresponds to x[2,2] and last to x[0,0]
                # it doesn't really matter but the results are more natural this way
                Y = matrix(0,9)

                #we do a loop over all rows of the matrix Y - ie over all linear equations
                for i in range(3):
                    for j in range(3):

                        #convert_index transforms an index in a 3x3 matrix into an index in a 1x9 vector form
                        m = convert_index(i,j)

                        #a loop over all columns of matrix Y
                        for k in range(3):
                            for ll in range(3):

                                #again converts an index from 3x3 matrix form to the 1x9 vector form, but in this case in the reversed order
                                n = convert_index_rev(k,ll)

                                #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                                Y_p = X[l][i,j]-X_trans[i,j]
                                for o in range(3):
                                    for p in range(3):
                                        if o == k and p == ll:
                                            Y_p = Y_p.subs(X[l].x[o,p],1)
                                        else:
                                            Y_p = Y_p.subs(X[l].x[o,p],0)

                                Y[m,n] = Y_p

                #this transforms the matrix into the Reduced row echelon form
                #piv are the indeces o the pivot columns
                [rref,piv] = Y.mat().rref()        

                if debug:
                    print ''
                    print 'Matrix representing the linear equation system that has to be satisfied: (right hand side is zero)'
                    sympy.pprint(Y.mat())
                    print ''
                    print 'Reduced row echelon form and indeces of the pivot columns:'
                    sympy.pprint([rref,piv])
                    print ''

                #a loop over all the pivots: it's the pivots that give interesting information
                for j in list(reversed(piv)):

                    
                    #find the row of pivot j
                    found = False
                    i = 8
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
                    for ll in range(j+1,9):
                        tmp = tmp - rref[i,ll]*X[l].x[inconvert_index_rev(ll)[0],inconvert_index_rev(ll)[1]]
                    X[l] = X[l].subs(X[l].x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]],tmp)

                    if debug:
                        print 'substituting ',
                        sympy.pprint(X[l].x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]])
                        print ' for ',
                        sympy.pprint(tmp)
                        print ''


                if debug:
                    print 'Current form of the matrix:'
                    sympy.pprint(X[l].mat())
                    print ''

    if debug:
        print ''
        print 'Symmetrized matrix even part:'
        sympy.pprint(X[0].mat())
        print ''
        print 'Symmetrized matrix odd part:'
        sympy.pprint(X[1].mat())
        print ''
        print '======= End symmetrizing ======='


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

