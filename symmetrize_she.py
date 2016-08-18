import re
import copy
import sys
import math
import itertools

import sympy
import numpy as np

from tensors import matrix, mat2ten
from read import transform_position
from funcs import *
from conv_index import *


def symmetr_3op(symmetries,op1,op2,op3,proj=-1,debug=False,T=None):
    """
    Returns a symmetrical form of a response tensor for a 3 operator linear response formula
            and given list of symmetries.

    Args:
        symmetries: A list of symmetry operations. In a format outputted by read.py.
        op1 (string): First operator type.
            's' for spin
            'v' for velocity operator
            'x' for position
        op2 (string): Second operator type. 
        op3 (string): Third operator type. 
        proj (Optional[int]): Determines projection on atom. Defaults to -1.
            If set to -1, there is no projection.
            If set to positive integer, it determines atom number.
            !!!I'm not sure if this option really has any meaning in this context.!!!
        debug (Optional[boolean]): Defaults to false. If set to true, additional debug output is printed.
        T (Optional[sympy matrix): A linear response matrix. Defaults to None
            If it is set, then the linear response matrix is used to transform the symmetry operations.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.


    Outputs:
        X ([X[0],X[1]]): A list which contains symmetrized form of the even and the odd part of the linear response tensor.
    """

    #this defines starting response matrix
    #we repeat it twice, once for the even part and once for the odd part
    #matrices are stored as tensor class, which is made using sympy
    #matrix is a sublcass of tensor
    
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
                X_trans = transform_tensor_3op(X[l],sym,op1,op2,op3,l,debug=debug,T=T)
                
                if debug:
                    print ''
                    print 'Current form of the tensor:'
                    print ''
                    print X[l]
                    print ''
                    print 'Transformed tensor:'
                    print ''
                    print X_trans

                #the tensor must be equal to the transformed tensor, this give us a system of 27 linear equations
                #matrix Y is a matrix that represents this system, ie the system X-X_trans = 0
                #we reverse the order of the rows - ie the first row corresponds to x[2,2,2] and last to x[0,0,0]
                # it doesn't really matter but the results are more natural this way
                Y = matrix(0,27)

                #we do a loop over all rows of the matrix Y - ie over all linear equations
                for i,j,q in itertools.product(range(3),range(3),range(3)):

                    #convert_index transforms an index in a 3x3x3 matrix into an index in a 1x27 vector form
                    m = convert_index_3(i,j,q)

                    #a loop over all columns of matrix Y
                    for k,ll,r in itertools.product(range(3),range(3),range(3)):

                        #again converts an index from 3x3x3 matrix form to the 1x27 vector form, but in this case in the reversed order
                        n = convert_index_rev_3(k,ll,r)

                        #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                        Y_p = X[l][i,j,q]-X_trans[i,j,q]
                        for o,p,qq in itertools.product(range(3),range(3),range(3)):
                                if o == k and p == ll and r == qq:
                                    Y_p = Y_p.subs(X[l].x[o,p,qq],1)
                                else:
                                    Y_p = Y_p.subs(X[l].x[o,p,qq],0)

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
                    i = 26
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
                    for ll in range(j+1,27):
                        invindx = inconvert_index_rev_3(ll)
                        tmp = tmp - rref[i,ll]*X[l].x[invindx[0],invindx[1],invindx[2]]
                    invindx = inconvert_index_rev_3(j)
                    X[l] = X[l].subs(X[l].x[invindx[0],invindx[1],invindx[2]],tmp)

                    if debug:
                        print 'substituting ',
                        sympy.pprint(X[l].x[invindx[0],invindx[1],invindx[2]])
                        print ' for ',
                        sympy.pprint(tmp)
                        print ''


                if debug:
                    print 'Current form of the tensor:'
                    print X[l]
                    print ''

    if debug:
        print ''
        print 'Symmetrized tensor even part:'
        print X[0]
        print ''
        print 'Symmetrized tensor odd part:'
        print X[1]
        print ''
        print '======= End symmetrizing ======='


    return X

