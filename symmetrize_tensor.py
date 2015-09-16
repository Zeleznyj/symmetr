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

def symmetr_ten(symmetries,op1,op2,order,proj=-1,T=None,debug=False,debug_Y=False,debug_time=False):
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
            #this transforms the tensor by the symmetry operation
            #even and odd parts transform differently
            if order % 2 == 0:
                X_trans = transform_tensor(X,sym,op1,op2,0,T=T,debug=debug)
            if order % 2 == 1:
                X_trans = transform_tensor(X,sym,op1,op2,1,T=T,debug=debug)

            
            if debug:
                print ''
                print 'Current form of the tensor:'
                print ''
                print X
                print ''
                print 'Transformed tensor:'
                print ''
                print X_trans

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
                print X
                print ''

    if debug:
        print 'Symmetrized tensor:'
        print X
        print ''
        print '======= End symmetrizing ======='

    return X

def print_tensor(ten):
    """
    Prints the expansion tensor in a nice form.

    Not tested for higher order than 1!!!
    """
    
    X = matrix(0,3)

    if ten.dim2 > 2:

        m = {}
        for i in range(3):
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

    X.pprint()

