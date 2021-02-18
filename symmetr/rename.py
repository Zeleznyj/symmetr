from __future__ import print_function
from __future__ import absolute_import
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from builtins import str
from builtins import range
import re

import sympy
import numpy
from .tensors import matrix, mat2ten, Tensor
from .conv_index import *


def should_rename(X,X_t):
    #looks at the element and the transformation of the element and decides if the element should be renamed

    if X - X_t == 0:
        return False
    else: 
        inds = re.findall("\[([0-9]), ([0-9])\]",str(X))
        if len(inds) > 1:
            return False
        elif len(inds) == 0:
            sys.exit('error in should_rename')
        else:
            ret = True
            inds_t = re.findall("\[([0-9]), ([0-9])\]",str(X_t))
            for ind in inds_t:
                if convert_index(inds[0][0],inds[0][1]) < convert_index(ind[0],ind[1]):
                    ret = False
            return ret


#renames the matrix so that it has the simplest form possible
#it looks for the relation between components so no information is lost
#debug is an optional parameter, if it's true, then the routine outputs lots of information
def rename(X,name,debug=False,ignore_ren_warning=False):

    ninds = len(list(set(re.findall(r'x[0-9]+',sympy.srepr(X)))))

    #v contains the symbolic variables:
    V = matrix('s',3,name)
    #Y will be the renamed matrix
    Y = matrix(0,3)

    if debug:
        print('')
        print('======= Start renaming =======')
        print('Input matrix:')
        sympy.pprint(X.mat())

        print('')
        print('independent components:')
        print(list(set(re.findall(r'x[0-9]+',sympy.srepr(X)))))
        print('')

    #a loop over all components of  the input matrix
    for i in range(3):
        for j in range(3):

            #if a components is zero we do nothing with it
            if X[i,j] == 0:
                Y[i,j] = 0
            else:

                if debug:
                    print('')
                    print('component: ', i, j, ' = ', end=' ')
                    sympy.pprint(X[i,j])

                pos = convert_index(i,j)
                if pos == 0: #we always rename the first component of the matrix (unless it's zero)
                    Y[i,j] = V.x[i,j]
                    if debug:
                        print('renaming to:', ' ', end=' ')
                        sympy.pprint(V.x[i,j])
                # if we don't have the first component we try if the component is not a linear combination
                # of the previous
                #components, if it's not we give it a new name, if it is then we set it equal to that combination
                else: 
                    #to find if the component is a linear combination of a previous components we set a linear equations system
                    #matrix of that system will be Z, last column of Z is the right hand side
                    Z = sympy.zeros(9,pos+1)
                    # al loop over all components of the Z matrix
                    for p in range(pos+1):
                        for o in range(9):

                            #we set the variable that corresponds to the column to 1 and the rest to 0
                            tmp = X[inconvert_index(p)[0],inconvert_index(p)[1]]
                            for m in range(9):
                                if m == o:
                                    tmp = tmp.subs(V.x[inconvert_index(m)[0],inconvert_index(m)[1]],1)
                                else:
                                    tmp = tmp.subs(V.x[inconvert_index(m)[0],inconvert_index(m)[1]],0)

                            Z[o,p] = tmp
                    
                    if debug:
                        print('linear equations system: (last column is the right hand side)')
                        sympy.pprint(Z)

                    #this will give a solution to the linear equation system
                    # if there is no solution it outputs none
                    #there is a sympy routine for this, but it outputs a general solution that may contain parameters
                    #solve_lin routines sets all arbitrary parameters to 0, ie if there is infinitely many solutions it outputs just one
                    solution = solve_lin(Z)
                    if solution:
                        for ii in range(len(solution)):
                            if abs(solution[ii]) < 1.e-14:
                                solution[ii] = 0
                        sol = False

                    if debug:
                        print('solution of the system:')
                        print(solution)

                    #if the component is not a linear combination of previous components we rename it
                    if not solution:
                        Y[i,j] = V.x[i,j]

                        if debug:
                            print('renaming to:', ' ', end=' ')
                            sympy.pprint(V.x[i,j])

                    #if it is we set it equal to the linear combination
                    else:
                        for r in range(len(solution)):
                            Y[i,j] = Y[i,j] + solution[r] * Y[inconvert_index(r)[0],inconvert_index(r)[1]]
                        if debug:
                            print('renaming to:', ' ', end=' ')
                            sympy.pprint(Y[i,j])
                    

    if debug:
        print('')
        print('Input matrix:')
        sympy.pprint(X.mat())
        print('')
        print('renamed matrix:')
        sympy.pprint(Y.mat())
        print('')
        print('======= end rename =======')
        print('') 

    ninds_new = len(list(set(re.findall(r'x[0-9]+',sympy.srepr(Y)))))

    if ninds != ninds_new and not ignore_ren_warning:
        print('!WARNING! Problem in renaming tensor components. Try --no-rename.') 
        print('old number of independent components:', ninds)
        print(list(set(re.findall(r'x[0-9]+',sympy.srepr(X)))))
        print('new number of independent components:', ninds_new)
        print(list(set(re.findall(r'x[0-9]+',sympy.srepr(Y)))))

    return Y


#outputs a solution to the linear equation system represented by matrix Z
#last column of Z is the right hand side
# if there is more than one solution, this chooses one solution by setting all arbitraty parameters to 0
def solve_lin(Z):

    #we need to loop over variables, this doesn't seem to be possible with sympy so this is a way around that
    b0 = sympy.Symbol('b0')
    b1 = sympy.Symbol('b1')
    b2 = sympy.Symbol('b2')
    b3 = sympy.Symbol('b3')
    b4 = sympy.Symbol('b4')
    b5 = sympy.Symbol('b5')
    b6 = sympy.Symbol('b6')
    b7 = sympy.Symbol('b7')

    def convert_b(i):
        if i == 0:
            return b0
        if i == 1:
            return b1
        if i == 2:
            return b2
        if i == 3:
            return b3
        if i == 4:
            return b4
        if i == 5:
            return b5
        if i == 6:
            return b6
        if i == 7:
            return b7

    Zr = sympy.zeros(Z.rows,Z.cols)
    for i in range(Z.rows):
        for j in range(Z.cols):
            Zr[i,j] = Z[i,j]

    #this solves  the system, not all the variabes may be needed, but luckily sympy doesn't complain about that
    solution = sympy.solve_linear_system(Zr,b0,b1,b2,b3,b4,b5,b6,b7,minimal=True,warn=True)

    #number of variables of the linear equation system
    dim = len(Z[1,:])-1
    #the list that will contain the solution
    out = []
    #if there is no solution we return None
    if not solution:
        return None
    else:
        #if there is a solution, we look at every component of the solution and set the arbitrary parameters to 0
        #the arbitrary parameters are not in the solution dict, that's how we find them
        for i in range(dim):
            x = convert_b(i)
            if x in solution:
                out_temp = solution[x]
                #we go over all variables and check if they are arbitrary
                for j in range(dim):
                    y = convert_b(j)
                    #if they are arbitrary we set them to zero
                    if y not in solution:
                        out_temp = out_temp.subs(y,0)
            else:
                out_temp = 0
            out.append(out_temp)

    return out






    
