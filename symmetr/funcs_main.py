# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
import re
import sys
import os

import symmetrize
import symmetrize_exp as st
import fslib
import funcs
import find_eq
import symT
import mham
from tensors import tensor,matrix

from sympy import sympify as spf

#this finds the location of the main.py file and ads this location to the path where modules are searched
#this way the modules have to be present only in the install directory and not in the run directory
dirname, filename = os.path.split(os.path.abspath(__file__))
sys.path.append(str(dirname))

def same_op_sym(X,T,debug=False,debug_rename=False):
    """
    If op1==op2 then additional requirement is that even part of the tensor is symmetric
    and the odd part antisymmetric.

    Args:
        X[matrix,matrix]: the tensors to be symmetrized
        T: transformation matrix to some cartesian basis

    Returns:
        X_S: the symmetrized (antisymmetrized) part of X
    """

    if debug:
        print ''
        print 'matrix in the original system before (anti-)symmetrizing'
        print 'even part'
        X[0].pprint()
        print 'odd part'
        X[1].pprint()

    X_C  = funcs.convert_X(X,T,debug=debug_rename)

    if debug:
        print ''
        print 'matrix in the cartesian system before (anti-)symmetrizing'
        print 'even part'
        X_C[0].pprint()
        print 'odd part'
        X_C[1].pprint()

    X_C[0] = funcs.sym_part(X_C[0])
    X_C[1] = funcs.asym_part(X_C[1])

    if debug:
        print ''
        print 'matrix in the cartesian system after (anti-)symmetrizing'
        print 'even part'
        X_C[0].pprint()
        print 'odd part'
        X_C[1].pprint()

    X_S = funcs.convert_X(X_C,T.inv(),debug=debug_rename,ignore_ren_warning=True)

    if debug:
        print ''
        print 'matrix in the original system after (anti-)symmetrizing'
        print 'even part'
        X_S[0].pprint()
        print 'odd part'
        X_S[1].pprint()

    return X_S


def sym_linres(opt,printit=False):
    """
    Finds the response tensors for linear response based on the input arguments.

    All the arguments are stored in opt. printit controls whether the output is printed.
    The number of output arguments depends on input parameters!

    Args:
        opt (class options): stores all the input arguments. Only some are used here.
        printit (optional[boolean]): if true this prints the output
    Returns:
        X([X1,X2]): where X1 and X2 are the two components of the linear response
        X_2[X1_2,X2_2): the two componets for second atom. Only outputed if opt['atom2'] is set.
        C(class confs): The linear response tensors for all equivalent magnetic configurations.
            Only outputed if opt['equiv'] is set.
    """

    op1 = opt['op1']
    op2 = opt['op2']
    op3 = opt['op3']

    eo =  symmetrize.even_odd(op1,op2,op3)
    
    #the symmetry operations given in the basis used by findsym
    syms = symT.get_syms(opt)
    #transformation matrix from the basis used by findsym to the user defined basis
    T = symT.get_T(opt)
    if opt['noso']:
        syms_noso = symT.get_syms_noso(opt)

    #This is a hacky code that will need to be replaced in the future.
    if not opt['ig_op1eqop2']:
        if op1 == op2 and op3 == None:
            
            lines = fslib.run_fs(opt['inp'])
            (vec_a,vec_b,vec_c) = fslib.r_basis(lines)
            fin = fslib.read_fs_inp(opt['inp'])

            T_m = symT.create_Tm(vec_a,vec_b,vec_c)
            T_i = symT.create_Ti(fin)
            T_cart1 = T_i*T_m
            T_cart2 = T_i*T_m*T.inv()

    #If this option is set then the tensors are symmetrized in the findsym basis and then transformed
    if opt['transform_result']:
        if not opt['noso']:
            X = symmetrize.symmetrize_linres(syms,op1,op2,op3=op3,proj=opt['atom'],\
                    debug=opt['debug_sym'],debug_time=opt['debug_time'],debug_Y=opt['debug_symY'])
        if opt['noso']:
            X = symmetrize.symmetrize_linres(syms_noso,op1,op2,op3=op3,proj=opt['atom'],\
                    debug=opt['debug_sym'],sym_format='mat',debug_time=opt['debug_time'],debug_Y=opt['debug_symY'])

        if not opt['ig_op1eqop2']:
            if op1 == op2 and op3 == None:

                X = same_op_sym(X,T_cart1,opt['debug_op1eqop2'],opt['debug_rename'])

        if op3 == None:
            if not opt['no_rename']:
                X_T = funcs.convert_X(X,T,debug=opt['debug_rename'])
            else:
                X_T = funcs.convert_X(X,T,ren=False,debug=opt['debug_rename'])
        else:
            X_T = []
            X_T.append(funcs.convert_tensor_3op(X[0],T))
            X_T.append(funcs.convert_tensor_3op(X[1],T))

    #If transform_result is not set then the tensors are symmetrized direclty in the basis chosen by user
    else:
        if not opt['noso']:
            X_T = symmetrize.symmetrize_linres(syms,op1,op2,op3=op3,proj=opt['atom'],T=T,\
                    debug=opt['debug_sym'],debug_time=opt['debug_time'],debug_Y=opt['debug_symY'])
        if opt['noso']:
            X_T = symmetrize.symmetrize_linres(syms_noso,op1,op2,op3=op3,proj=opt['atom'],T=T,\
                    sym_format='mat',debug=opt['debug_sym'],debug_time=opt['debug_time'],debug_Y=opt['debug_symY'])

        if not opt['ig_op1eqop2']:
            if op1 == op2 and op3 == None:

                X_T = same_op_sym(X_T,T_cart2,opt['debug_op1eqop2'],opt['debug_rename'])

    if printit:
        if op3 == None:
            print 'Symmetry restricted form of the tensor, %s part' % eo[0]
            X_T[0].pprint()
            if opt['latex']:
                X_T[0].pprint(latex=opt['latex'])
            print ''

            print 'Symmetry restricted form of the tensor, %s part' % eo[1]
            X_T[1].pprint()
            if opt['latex']:
                X_T[1].pprint(latex=opt['latex'])
            print ''

        else:
            spins=['x','y','z']
            print '%s part:' % eo[0]
            for i in range(3):
                print 'op1=',spins[i]
                X_T[0].reduce(0,i).pprint(latex=opt['latex'])
                print ''

            print '%s part:' % eo[1]
            for i in range(3):
                print 'op1=',spins[i]
                X_T[1].reduce(0,i).pprint(latex=opt['latex'])
                print ''

    #If atom2 is set then the response tensors for atom1 are transformed to atom2
    if opt['atom2'] != -1 and op3 == None:

        X_T_2 = symmetrize.symmetr_AB(syms,X_T,op1,op2,opt['atom'],opt['atom2'],T=T)

        if X_T_2 == None:
            print 'no relation with atom %s found' % opt['atom2']
        else:
            if printit:
                print 'Symmetry restricted form of the tensor, %s part, atom %s' % (eo[0],opt['atom2'])
                X_T_2[0].pprint()
                if opt['latex']:
                    X_T_2[0].pprint(latex=True)
                print 'Symmetry restricted form of the tensor, %s part, atom %s' % (eo[1],opt['atom2'])
                X_T_2[1].pprint()
                if opt['latex']:
                    X_T_2[1].pprint(latex=True)
                print ''

    #if equiv is set then we transform the tensor to all equivalent magnetic configurations
    if opt['equiv']:
        fin_c = fslib.read_fs_inp(opt['inp'])
        mags = fslib.r_mag_fin(fin_c)
        lines = fslib.run_fs(opt['inp'])
        lines_nm = fslib.run_fs_nonmag(opt['inp'])
        [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
        [vec_a_nm,vec_b_nm,vec_c_nm] = fslib.r_basis(lines_nm)
        syms_nm = symT.get_syms_nonmag(opt)

        #reads the magnetic moments from the input file
        mags = fslib.r_mag_fin(fin_c)
        #transformation matrix from the magnetic basis to the input one
        Tm = symT.create_Tm(vec_a,vec_b,vec_c)

        #transformation matrix from the nonmagnetic basis to the input one
        Tnm = symT.create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)

        #convert the magnetic moments to the selected basis
        #T*Tm.inv() is a matrix that transforms the magnetic moments to the selected basis because:
        #Tm.inv() converts to the magnetic basis and then T converts to the selected
        mags_T = funcs.convert_vecs(mags,T*Tm.inv())
        #convert the non-magnetic symmetry operations to the chosen basis
        syms_nm_T = funcs.convert_sym_mat(syms_nm,T*Tm.inv()*Tnm,sym_format='findsym')

        #this outputs all the equaivalen configurations
        #C is a conf class, it contains both the configurations and the transformed tensors
        C = find_eq.find_equiv(X_T,opt['op1'],opt['op2'],opt['atom'],syms_nm_T,mags_T,op3=opt['op3'],\
                debug=opt['debug_equiv'])
        if printit:
            print ''
            print 'Equivalent configurations:'
            C.pprint(eo,latex=opt['latex'])

    #based on the chosen options different variables are returned
    if not opt['equiv']:
        if opt['atom2'] == -1:
            return X_T
        else:
            return X_T,X_T_2
    else:
        if opt['atom2'] == -1:
            return X_T,C
        else:
            return X_T,X_T_2,C

def sym_exp(opt,printit=False):
    """
    Finds the tensor describing expansion term of linear response in the direction of magnetization.

    Args:
        opt (class options): stores all the input arguments. Only some are used here.
        printit (optional[boolean]): if true this prints the output
    """

    if opt['group']:
        print '!!!The input group must be one of the nonmagnetic point groups, otherwise the ouput will be wrong.!!!' 
        syms_nm = symT.get_syms(opt)
        T = symT.get_T(opt)
    else:
        T = symT.get_T(opt,nonmag=True)
        syms_nm = symT.get_syms_nonmag(opt)

    X = st.symmetrize_exp(syms_nm,opt['op1'],opt['op2'],opt['exp'],proj=opt['atom'],T=T,\
            debug=opt['debug_sym'],debug_Y=opt['debug_symY'],debug_time=opt['debug_time'])

    if printit:
        st.print_tensor(X)
        if opt['latex']:
          st.print_tensor(X,latex=opt['latex'])

    return X

def sym_res(opt,printit=False):
    """A wrapper function that returns the appropriate response tensor based on the input options.

    Args:
        opt (class options): stores all the input arguments. Only some are used here.
        printit (optional[boolean]): if true this prints the output
    """
    if opt['exp'] == -1:
        return sym_linres(opt,printit=printit)
    else:
        return sym_exp(opt,printit=printit)

def sym_mham(opt,printit=False):
    T = symT.get_T(opt,nonmag=True)
    syms = symT.get_syms_nonmag(opt)
    if opt['transform_syms']:
        H_T = mham.sym_mag_ham(opt['sites'],syms,T=T,debug=opt['debug_sym'])
    else:
        H = mham.sym_mag_ham(opt['sites'],syms,T=None,debug=opt['debug_sym'])
        H_T = mham.convert_mag_ham(H,T)
    if opt['equiv']:
        H_E = mham.equiv(H_T,opt['sites'],syms,T)
    if printit:
        if H_T.dim2 == 2:
            print 'Hamiltonian term in matrix form:'
            H_T.pprint(latex=opt['latex'])
            print ''
        mham.print_Ham(H_T,opt['sites'],latex=opt['latex'])
        if opt['equiv']:
            print ''
            print 'Hamiltonian terms for all equivalent combinations of sites:'
            for sites in H_E:
                print str(sites)+':'
                mham.print_Ham(H_E[sites],sites,latex=opt['latex'])
                print ''
    if not opt['equiv']:
        return H_T
    else:
        return H_T,H_E

