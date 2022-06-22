# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Contains functions for obtaining list of symmetry operations and transformation matrix
based on user input.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import range
from past.utils import old_div
from copy import deepcopy
from itertools import product

from . import fslib
import sympy
from sympy import sympify as spf
from .groups import group_sym
from .noso import noso_syms
from .noso_new import NosoSymFinder
from .symmetry import findsym2sym, matsym2sym, Symmetry, create_I, create_P, create_T
from fractions import Fraction
import numpy as np
from mpmath import lu_solve
from numpy.linalg import norm

from .hsnf import smith_normal_form

def make_rational(mat):
    """Converts sympy matrix to a rational number form.

    This is useful since sympy can work exactly with rational numbers so there are no floating point accuracy issues.
    """

    ncols = mat.cols
    nrows = mat.rows
    mat_r = sympy.zeros(mat.rows,mat.cols)

    for i in range(nrows):
        for j in range(ncols):
            if type(mat[i,j]) == str and '\\' in mat[i,j]:
                mat_r[i,j] = sympy.sympify(Fraction(mat[i,j]))
            else:
                mat_r[i,j] = sympy.sympify(Fraction(float(mat[i,j])))

    return mat_r

def create_Tm(vec_a,vec_b,vec_c):

    T = sympy.zeros(3)

    T[0,0] = vec_a[0]
    T[1,0] = vec_a[1]
    T[2,0] = vec_a[2]

    T[0,1] = vec_b[0]
    T[1,1] = vec_b[1]
    T[2,1] = vec_b[2]

    T[0,2] = vec_c[0]
    T[1,2] = vec_c[1]
    T[2,2] = vec_c[2]

    T = make_rational(T)

    return T

def create_Ti(fin):

    T_i = sympy.zeros(3)

    inp_type = int(fin[2])

    if inp_type == 1:

        vec_1 = fin[3].split()
        vec_2 = fin[4].split()
        vec_3 = fin[5].split()

        T_i[0,0] =  spf(vec_1[0])
        T_i[1,0] =  spf(vec_1[1])
        T_i[2,0] =  spf(vec_1[2])

        T_i[0,1] =  spf(vec_2[0])
        T_i[1,1] =  spf(vec_2[1])
        T_i[2,1] =  spf(vec_2[2])

        T_i[0,2] =  spf(vec_3[0])
        T_i[1,2] =  spf(vec_3[1])
        T_i[2,2] =  spf(vec_3[2])

    if inp_type == 2:

        a,b,c,al,bet,gam = fin[3].split()
        al = al+'*2*pi/360'
        bet = bet+'*2*pi/360'
        gam = gam+'*2*pi/360'
        gam2 = 'acos((cos({gam})-cos({bet})*cos({al}))/(sin({al})*sin({al})))'.format(gam=gam,al=al,bet=bet)        

        
        T_i[0,0] =  spf('{a}*sin({gam2})*sin({bet})'.format(a=a,gam2=gam2,bet=bet))
        T_i[1,0] =  spf('{a}*cos({gam2})*sin({bet})'.format(a=a,gam2=gam2,bet=bet))
        T_i[2,0] =  spf('{a}*cos({bet})'.format(a=a,gam2=gam2,bet=bet))

        T_i[0,1] =  spf(0)
        T_i[1,1] =  spf('{b}*sin({al})'.format(b=b,al=al))
        T_i[2,1] =  spf('{b}*cos({al})'.format(b=b,al=al))

        T_i[0,2] =  spf(0)
        T_i[1,2] =  spf(0)
        T_i[2,2] =  spf(c)
    
    return T_i

def is_hex(lines):
    for i,line in enumerate(lines):
        if 'Values of a,b,c,alpha,beta,gamma:' in line:
            pos = i
    angles = lines[pos+1].split()[3:6]
    hexag= False
    for angle in angles:
        if int(round(float(angle))) != 90:
            if int(round(float(angle))) == 120:
                hexag = True
            else:
                for i,line in enumerate(lines):
                    if 'Space Group' in line:
                        sg = line.split()[-1]
                if sg == 'P1' or sg == 'P-1':
                    hexag = False
                    break
                raise Exception('one of the angles in findsym output is neither 90 nor 120.')
    return hexag

def get_syms(opt):
    if opt['inp']:
        #runs findsym and reads the output
        lines = fslib.run_fs(opt['inp'])
        syms = fslib.r_sym(lines,num_prec = opt['pos_prec'])

    if opt['group']:
        atom = -1
        _,syms=group_sym(opt['group'],debug=False)

    #this selects some of the symmetries 
    if opt['syms_sel'] != -1:

        syms_new = []
        for i in range(len(syms)):
            if i+1 in opt['syms_sel']:
                syms_new.append(syms[i])

        syms = syms_new

    syms_g = []
    for sym in syms:
        syms_g.append(findsym2sym(sym))
    syms = syms_g

    if opt['simplify_syms']:
        idxs,syms = simplify_symmetry_operations(syms,
                                                 generators=opt['generators'],
                                                 remove_P=opt['remove_P'],
                                                 remove_T=opt['remove_T']
                                                 )
        idxs1 = [i+1 for i in idxs]
        print('Using symmetry operations: {}'.format(idxs1))

    return syms

def get_T(opt,nonmag=False):
    """
    Returns transformation matrix from the conventional coordinate system (used by findsym)
    to the user selected basis.

    There are two different behaviors:
        nonmag=False: conventional coordinate system for the magnetic system
        nonmag=True: conventional coordinate system for the nonmagnetic system
    """
    if opt['inp'] is not None:
        #runs findsym and reads the output
        fin = fslib.read_fs_inp(opt['inp'],clean=False)
        lines = fslib.run_fs(opt['inp'])
        if nonmag:
            lines_nm = fslib.run_fs_nonmag(opt['inp'])

        #reads the input
        #vec_a,b,c are needed to know the basis transformation
        #syms contain the symmetries in the form that is needed by symmetr
        if not nonmag:
            [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
        else:
            [vec_a,vec_b,vec_c] = fslib.r_basis(lines_nm)

        #construct the transformation matrix from the findsym basis to the selected basis
        if 'abc' == opt['basis']:

            T = sympy.Matrix(sympy.Identity(3))

        if 'i' == opt['basis']:

            T = create_Tm(vec_a,vec_b,vec_c)

        if 'cart' == opt['basis']:

            T_m = create_Tm(vec_a,vec_b,vec_c)

            T_i = create_Ti(fin)

            T = T_i*T_m

        if 'custom' == opt['basis']:
            
            #T_i is a transformation matrix from the input coordinate system to the cartesian one
            T_i = create_Ti(fin)
                    
            for i in range(len(fin)):
                if 'axes:' in fin[i]:
                    loc = i
            
            vec_1 = fin[loc+1].split()
            vec_2 = fin[loc+2].split()
            vec_3 = fin[loc+3].split()

            #T_c is the transformation matrix from the used-defined basis to the cartesian one
            T_c = sympy.zeros(3)

            T_c[0,0] =  spf(vec_1[0])
            T_c[1,0] =  spf(vec_1[1])
            T_c[2,0] =  spf(vec_1[2])

            T_c[0,1] =  spf(vec_2[0])
            T_c[1,1] =  spf(vec_2[1])
            T_c[2,1] =  spf(vec_2[2])

            T_c[0,2] =  spf(vec_3[0])
            T_c[1,2] =  spf(vec_3[1])
            T_c[2,2] =  spf(vec_3[2])

            normalize = True
            if normalize == True:
                norm = sympy.sqrt(T_c[0,0]**2 + T_c[1,0]**2 + T_c[2,0]**2)
                T_c[0,0] = old_div(T_c[0,0], norm)
                T_c[1,0] = old_div(T_c[1,0], norm)
                T_c[2,0] = old_div(T_c[2,0], norm)

                norm = sympy.sqrt(T_c[0,1]**2 + T_c[1,1]**2 + T_c[2,1]**2)
                T_c[0,1] = old_div(T_c[0,1], norm)
                T_c[1,1] = old_div(T_c[1,1], norm)
                T_c[2,1] = old_div(T_c[2,1], norm)

                norm = sympy.sqrt(T_c[0,2]**2 + T_c[1,2]**2 + T_c[2,2]**2)
                T_c[0,2] = old_div(T_c[0,2], norm)
                T_c[1,2] = old_div(T_c[1,2], norm)
                T_c[2,2] = old_div(T_c[2,2], norm)

            T_m = create_Tm(vec_a,vec_b,vec_c)

            T = T_c.inv()*T_i*T_m

        if 'abc_c' == opt['basis']:

            T = sympy.zeros(3)

            if not nonmag:
                abc = fslib.r_abc(lines)
            else:
                abc = fslib.r_abc(lines_nm)

            a = abc[0]
            b = abc[1]
            c = abc[2]
            al = abc[3]
            bet = abc[4]
            gam = abc[5]

            al = al+'*2*pi/360'
            bet = bet+'*2*pi/360'
            gam = gam+'*2*pi/360'
            gam2 = 'acos((cos({gam})-cos({bet})*cos({al}))/(sin({al})*sin({al})))'.format(gam=gam,al=al,bet=bet)        
            
            T[0,0] =  spf('{a}*sin({gam2})*sin({bet})'.format(a=a,gam2=gam2,bet=bet))
            T[1,0] =  spf('{a}*cos({gam2})*sin({bet})'.format(a=a,gam2=gam2,bet=bet))
            T[2,0] =  spf('{a}*cos({bet})'.format(a=a,gam2=gam2,bet=bet))

            T[0,1] =  spf(0)
            T[1,1] =  spf('{b}*sin({al})'.format(b=b,al=al))
            T[2,1] =  spf('{b}*cos({al})'.format(b=b,al=al))

            T[0,2] =  spf(0)
            T[1,2] =  spf(0)
            T[2,2] =  spf(c)

    if opt['group'] is not None:
        atom = -1
        hex_group,_=group_sym(opt['group'],debug=False)

        if 'i' == opt['basis'] or 'abc' == opt['basis']:
            print('Using the conventional coordinate system!')
            T = sympy.Matrix(sympy.Identity(3))

        if 'cart' == opt['basis']:

            print('Using a cartesian coordinate system')
            if hex_group:
                T = sympy.zeros(3)
                T[0,0] = 1
                T[0,1] = sympy.sympify(Fraction(-0.5))
                T[0,2] = 0
                T[1,0] = 0 
                T[1,1] = old_div(sympy.sqrt(3),2)
                T[1,2] = 0
                T[2,0] = 0
                T[2,1] = 0
                T[2,2] = 1
            else:
                T = sympy.Matrix(sympy.Identity(3))

    return T

def get_syms_nonmag(opt):
    lines = fslib.run_fs_nonmag(opt['inp'])
    syms = fslib.r_sym(lines)
    
    syms_g = []
    for sym in syms:
        syms_g.append(findsym2sym(sym))
    syms = syms_g

    if opt['syms_sel'] != -1:

        syms_new = []
        for i in range(len(syms)):
            if i+1 in opt['syms_sel']:
                syms_new.append(syms[i])

        syms = syms_new

    if opt['simplify_syms']:
        idxs,syms = simplify_symmetry_operations(syms,
                                                 generators=opt['generators'],
                                                 remove_P=opt['remove_P'],
                                                 remove_T=opt['remove_T']
                                                 )
        idxs1 = [i+1 for i in idxs]
        print('Using symmetry operations: {}'.format(idxs1))

    return syms

def check_nonmag_syms(syms, prec=1e-5):
    n_syms = len(syms)
    all_good = True

    syms_fix = []
    for sym in syms:
        syms_fix.append(sym.copy())

    if n_syms % 2 != 0:
        print('odd number of symmetries!!!')
        return False, None

    sym0 = syms[0]
    if sym0.has_T:
        print('first symmetry has time-reversal, that should not happen')
        return False, None

    found_shift = False
    for i, sym in enumerate(syms):
        if ((sym0.R - sym.R).norm() < prec and
                sym0.permutations == sym.permutations and
                sym.has_T):
            found_shift = True
            shift = i
            break
    if not found_shift:
        print('Time-reversal variant of the first symmetry not found!')
        return False, None

    for i, sym in enumerate(syms):
        if i + shift == len(syms):
            break
        if sym.has_T:
            continue
        sym_op = syms[i + shift]
        if not sym_op.has_T:
            print(i,shift)
            print(sym)
            print(sym_op)
            print("Shifted symmetry does not have time-reversal, don't know what to do!")
            return False, None
        if (sym.R - sym_op.R).norm() > prec or sym.permutations != sym_op.permutations:
            print('Findsym problem with non-relativistic symmetry {}, fixing.'.format(i + shift))
            all_good = False
            syms_fix[i + shift] = sym.copy()
            syms_fix[i + shift].has_T = True
            syms_fix[i + shift].Rs = -sym.Rs
    return all_good, syms_fix

def get_syms_noso(opt):
    fin_c = fslib.read_fs_inp(opt['inp'])
    mags = fslib.r_mag_fin(fin_c)
    lines = fslib.run_fs(opt['inp'])
    lines_nm = fslib.run_fs_nonmag(opt['inp'])
    [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
    [vec_a_nm, vec_b_nm, vec_c_nm] = fslib.r_basis(lines_nm)

    Tm = create_Tm(vec_a,vec_b,vec_c)
    Tnm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)
    Ti = create_Ti(fin_c)
    for mag in mags:
        for i in range(3):
            mag[i] = mag[i] / Ti[:,i].norm()
    mags_T = convert_vecs(mags,Ti)
    mags_T_np = []
    for mag in mags_T:
        mags_T_np.append(np.array(mag, dtype=float).squeeze())

    syms_nm = get_syms_nonmag(opt)

    for sym in syms_nm:
        sym.convert(Ti*Tnm,in_place=True)

    if len(mags) > len(syms_nm[0].permutations):
        perms = get_full_permutations(lines_nm,debug=False)
        syms_nm_full = []
        for i, sym in enumerate(syms_nm):
            for j in range(len(perms[i])):
                sym_n = sym.copy()
                sym_n.permutations = perms[i][j]
                syms_nm_full.append(sym_n)
        syms_nm = syms_nm_full

    #there is a bug in findsym, where sometimes in a non-magnetic case, there is a wrong symmetry
    #There should be a symmetry with time-reversal for every symmetry without time-reversal, but
    #sometimes this is not the case
    #This can fix it, lathough it's more of a hack and ideally we should transition to newer version
    #of findsym, where this shouldn't happen.
    #This has to be run after the get_full_permutations!
    syms_good, syms_nm_fix = check_nonmag_syms(syms_nm)
    if not syms_good:
        if syms_nm_fix is None:
            raise Exception('Problem with non-magnetic symmetries.')
        else:
            print('fixing')
            syms_nm = syms_nm_fix

    symfinder = NosoSymFinder(prec=opt['noso_prec'],moment_zero=opt['noso_moment_zero'],debug=opt['noso_debug'])
    syms_noso = symfinder.find_noso_syms(syms_nm,mags_T_np)

    #we convert to the magnetic findsym basis, since this is what's used in the code elsewhere
    #this could be avoided though since we will just transform straight back in most situations
    for sym in syms_noso:
        sym.convert(Tm.inv()*Ti.inv(),in_place=True)

    return syms_noso

def get_syms_noso_old(opt):

    fin_c = fslib.read_fs_inp(opt['inp'])
    mags = fslib.r_mag_fin(fin_c)
    lines = fslib.run_fs(opt['inp'])
    lines_nm = fslib.run_fs_nonmag(opt['inp'])
    [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
    [vec_a_nm,vec_b_nm,vec_c_nm] = fslib.r_basis(lines_nm)

    Tm = create_Tm(vec_a,vec_b,vec_c)
    Tnm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)
    mags_T = convert_vecs(mags,Tm.inv())
    syms_nm = get_syms_nonmag(opt)
    for sym in syms_nm:
        sym.convert(Tm.inv()*Tnm,in_place=True)

    perms = get_full_permutations(lines_nm)
    for i,sym in enumerate(syms_nm):
        sym.permutations = perms[i]

    hexag = is_hex(lines)
    syms_noso = noso_syms(syms_nm,mags_T,hexag,debug=opt['debug_noso'])

    if opt['syms_sel_noso'] != -1:
        syms_sel_noso = opt['syms_sel_noso'].split(',')
        syms_sel_noso2 = []
        for i in range(len(syms_sel_noso)):
            if '-' in syms_sel_noso[i]:
                s = syms_sel_noso[i].split('-')
                syms_sel_noso2 += list(range(int(s[0]),int(s[1])+1))
            else:
                syms_sel_noso2.append(int(syms_sel_noso[i]))

        syms_noso_new = []
        for i in range(len(syms_noso)):
            if i+1 in syms_sel_noso2:
                syms_noso_new.append(syms_noso[i])

        syms_noso = syms_noso_new

    return syms_noso

def get_mags(inp):

    fin_c = fslib.read_fs_inp(inp)
    mags = fslib.r_mag_fin(fin_c)

    return mags

def convert_vec(vec,T):
    vec_T = T*vec.T 
    vec_T = vec_T.T
    return vec_T

def convert_vecs(mags,T):
    """Converts a list of vectors by transform matrix T.
    """

    mags_T = []
    for i,mag in enumerate(mags): 
        mag_T = convert_vec(mag,T)
        mags_T.append(mag_T)

    return mags_T

def get_Tm(inp,nonmag=False):

    if not nonmag:
        lines = fslib.run_fs(inp)
    else:
        lines = fslib.run_fs_nonmag(inp)

    [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
    Tm = create_Tm(vec_a,vec_b,vec_c)

    return Tm

def get_Tm_fin(inp_list,nonmag=False):
    if not nonmag:
        lines = fslib.run_fs_fin(inp_list)
    else:
        lines = fslib.run_fs_fin(fslib.make_fsinp_nonmag(inp_list))

    [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
    Tm = create_Tm(vec_a,vec_b,vec_c)

    return Tm

def get_metric(opt,debug=False,nonmag=False):

    if opt['group'] is not None:
        #the metric tensor is trivial in the cartesian coordinate system
        G = sympy.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        opt2 = deepcopy(opt)
        T = get_T(opt2)
        G = T.inv() * G * T.inv().T

    else:
        if debug:
            print('')
            print('Get_metric debug output')
            print('================================')
        #read the basis vectors used by findsym
        #they are given in the basis used in the findsym input
        if not nonmag:
            lines = fslib.run_fs(opt['inp'])
        else:
            lines = fslib.run_fs_nonmag(opt['inp'])
        [vec_a,vec_b,vec_c] = fslib.r_basis(lines)

        vec_a = sympy.Matrix(vec_a)
        vec_b = sympy.Matrix(vec_b)
        vec_c = sympy.Matrix(vec_c)

        #create transformation matrix from the input basis to a cartesian basis
        fin = fslib.read_fs_inp(opt['inp'],clean=False)
        T_i = create_Ti(fin)

        #Transform to a cartesian coordinate system
        vec_a = T_i * vec_a
        vec_b = T_i * vec_b
        vec_c = T_i * vec_c

        #calculate the dual basis vectors
        V = vec_a.dot(vec_b.cross(vec_c))
        vec_A = old_div(vec_b.cross(vec_c),V)
        vec_B = old_div(vec_c.cross(vec_a),V)
        vec_C = old_div(vec_a.cross(vec_b),V)

        if debug:

            print('basis vectors')
            print(vec_a,vec_b,vec_c)
            print('dual basis vectors')
            print(vec_A,vec_B,vec_C)

        A_l = vec_a.row_join(vec_b)    
        A_l = A_l.row_join(vec_c)

        A_L = vec_A.row_join(vec_B)    
        A_L = A_L.row_join(vec_C)

        if debug:
            print('A_l')
            print(A_l)
            print('A_L')
            print(A_L)

        if opt['num_prec'] is None:
            G = A_l.LUsolve(A_L)
        else:
            A_ln = np.array(A_l).astype(np.float64)
            A_Ln = np.array(A_L).astype(np.float64)
            if debug:
                print('A_ln')
                print(A_ln.dtype)
                print('A_Ln')
                print(A_Ln.dtype)
            G = np.linalg.solve(A_ln,A_Ln)
            G = sympy.Matrix(G)

        if debug:
            print('The metrix tensor')
            sympy.pprint(G)
            print('The inverse of the metrix tensor')
            sympy.pprint(G.inv())

        if debug:
            print('metric tensor test')
            for i in range(3):
                a_I = 0*vec_a
                for j in range(3):
                    a_I += G.T[i,j]*A_l[:,j]
                print(a_I)

            print('Get metric debug output end')
            print('================================')
            print('')

    return G


def get_generators(syms):
    gens = [syms[1]]
    gens_i = [1]
    for i, sym in enumerate(syms):
        is_generated = False
        for gen1 in gens:
            if sym == gen1:
                is_generated = True
                break
            if sym == gen1.inv():
                is_generated = True
                break
            for gen2 in gens:
                if gen1 * gen2 == sym:
                    is_generated = True
                    break
                if gen1 * gen2.inv() == sym:
                    is_generated = True
                    break
            if is_generated:
                break
        if not is_generated:
            for gen1 in gens:
                for gen2 in gens:
                    for gen3 in gens:
                        if gen1 * gen2 * gen3 == sym:
                            is_generated = True
        if not is_generated:
            gens.append(sym)
            gens_i.append(i)

    return gens_i, gens

def remove_symmetry(syms,R_sym,mode=0):
    if not isinstance(R_sym,Symmetry):
        if R_sym == 'P':
            R_sym = create_P()
        elif R_sym == 'I':
            R_sym = create_I()
        elif R_sym == 'T':
            R_sym = create_T()
        else:
            raise Exception('Wrong R_sym format!')
    #If R_sym has no permutations we assume it doesn't permute
    #the atoms. This is correct for I and T, but not necessarily for P!!!
    if R_sym.permutations is None:
        R_sym.permutations = {}
        for key in syms[0].permutations:
            R_sym.permutations[key] = key
    syms_out_i = []
    syms_out = []
    if mode == 0:
        for i,sym in enumerate(syms):
            sym2 = sym * R_sym
            if sym2 not in syms_out:
                syms_out.append(sym)
                syms_out_i.append(i)
    elif mode == 1:
        for i,sym in enumerate(syms):
            if sym != R_sym:
                syms_out.append(R_sym)
                syms_out_i.append(i)
    return syms_out_i,syms_out

def simplify_symmetry_operations(syms,generators=True,remove_P=False,remove_T=False):

    idxs = list(range(len(syms)))
    syms_g = syms
    if remove_P:
        i1, syms_g = remove_symmetry(syms_g, 'P')
        idxs = [idxs[i] for i in i1]
    if remove_T:
        i1, syms_g = remove_symmetry(syms_g, 'T')
        idxs = [idxs[i] for i in i1]
    if generators:
        i1, syms_g = get_generators(syms_g)
        idxs = [idxs[i] for i in i1]

    return idxs,syms_g

def get_full_permutations(lines,prec=3,debug=False):
    """
    This function determines permutations for all atoms. This is necessary in cases, where
    the non-magnetic unit cell is smaller than the magnetic and then the conventional method
    only determines the permutation for a subset of atoms, which is then a problem for noso.

    The way this works is:
    1. Get a list of atoms in the input, which are given in teh findsym basis.
    2. Convert those to the crystallographic basis used in findsym output.
    3. Transform these by symmetries.
    4. Transform back to the findsym basis.
    5. Check the transformations.

    In addition, when the non-magnetic unit cell is smaller we need to add additional translations to
    the original symmetry operations. These are identified by the find_translations function and are
    added to the non-magnetic symmetry operations.

    Args:
        lines: the findsym output, as a list of lines
        prec: the rounding precision, 3 means 3 digits of rounding and then atomic positions
            closer than 1e-2 are taken as equal.

    Returns:
        out: list of lists such that out[i] contains all the permutations for the non-magnetic symmetry
        i, that is the permutation for the original non-magnetic symmetry and the permutatinos for all
        the translations
    """

    prec2 = 1e-4

    if debug:
        print('get_full_permutations starting')

    origin = np.array(fslib.r_origin(lines),dtype=float)
    [vec_a,vec_b,vec_c] = fslib.r_basis(lines)
    Tm = np.array(create_Tm(vec_a,vec_b,vec_c),dtype=float)
    Tmi = np.array(np.linalg.inv(Tm),dtype=float)

    if norm(np.rint(Tmi) - Tmi) > prec2:
        raise Exception('Non-integerer transformation matrix: rotation of coordinate system present.')

    ts = find_translations(Tm)
    if debug:
        print(ts)

    start = False
    end = False
    pos_o = []
    for line in lines:
        if 'Position and magnetic moment of each atom (dimensionless coordinates)' in line:
            start = True
            continue
        if start and '------------------------------------------' in line:
            end = True
            continue
        if start and not end:
            pos_o.append([float(i) for i in line.split()[1:4]])          

    if debug:
        print('pos_o')
        for p in pos_o:
            print(p)
        print('')
    
    pos_T = []
    for p in pos_o:
        pos_T.append(np.dot(Tmi,p[0:3]-origin))

    if debug:
        print('pos_T')
        for p in pos_T:
            print(p)
        print('')
        
    syms = fslib.r_sym(lines,syms_only=False)
    
    out = []
    for isym,sym in enumerate(syms):
        if debug:
            print('')
            print(isym)
            print(sym[1])
            print(sym[4])
            print({x[0]: x[1] for x in sorted(sym[4])})
        out_t = []
        for trans in ts:
            if debug:
                print(trans)
            pos_TS = []
            for pT in pos_T:
                pTS = []
                for i in range(3):
                    pTS.append(sym[1][i])
                    for j,s in enumerate(['x','y','z']):
                        pTS[i] = pTS[i].replace(s,str(pT[j]))
                pTS = np.array(sympy.sympify(pTS),dtype=float)+trans
                pos_TS.append(pTS)

            if debug:
                print('pos_TS')
                for p in pos_TS:
                    print(p)
                print('')

            pos_TSTi = []
            for p in pos_TS:
                pos_TSTi.append( np.round( np.dot( np.array(Tm,dtype=float), p) + origin, prec) % 1)

            pos_TSTi_test = []
            for p in pos_TS:
                pos_TSTi_test.append( np.dot( np.array(Tm,dtype=float), p) + origin)

            if debug:
                print('pos_TSTi')
                for p in pos_TSTi:
                    print(p)
                print('')

            if debug:
                print('pos_TSTi_test')
                for p in pos_TSTi_test:
                    print(p)
                print('')

            pos_or = np.round(pos_o,prec) % 1
            if debug:
                print('pos_or')
                print(pos_or)

            permutations = {}
            for i,p in enumerate(pos_TSTi):
                found = False
                for j,po in enumerate(pos_or):
                    if np.linalg.norm(p-po) < 1.5/(10**(prec)):
                        #if debug:
                        #    print(i,p)
                        #    print(j,po)
                        #    print('')
                        perm = j
                        found = True
                        break
                if not found:
                    print(p)
                    print(pos_or)
                    raise Exception('Permutation not found')
                else:
                    permutations[i+1] = perm+1
            #print(sym[4])
            if debug:
                print(permutations)
            out_t.append(permutations)
        out.append(out_t)
        
    return out

def get_integer_solution(A, b, debug=False):
    """
    Finds an integer solution to a system of linear equtions.
    If there is no solution returs None.
    """
    B, U, V = smith_normal_form(A)
    if debug:
        print(B)
        print(U)
        print(V)
        print(np.dot(U, np.dot(A, V)))
    c = np.dot(U, b)
    y = np.zeros(B.shape[1])

    has_solution = True
    for i in range(B.shape[0]):
        if B[i, i] == 0:
            continue
        if c[i] % B[i, i] == 0:
            y[i] = c[i] // B[i, i]
        else:
            has_solution = False
            break
    if not has_solution:
        return None
    else:
        return np.dot(V, y)

def find_translations(T):
    """
    For matrix T that defines a supercell it finds all lattice translations that needs to be added to the symmetries
    of the supercell. This is because these translations are Bravais lattice translations in the smaller unit cell, but
    not in the supercell.

    This works by trying different combinations of the small unit cell lattice vectors and trying if these are an
    integer linear combination of the large unit cell lattice vectors (that are defined by T) or if they differ from
    already identified translation by an integer linear combination of the large unit cell lattice vectors.

    If det(T) = n that means the supercell is n-times larger and that should mean there is n translations.
    """

    prec = 1e-4
    Ti_o = np.linalg.inv(T)
    Ti = np.rint(Ti_o).astype(int)

    if norm(Ti_o - Ti) > prec:
        raise Exception('Not integer matrix!')

    det = np.rint(np.linalg.det(Ti))
    abs_det = abs(int(det))

    n_trans = 0

    it = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5]
    itt = list(product(it, it, it))
    itt = sorted(itt, key=norm)

    A = Ti.copy()

    ts = [np.array([0, 0, 0])]

    i = 0
    # while len(ts) < abs(det):
    while len(ts) < abs_det and i < len(itt):
        i += 1
        t = np.array(itt[i])
        # print(i)
        # print(t)

        found = False
        for ti in ts:
            x = get_integer_solution(A, t - ti)
            if x is not None:
                found = True
                break
        if not found:
            ts.append(t)
    if len(ts) < abs_det:
        print(i,len(itt))
        print(abs_det)
        print(ts)
        raise Exception('Did not find all transformations.')
    return ts

def get_full_permutations_old(lines, prec=3, debug=False):
    """
    This function determines permutations for all atoms. This is necessary in cases, where
    the non-magnetic unit cell is smaller than the magnetic and then the conventional method
    only determines the permutation for a subset of atoms, which is then a problem for noso.

    The way this works is:
    1. Get a list of atoms in the input, which are given in teh findsym basis.
    2. Convert those to the crystallographic basis used in findsym output.
    3. Transform these by symmetries.
    4. Transform back to the findsym basis.
    5. Check the transformations.

    Args:
        lines: the findsym output, as a list of lines
        prec: the rounding precision, 3 means 3 digits of rounding and then atomic positions
            closer than 1e-2 are taken as equal.

    Returns:
        list of permutations for all symmetry operations
    """

    if debug:
        print('get_full_permutations starting')

    origin = np.array(fslib.r_origin(lines), dtype=float)
    [vec_a, vec_b, vec_c] = fslib.r_basis(lines)
    Tm = create_Tm(vec_a, vec_b, vec_c)
    Tmi = np.array(Tm.inv().evalf(), dtype=float)

    start = False
    end = False
    pos_o = []
    for line in lines:
        if 'Position and magnetic moment of each atom (dimensionless coordinates)' in line:
            start = True
            continue
        if start and '------------------------------------------' in line:
            end = True
            continue
        if start and not end:
            pos_o.append([float(i) for i in line.split()[1:4]])

    if debug:
        print('pos_o')
        for p in pos_o:
            print(p)
        print('')

    pos_T = []
    for p in pos_o:
        print('wtttf')
        pos_T.append(np.dot(Tmi, p[0:3] - origin))

    if debug:
        print('pos_T')
        for p in pos_T:
            print(p)
        print('')

    syms = fslib.r_sym(lines, syms_only=False)

    out = []
    for isym, sym in enumerate(syms):
        if debug:
            print('')
            print(isym)
            print(sym[1])
            print(sym[4])
            print({x[0]: x[1] for x in sorted(sym[4])})
        pos_TS = []
        for pT in pos_T:
            pTS = []
            for i in range(3):
                pTS.append(sym[1][i])
                for j, s in enumerate(['x', 'y', 'z']):
                    pTS[i] = pTS[i].replace(s, str(pT[j]))
            pTS = np.array(sympy.sympify(pTS), dtype=float)
            pos_TS.append(pTS)

        if debug:
            print('pos_TS')
            for p in pos_TS:
                print(p)
            print('')

        pos_TSTi = []
        for p in pos_TS:
            pos_TSTi.append(np.round(np.dot(np.array(Tm, dtype=float), p), prec) % 1)

        pos_TSTi_test = []
        for p in pos_TS:
            pos_TSTi_test.append(np.dot(np.array(Tm, dtype=float), p))

        if debug:
            print('pos_TSTi')
            for p in pos_TSTi:
                print(p)
            print('')

        if debug:
            print('pos_TSTi_test')
            for p in pos_TSTi:
                print(p)
            print('')

        pos_or = np.round(pos_o, prec) % 1
        if debug:
            print('pos_or')
            print(pos_or)

        permutations = {}
        for i, p in enumerate(pos_TSTi):
            found = False
            for j, po in enumerate(pos_or):
                if np.linalg.norm(p - po) < 1 / (10 ** (prec - 1)):
                    # if debug:
                    #    print(i,p)
                    #    print(j,po)
                    #    print('')
                    perm = j
                    found = True
                    break
            if not found:
                raise Exception('Permutation not found')
            else:
                permutations[i + 1] = perm + 1
        # print(sym[4])
        if debug:
            print(permutations)
        out.append(permutations)

    return out


