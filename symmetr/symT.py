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

from . import fslib
import sympy
from sympy import sympify as spf
from .groups import group_sym
from .noso import noso_syms
from .symmetry import findsym2sym,matsym2sym
from fractions import Fraction

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
                sys.exit('one of the angles in findsym output is neither 90 nor 120.')
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
        syms_sel = opt['syms_sel'].split(',')
        syms_sel2 = []
        for i in range(len(syms_sel)):
            if '-' in syms_sel[i]:
                s = syms_sel[i].split('-')
                syms_sel2 += list(range(int(s[0]),int(s[1])+1))
            else:
                syms_sel2.append(int(syms_sel[i]))

        syms_new = []
        for i in range(len(syms)):
            if i+1 in syms_sel2:
                syms_new.append(syms[i])

        syms = syms_new

    syms_g = []
    for sym in syms:
        syms_g.append(findsym2sym(sym))
    syms = syms_g

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

    return syms

def get_syms_noso(opt):

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

def get_metric(opt,debug=False):

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
        lines = fslib.run_fs(opt['inp'])
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

        G = A_l.LUsolve(A_L)
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

