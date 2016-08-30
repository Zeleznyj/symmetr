#!/usr/bin/python
import re
import sys
import os
import subprocess

import symmetrize_sympy
import symmetrize_she
import symmetrize_tensor as st
import read
from tensors import matrix, mat2ten
import funcs
from rename import rename
from groups import group_sym

import find_eq

import sympy
from sympy import sympify as spf
import numpy as np

from mpmath import cos as mcos
from mpmath import sin as msin
from mpmath import acos as macos
from mpmath import radians as mradians

from fractions import Fraction

import argparse

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

    T = funcs.make_rational(T)

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

def fs_nonmag(fin_c):
    
    #replaces the magnetic moments by 0
    start = False
    fin_cnm = []
    for i in range(len(fin_c)):
        if 'magnetic' in fin_c[i]:
            start = True
        if start:
            fin_cnm.append(re.sub(r'([0-9\.\-]+ +[0-9\.\-]+ +[0-9\.\-]+).+',r'\1 0 0 0',fin_c[i],count=1))
        else:
            fin_cnm.append(fin_c[i])
    
    #sends the nonmagnetic input file to findsym
    fs = subprocess.Popen([dirname+'/findsym'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    out_nm = fs.communicate(input=''.join(fin_cnm))[0]
    lines_nm = out_nm.split('\n')

    return lines_nm

def same_op_sym(X,T):
    """
    If op1==op2 then additional requirement is that even part of the tensor is symmetric
    and the odd part antisymmetric.

    Args:
        X[matrix,matrix]: the tensors to be symmetrized
        T: transformation matrix to some cartesian basis

    Returns:
        X_S: the symmetrized (antisymmetrized) part of X
    """

    X_C  = funcs.convert_X(X,T)
    X_C[0] = funcs.sym_part(X_C[0])
    X_C[1] = funcs.asym_part(X_C[1])
    X_S = funcs.convert_X(X_C,T.inv())

    return X_S

#this finds the location of the main.py file and ads this location to the path where modules are searched
#this way the modules have to be present only in the install directory and not in the run directory
dirname, filename = os.path.split(os.path.abspath(__file__))
sys.path.append(str(dirname))

#parses the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('op1', help='Type of the first operator')
parser.add_argument('op2', help='Type of the second operator')
parser.add_argument('-p','--projection',help='Sets a projection on an atom.',default=-1)
parser.add_argument('-p2','--projection2',help='Sets a projection on a second atom. Tries to find a relation between tensors on the first \
        atom and on the second atom.',default=-1)
parser.add_argument('-f','--findsym',help='Findsym input file',default=None)
parser.add_argument('-g','--group',help='group name',default=None)
parser.add_argument('-op3',help='third operator in the linear response formula',default=None)
parser.add_argument('-b','--basis',help='Sets a coordinate basis: abc for conventional crystallographic basis, i for the one used in input \
(default). cart for a cartesian basis in which the input basis is define. \
        abc_c for orthogonalized crystalographic basis (not tested much).',default='i')
parser.add_argument('-e','--equivalent',action='store_true',help='finds response matrices for equivalent magnetic configurations. Needs output of finddsym with\
        zero moments as an input.')
parser.add_argument('--no-rename',action='store_true')
parser.add_argument('--debug',help='Controls if debug output is printed. all means all debug output is printed, symmetrize means debug\
        output for symmetrizing, rename for renaming, equiv for finding the equivalent configurations',default='')
parser.add_argument('--latex',action='store_const',const=True,default=False,help='If set, the matrices are printed also in a latex format.')
parser.add_argument('--exp',default=-1)
parser.add_argument('--print-syms',action='store_const',const=True,default=False,help='Prints all symmetry operations.')
parser.add_argument('--transform-result',action='store_const',const=True,default=False,help='By default, the symmetry operations are \
        transformed to the correct basis. If this option is chosen, the symmetry operations are not transformed and instead the \
        result is transformed. Only works for the three operators.')
parser.add_argument('--syms',default=-1,help='Choose which symmetry operations to take, the rest is ignored. Insert symmetry operation\
         numbers separated by commas with no spaces. They are numbered as they appear in the findsym output file.\
         Also can include ranges. Example: 1-3,7,9-12')
args = parser.parse_args()

op1=args.op1 #type of the first operator
op2=args.op2 #type of the second operator

atom=int(args.projection) #number of atom on which projection is done, -1 means no projection
atom2=int(args.projection2)

basis=args.basis #a list specifying the basis to be used
equiv = args.equivalent #contains name of the nonmagnetic findsym output
#debug controls whether debug output is printed
debug = args.debug.split(',')
latex = args.latex
inp = args.findsym
no_rename = args.no_rename
exp = int(args.exp)
print_syms = args.print_syms
group = args.group
op3 = args.op3
transform_result = args.transform_result
syms_sel = args.syms

if atom2 != -1:
    if atom == -1:
        sys.exit('projection2 can be setonly if projection1 is set')

if inp and group:
    sys.exit('You cannot specify both the symmetry group and Findsym input file')

if ( atom != -1 or atom2 !=-1 ) and group:
    sys.exit('Projections not possible with group name input. Use Findsym input instead.')

if equiv and group:
    sys.exit('Equivalent configurations are not possible with group name input. Use findsym input file.')

if op3 and equiv:
    sys.exit('Equivalent configurations not implemented for three operators.')

if op3 and (exp != -1):
    sys.exit('Expansions are not implemented for three operators.')

debug_sym = False
debug_rename = False
debug_equiv = False
debug_tensor = False
debug_time = False
debug_Y = False

if 'symmetrize' in debug or 'all' in debug:
    debug_sym = True
if 'rename' in debug or 'all' in debug:
    debug_rename = True
if 'equiv' in debug or 'all' in debug:
    debug_equiv = True
if 'exp' in debug or 'all' in debug:
    debug_tensor = True
if 'time' in debug or 'all' in debug:
    debug_time = True
if 'Y' in debug:
    debug_Y = True

if inp:
    #runs findsym and reads the output
    with open(inp,'r') as f:

        fin = f.readlines()
        
        #fin_c cleans definition of axes from the findsym input
        #otherwise findsym crashes
        fin_c = []
        found = False
        i = 0
        for i in range(len(fin)):
            if 'axes:' in fin[i]:
                found = True
            if not found:
                fin_c.append(fin[i])

        try:  
            fs = subprocess.Popen([dirname+'/findsym'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out = fs.communicate(input=''.join(fin_c))[0]
            lines = out.split('\n')
        except:
            sys.exit('Error in findsym input') 

    #reads the input
    #vec_a,b,c are needed to know the basis transformation
    #syms contain the symmetries in the form that is needed by symmetr
    [vec_a,vec_b,vec_c] = read.r_basis(lines)
    syms = read.r_sym(lines)

    #construct the transformation matrix from the magnetic basis to the selected basis
    if 'abc' == basis:

        T = sympy.Matrix(sympy.Identity(3))

    if 'i' == basis:

        T = create_Tm(vec_a,vec_b,vec_c)

    if 'cart' == basis:

        T_m = create_Tm(vec_a,vec_b,vec_c)
        T_i = create_Ti(fin)

        T = T_i*T_m

    if 'custom' == basis:
        
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
            T_c[0,0] = T_c[0,0] / norm
            T_c[1,0] = T_c[1,0] / norm
            T_c[2,0] = T_c[2,0] / norm

            norm = sympy.sqrt(T_c[0,1]**2 + T_c[1,1]**2 + T_c[2,1]**2)
            T_c[0,1] = T_c[0,1] / norm
            T_c[1,1] = T_c[1,1] / norm
            T_c[2,1] = T_c[2,1] / norm

            norm = sympy.sqrt(T_c[0,2]**2 + T_c[1,2]**2 + T_c[2,2]**2)
            T_c[0,2] = T_c[0,2] / norm
            T_c[1,2] = T_c[1,2] / norm
            T_c[2,2] = T_c[2,2] / norm

        T_m = create_Tm(vec_a,vec_b,vec_c)

        T = T_c.inv()*T_i*T_m

    if 'abc_c' == basis:

        T = sympy.zeros(3)

        abc = read.r_abc(lines)

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

if group:
    atom = -1
    print group
    hex_group,syms=group_sym(group,dirname=str(dirname),debug=False)

    if 'i' == basis or 'abc' == basis:
        print 'Using the conventional coordinate system!'
        T = sympy.Matrix(sympy.Identity(3))

    if 'cart' == basis:

        print 'Using a cartesian coordinate system'
        if hex_group:
            T = sympy.zeros(3)
            T[0,0] = 1
            T[0,1] = sympy.sympify(Fraction(-0.5))
            T[0,2] = 0
            T[1,0] = 0 
            T[1,1] = sympy.sqrt(3)/2
            T[1,2] = 0
            T[2,0] = 0
            T[2,1] = 0
            T[2,2] = 1
        else:
            T = sympy.Matrix(sympy.Identity(3))

if syms_sel != -1:
    syms_sel = syms_sel.split(',')
    syms_sel2 = []
    for i in range(len(syms_sel)):
        if '-' in syms_sel[i]:
            s = syms_sel[i].split('-')
            syms_sel2 += range(int(s[0]),int(s[1])+1)
        else:
            syms_sel2.append(int(syms_sel[i]))

    syms_new = []
    for i in range(len(syms)):
        if i+1 in syms_sel2:
            syms_new.append(syms[i])

    syms = syms_new

if print_syms:
    print 'Symmetry operations:'
    print 'Format: Number, space transformation, magnetic moment transformation, time-reversal, transformation of the sublattices'
    for sym in syms:
        print sym


if exp == -1:

    #this returns the symmetrical form of spin-response tensor for atom with index atom
    #if atom is -1 no projections are done
    #operator types are given by op1 and op2
    if not op3:
        X = symmetrize_sympy.symmetr(syms,op1,op2,atom,debug_sym)
        
        if op1 == op2:
            T_m = create_Tm(vec_a,vec_b,vec_c)
            T_i = create_Ti(fin)
            T_cart = T_i*T_m

            X = same_op_sym(X,T_cart)

        if not no_rename:
            X_T = funcs.convert_X(X,T,debug=debug_rename)
        else:
            X_T = funcs.convert_X(X,T,ren=False,debug=debug_rename)

        print 'Symmetry restricted form of the tensor, even part'
        X_T[0].pprint()
        if latex:
            X_T[0].pprint(latex=latex)
        print ''

        print 'Symmetry restricted form of the tensor, odd part'
        X_T[1].pprint()
        if latex:
            X_T[1].pprint(latex=latex)
        print ''

        if atom2 != -1:

            X_T_2 = symmetrize_sympy.symmetr_AB(syms,X_T,op1,op2,atom,atom2,T=T)

            if X_T_2 == None:
                print 'no relation with atom %s found' % atom2
            else:
                print 'Symmetrized matrix in the input basis even part, atom %s' % atom2
                X_T_2[0].pprint()
                if latex:
                    X_T_2[0].pprint(latex=latex)
                print 'Symmetrized matrix in the input basis odd part, atom %s' % atom2
                X_T_2[1].pprint()
                if latex:
                    X_T_2[1].pprint(latex=latex)
                print ''

    else:
        if transform_result == False:
            X = symmetrize_she.symmetr_3op(syms,op1,op2,op3,atom,T=T)
        else:
            X = symmetrize_she.symmetr_3op(syms,op1,op2,op3,atom)
            X_T = []
            X_T.append(funcs.convert_tensor_3op(X[0],T))
            X_T.append(funcs.convert_tensor_3op(X[1],T))
            X = X_T

        spins=['x','y','z']
        print 'even part:'
        for i in range(3):
            print 'op1=',spins[i]
            sympy.pprint(X[0].reduce(0,i).mat())
            print ''

        print 'odd part:'
        for i in range(3):
            print 'op1=',spins[i]
            sympy.pprint(X[1].reduce(0,i).mat())
            print ''

if exp != -1:

    if inp:

        #findsym output for nonmagnetic structure
        lines_nm = fs_nonmag(fin_c)

        #reads the nonmagnetic findsym output
        syms_nm = read.r_sym(lines_nm) 
        [vec_a_nm,vec_b_nm,vec_c_nm] = read.r_basis(lines_nm)
        Tnm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)

        if 'i' == basis:
            T_exp = Tnm

    if group:
        print '!!!The input group must be one of the nonmagnetic point groups, otherwise the ouput will be wrong.!!!' 
        syms_nm = syms

    X = st.symmetr_ten(syms_nm,op1,op2,exp,proj=atom,T=T,debug=debug_tensor,debug_Y=debug_Y,debug_time=debug_time)
    st.print_tensor(X)
    if latex:
      st.print_tensor(X,latex=latex)
    
if equiv and ( exp != -1 ):
    print 'You cannot use --equiv and --exp together. Equivalent configurations not supported for expansions.'

if equiv and ( exp == -1):
    #outputs also the form of the tensor for all equivalent magnetic configurations

    lines_nm = fs_nonsym(fin_c)

    #reads the nonmagnetic findsym output
    syms_nm = read.r_sym(lines_nm) 
    [vec_a_nm,vec_b_nm,vec_c_nm] = read.r_basis(lines_nm)

    #transformation matrix from the magnetic basis to the input one
    T_m = create_Tm(vec_a,vec_b,vec_c)

    #transformation matrix from the nonmagnetic basis to the input one
    T_nm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)

    T_nmX = T*T_m.inv()*T_nm

    #o_nm = np.array(read.r_origin(lines_nm))

    #the shift from the magnetic to the nonmagnetic
    #shift = np.dot(np.linalg.inv(T_nm),o_nm)

    #atomic positions including magnetic moments
    #we need the moments in the correct basis, for the the length of the vectos a,b,c is needed:
    fix_m = [np.linalg.norm(vec_a),np.linalg.norm(vec_b),np.linalg.norm(vec_c)]
    fix_m = funcs.make_rational(sympy.Matrix([[fix_m[0],fix_m[1],fix_m[2]]]))
    pos = read.r_pos(lines,fix_m)
    #converted magnetic moments to the selected basis
    mag = funcs.convert_mag(pos,T)

    if debug_equiv:
        print ''
        print 'positions in the magnetic basis:'
        for p in pos:
            print p
        print ''
        print 'transformation matrix from the magnetic to the selected basis:'
        print T
        print ''
        print 'positions in the selected basis:'
        for p in mag_t:
            print p

    #this outputs all the equaivalen configurations
    #C is a conf class, it contains both the configurations and the transformed tensors
    C = find_eq.find_equiv(X_T,op1,op2,atom,syms_nm,mag,T_nmX,debug_equiv)
    print ''
    print 'Equivalent configurations:'
    C.pprint(latex=latex)
