import re
import sys
import os
import subprocess

import symmetrize
import symmetrize_exp as st
import read
from tensors import matrix, mat2ten
import funcs
from rename import rename
from groups import group_sym
from noso import noso_syms
import input as inp

import find_eq

import sympy
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

def read_fs_inp(inp,clean=True):
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

    if clean:
        return fin_c
    else:
        return fin

def run_fs(inp):
    dirname, filename = os.path.split(os.path.abspath(__file__))
    fin_c = read_fs_inp(inp)
    try:  
        fs = subprocess.Popen([dirname+'/../findsym/findsym'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        out = fs.communicate(input=''.join(fin_c))[0]
        lines = out.split('\n')
    except:
        sys.exit('Error in findsym input') 
    return lines

def run_fs_nonmag(inp):
    
    #replaces the magnetic moments by 0
    fin_c = read_fs_inp(inp)
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
    dirname, filename = os.path.split(os.path.abspath(__file__))
    fs = subprocess.Popen([dirname+'/../findsym/findsym'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    out_nm = fs.communicate(input=''.join(fin_cnm))[0]
    lines_nm = out_nm.split('\n')

    return lines_nm

def get_syms(opt):
    if opt['inp']:
        #runs findsym and reads the output
        lines = run_fs(opt['inp'])
        syms = read.r_sym(lines)

    if opt['group']:
        atom = -1
        print opt['group']
        _,syms=group_sym(opt['group'],dirname=str(dirname),debug=False)

    #this selects some of the symmetries 
    if opt['syms_sel'] != -1:
        syms_sel = opt['syms_sel'].split(',')
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

    return syms

def get_T(opt,nonmag=False):
    """
    Returns transformation matrix from the conventional coordinate system (used by findsym)
    to the user selected basis.

    There are two different behaviors:
        nonmag=False: conventional coordinate system for the magnetic system
        nonmag=True: conventional coordinate system for the nonmagnetic system
    """
    if opt['inp']:
        #runs findsym and reads the output
        fin = read_fs_inp(opt['inp'],clean=False)
        lines = run_fs(opt['inp'])
        if nonmag:
            lines_nm = run_fs_nonmag(opt['inp'])

        #reads the input
        #vec_a,b,c are needed to know the basis transformation
        #syms contain the symmetries in the form that is needed by symmetr
        if not nonmag:
            [vec_a,vec_b,vec_c] = read.r_basis(lines)
        else:
            [vec_a,vec_b,vec_c] = read.r_basis(lines_nm)

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

            T = T_c.inv()*T_i*T_m

        if 'abc_c' == opt['basis']:

            T = sympy.zeros(3)

            if not nonmag:
                abc = read.r_abc(lines)
            else:
                abc = read.r_abc(lines_nm)

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

    if opt['group']:
        atom = -1
        print opt['group']
        hex_group,_=group_sym(opt['group'],dirname=str(dirname),debug=False)

        if 'i' == opt['basis'] or 'abc' == opt['basis']:
            print 'Using the conventional coordinate system!'
            T = sympy.Matrix(sympy.Identity(3))

        if 'cart' == opt['basis']:

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

    return T

def get_syms_nonmag(opt):
    lines = run_fs_nonmag(opt['inp'])
    syms = read.r_sym(lines)
    return syms

def get_syms_noso(opt):

    fin_c = read_fs_inp(opt['inp'])
    mags = read.r_mag_fin(fin_c)
    lines = run_fs(opt['inp'])
    lines_nm = run_fs_nonmag(opt['inp'])
    [vec_a,vec_b,vec_c] = read.r_basis(lines)
    [vec_a_nm,vec_b_nm,vec_c_nm] = read.r_basis(lines_nm)

    Tm = create_Tm(vec_a,vec_b,vec_c)
    Tnm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)
    mags_T = funcs.convert_vecs(mags,Tm.inv())
    syms_nm = get_syms_nonmag(opt)
    syms_nm_T = funcs.convert_sym_mat(syms_nm,Tm.inv()*Tnm,sym_format='findsym')
    hexag = is_hex(lines)
    syms_noso = noso_syms(syms_nm_T,mags_T,hexag,debug=opt['debug_noso'])

    if opt['syms_sel_noso'] != -1:
        syms_sel_noso = opt['syms_sel_noso'].split(',')
        syms_sel_noso2 = []
        for i in range(len(syms_sel_noso)):
            if '-' in syms_sel_noso[i]:
                s = syms_sel_noso[i].split('-')
                syms_sel_noso2 += range(int(s[0]),int(s[1])+1)
            else:
                syms_sel_noso2.append(int(syms_sel_noso[i]))

        syms_noso_new = []
        for i in range(len(syms_noso)):
            if i+1 in syms_sel_noso2:
                syms_noso_new.append(syms_noso[i])

        syms_noso = syms_noso_new

    return syms_noso

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

    #the symmetry operations given in the basis used by findsym
    syms = get_syms(opt)
    #transformation matrix from the basis used by findsym to the user defined basis
    T = get_T(opt)
    if opt['noso']:
        syms_noso = get_syms_noso(opt)

    #This is a hacky code that will need to be replaced in the future.
    if not opt['ig_op1eqop2']:
        if op1 == op2 and op3 == None:
            
            lines = run_fs(opt['inp'])
            (vec_a,vec_b,vec_c) = read.r_basis(lines)
            fin = read_fs_inp(opt['inp'])

            T_m = create_Tm(vec_a,vec_b,vec_c)
            T_i = create_Ti(fin)
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
            X_T = symmetrize.symmetrize_linres(syms_noso,op1,op2,op3=op3,proj=atom,T=T,\
                    sym_format='mat',debug=opt['debug_sym'],debug_time=opt['debug_time'],debug_Y=opt['debug_symY'])

        if not opt['ig_op1eqop2']:
            if op1 == op2 and op3 == None:

                X_T = same_op_sym(X_T,T_cart2,opt['debug_op1eqop2'],opt['debug_rename'])

    if printit:
        if op3 == None:
            print 'Symmetry restricted form of the tensor, even part'
            X_T[0].pprint()
            if opt['latex']:
                X_T[0].pprint(latex=opt['latex'])
            print ''

            print 'Symmetry restricted form of the tensor, odd part'
            X_T[1].pprint()
            if opt['latex']:
                X_T[1].pprint(latex=opt['latex'])
            print ''

        else:
            spins=['x','y','z']
            print 'even part:'
            for i in range(3):
                print 'op1=',spins[i]
                X_T[0].reduce(0,i).pprint(latex=opt['latex'])
                print ''

            print 'odd part:'
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
                print 'Symmetrized matrix in the input basis even part, atom %s' % opt['atom2']
                X_T_2[0].pprint()
                if opt['latex']:
                    X_T_2[0].pprint(latex=True)
                print 'Symmetrized matrix in the input basis odd part, atom %s' % opt['atom2']
                X_T_2[1].pprint()
                if opt['latex']:
                    X_T_2[1].pprint(latex=True)
                print ''

    #if equiv is set then we transform the tensor to all equivalent magnetic configurations
    if opt['equiv']:
        fin_c = read_fs_inp(opt['inp'])
        mags = read.r_mag_fin(fin_c)
        lines = run_fs(opt['inp'])
        lines_nm = run_fs_nonmag(opt['inp'])
        [vec_a,vec_b,vec_c] = read.r_basis(lines)
        [vec_a_nm,vec_b_nm,vec_c_nm] = read.r_basis(lines_nm)
        syms_nm = get_syms_nonmag(opt)

        #reads the magnetic moments from the input file
        mags = read.r_mag_fin(fin_c)
        #transformation matrix from the magnetic basis to the input one
        Tm = create_Tm(vec_a,vec_b,vec_c)

        #transformation matrix from the nonmagnetic basis to the input one
        Tnm = create_Tm(vec_a_nm,vec_b_nm,vec_c_nm)

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
            C.pprint(latex=opt['latex'])

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
        syms_nm = get_syms(opt)
    else:
        T = get_T(opt,nonmag=True)
        syms_nm = get_syms_nonmag(opt)

    X = st.symmetrize_exp(syms_nm,opt['op1'],opt['op2'],opt['exp'],proj=opt['atom'],T=T,\
            debug=opt['debug_sym'],debug_Y=opt['debug_symY'],debug_time=opt['debug_time'])

    if printit:
        st.print_tensor(X)
        if opt['latex']:
          st.print_tensor(X,latex=opt['latex'])

    return X
