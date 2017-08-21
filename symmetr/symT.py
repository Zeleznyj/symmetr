"""Contains functions for obtaining list of symmetry operations and transformation matrix
based on user input.
"""

import fslib
import sympy
import funcs
from sympy import sympify as spf

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

def get_syms(opt):
    if opt['inp']:
        #runs findsym and reads the output
        lines = fslib.run_fs(opt['inp'])
        syms = fslib.r_sym(lines)

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
    lines = fslib.run_fs_nonmag(opt['inp'])
    syms = fslib.r_sym(lines)
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
