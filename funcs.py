import re
import copy
import sys
import math

import sympy
import numpy as np


from read import transform_position
from tensors import matrix, mat2ten,tensor
from rename import rename

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


def convert_op(sym,op_type):
    """
    Transforms operator component by a symmetry operation.

    Args:
        sym: The symmetry operation.
        op_type: Determines the operator type and operator component to be transformed.
            [operator type,component index(0,1 or 2)]

    Returns:
        out: A list of tuples. For example [(0,1),(1,-1)].
            First component means operator index. Second component means sign. The tuples are to be summed up.
            The example means: op_x-op_y
    """

    #velocity operator
    if op_type[0] == 'v':

        s = sym[1][op_type[1]]

        op = re.sub('-','+-',s)
        op = re.split('\+',op)
        op = filter(None,op) #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            match = False
            if re.match('^x',op[j]):
                t = (0,1)
                match = True
            if re.match('^-x',op[j]):
                t = (0,-1)
                match = True
            if re.match('^y',op[j]):
                t = (1,1)
                match = True
            if re.match('^-y',op[j]):
                t = (1,-1)
                match = True
            if re.match('^z',op[j]):
                t = (2,1)
                match = True
            if re.match('^-z',op[j]):
                t = (2,-1)
                match = True
            #if there is a time-reversal, v has a minus compared to space transformation
            if match:
                if sym[3] == '-1':
                    t = (t[0],-1*t[1])
                out.append(t)

        return out

    #spin operator
    if op_type[0] == 's':

        s = sym[2][op_type[1]]

        op = re.sub('-','+-',s)
        op = re.split('\+',op)
        op = filter(None,op) #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            if re.match('^mx',op[j]):
                t = (0,1)
                out.append(t)
            if re.match('^-mx',op[j]):
                t = (0,-1)
                out.append(t)
            if re.match('^my',op[j]):
                t = (1,1)
                out.append(t)
            if re.match('^-my',op[j]):
                t = (1,-1)
                out.append(t)
            if re.match('^mz',op[j]):
                t = (2,1)
                out.append(t)
            if re.match('^-mz',op[j]):
                t = (2,-1)
                out.append(t)

        return out

    #torque operator
    if op_type[0] == 't':

        s = sym[2][op_type[1]]

        op = re.sub('-','+-',s)
        op = re.split('\+',op)
        op = filter(None,op) #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            match = False
            if re.match('^mx',op[j]):
                t = (0,1)
                out.append(t)
                match = True
            if re.match('^-mx',op[j]):
                t = (0,-1)
                out.append(t)
                match = True
            if re.match('^my',op[j]):
                t = (1,1)
                out.append(t)
                match = True
            if re.match('^-my',op[j]):
                t = (1,-1)
                out.append(t)
                match = True
            if re.match('^mz',op[j]):
                t = (2,1)
                out.append(t)
                match = True
            if re.match('^-mz',op[j]):
                t = (2,-1)
                out.append(t)
                match = True
            if match:
                if sym[3] == '-1':
                    t = (t[0],-1*t[1])
                out.append(t)

        return out
    
    #transformation of a position operator
    if op_type[0] == 'x':

        s = sym[1][op_type[1]]

        op = re.sub('-','+-',s)
        op = re.split('\+',op)
        op = filter(None,op) #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            match = False
            if re.match('^x',op[j]):
                t = (0,1)
                match = True
            if re.match('^-x',op[j]):
                t = (0,-1)
                match = True
            if re.match('^y',op[j]):
                t = (1,1)
                match = True
            if re.match('^-y',op[j]):
                t = (1,-1)
                match = True
            if re.match('^z',op[j]):
                t = (2,1)
                match = True
            if re.match('^-z',op[j]):
                t = (2,-1)
                match = True
            if match:
                out.append(t)

        return out
    
    #finds a translational component of a transformation
    if op_type[0] == 'translation':
        
        s = sym[1][op_type[1]]

        op = re.sub('-','+-',s)
        op = re.split('\+',op)
        op = filter(None,op) #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            match = False
            if re.match('-?[0-9]*[./]?[0-9]+',op[j]):
                if re.match('-?[0-9.]+/[0-9.]+',op[j]):
                    op_s = op[j].split('/')
                    op[j] = float(op_s[0]) / float(op_s[1])
                    match = True
            if match:
                out.append(op[j])

        return out


def transform_matrix(matrix_current,sym,op1,op2,l,debug=False,T=None):
    """
    Transforms matrix of linear response by symmetry operation sym.

    Args:
        matrix_current (matrix, tensor class): The linear response matrix to be transformed.
        sym (Format as outputed by read.py): The symmetry transformation.
        op1 (str): First operator type.
        op2 (str): Second operator type.
        l (int): If this is 0, it is transformed as an even part. If it is 1, it is transformed as an odd part.
        T (matrix): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.

    Returns:
       trans_current: The transformed linear response matrix. 
    """
    sym_mat1 = sym2mat(sym,op1) #matrix representing the symmetry transformation of op1
    sym_mat2 = sym2mat(sym,op2) #matrix representing the symmetry transformation of op2
    #sym_mat2i = np.linalg.inv(sym_mat2) #inversion of sym_mat2
    sym_mat2i = sym_mat2.inv() #inversion of sym_mat2

    if debug:
        print ''
        print 'Symmetry matrix for operator 1:'
        sympy.pprint(sym_mat1)
        print 'Symmetry matrix for operator 2:'
        sympy.pprint(sym_mat2)

    if T != None:
        Ti = T.inv()
        #sym_mat1 =  np.dot(T,np.dot(sym_mat1,Ti)).round(8)
        #sym_mat2i =  np.dot(T,np.dot(sym_mat2i,Ti)).round(8)
        sym_mat1 = T * sym_mat1 * Ti
        sym_mat2i = T * sym_mat2i * Ti

        if debug:
            print 'Transformed symmetry matrix 1:'
            print sym_mat1
            print 'Transformed symmetry matrix 2 (inverted):'
            print sym_mat2i

    if sym[3] == '-1' and l == 1:
        trans_current = -mat2ten(sym_mat1) * matrix_current * mat2ten(sym_mat2i)
    else:
        trans_current = mat2ten(sym_mat1) * matrix_current * mat2ten(sym_mat2i)
    
    return trans_current

    #old code (remove once the new code is tested):
    if 0 == 1:
        trans_current = matrix(0,3)

        for i in range(3):
            for j in range(3):

                op1_trans = list(convert_op(sym,[op1,i]))
                op2_trans = list(convert_op(sym,[op2,j]))

                for v in op2_trans:
                    for s in op1_trans:
                        if sym[3] == '-1' and l == 1:
                            trans_current[i,j] = trans_current[i,j] - v[1] * s[1] * matrix_current[s[0],v[0]]
                        else:
                            trans_current[i,j] = trans_current[i,j] + v[1] * s[1] * matrix_current[s[0],v[0]]

        return trans_current

def transform_tensor(tensor_current,sym,op1,op2,l,T=None,debug=False):
    """
    Transforms a tensor by a symmetry operation.

    The tensor is an expansion term in the expansion of linear response tensor in powers of magnetization.
    First two indeces are transformed according to op types the rest correspons to expansion in powers of magnetization 

    Args:
        tensor_current: The tensor to be transformed.
        sym: The symmetry operation. Format as outputted by read.py.
        op1 (str): First operator type.
        op2 (str): Second operator type.
        l (int): If this is 0, it is transformed as an even part. If it is 1, it is transformed as an odd part.
        T (matrix): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.

    Returns:
        trans: The transformed tensor.

    Notes:
        This really only works for this particular tensors. However, it could easily be adapted for transformation of different
            tensor types.
    """

    mat1 = sym2mat(sym,op_type=op1)
    mat2 = sym2mat(sym,op_type=op2)
    mat3 = sym2mat(sym,op_type='s')

    if T!= None:

        if debug:
            print ''
            print 'symmetry matrix 1 (before coordinate transformation):'
            sympy.pprint(mat1)
            print ''
            print 'symmetry matrix 2 (before coordinate transformation):'
            sympy.pprint(mat2)
            print ''
            print 'symmetry matrix 3 (before coordinate transformation):'
            sympy.pprint(mat3)

        Ti = T.inv()
        mat1 = T * mat1 * Ti
        mat2 = T * mat2 * Ti
        mat3 = T * mat3 * Ti

        if debug:
            print''
            print 'The transformation matrix:'
            sympy.pprint(T)
            print ''
            print 'symmetry matrix 1 (after coordinate transformation):'
            sympy.pprint(mat1)
            print ''
            print 'symmetry matrix 2 (after coordinate transformation):'
            sympy.pprint(mat2)
            print ''
            print 'symmetry matrix 3 (after coordinate transformation):'
            sympy.pprint(mat3)


    mat2_iT = mat2.inv().T
    mat3_iT = mat3.inv().T

    if debug:
        print ''
        print 'transformation matrix 1 (first index transformation):'
        sympy.pprint(mat1)
        print '' 
        print 'transformation matrix 2 (second index transformation):'
        sympy.pprint(mat2_iT)
        print '' 
        print 'transformation matrix 3 (other indeces transformation):'
        sympy.pprint(mat3_iT)
        print '' 

    trans = tensor(0,tensor_current.dim1,tensor_current.dim2)
    
    for ind1 in trans:
        for ind2 in tensor_current:
            factor = 1
            for i in range(tensor_current.dim2):
                if i == 0:
                    factor *= mat1[ind1[i],ind2[i]]
                elif i == 1:
                    factor *= mat2_iT[ind1[i],ind2[i]]
                else:
                    factor *= mat3_iT[ind1[i],ind2[i]]
            if sym[3] == '-1' and l == 1:
                factor *= -1
            trans[ind1] += factor*tensor_current[ind2]

    return trans

def transform_tensor_3op(tensor_current,sym,op1,op2,op3,l,T=None,debug=False):
    """
    Transforms a tensor by a symmetry operation.

    The tensor describes a 3 operator linear response.

    Args:
        tensor_current: The tensor to be transformed.
        sym: The symmetry operation. Format as outputted by read.py.
        op1 (str): First operator type.
        op2 (str): Second operator type.
        op3 (str): Third operator type.
        l (int): If this is 0, it is transformed as an even part. If it is 1, it is transformed as an odd part.
        T (matrix): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.

    Returns:
        trans: The transformed tensor.

    Notes:
        This really only works for this particular tensors. However, it could easily be adapted for transformation of different
            tensor types.
    """

    mat1 = sym2mat(sym,op_type=op1)
    mat2 = sym2mat(sym,op_type=op2)
    mat3 = sym2mat(sym,op_type=op3)

    if T!= None:

        if debug:
            print ''
            print 'symmetry matrix 1 (before coordinate transformation):'
            sympy.pprint(mat1)
            print ''
            print 'symmetry matrix 2 (before coordinate transformation):'
            sympy.pprint(mat2)
            print ''
            print 'symmetry matrix 3 (before coordinate transformation):'
            sympy.pprint(mat3)

        Ti = T.inv()
        mat1 = T * mat1 * Ti
        mat2 = T * mat2 * Ti
        mat3 = T * mat3 * Ti

        if debug:
            print''
            print 'The transformation matrix:'
            sympy.pprint(T)
            print ''
            print 'symmetry matrix 1 (after coordinate transformation):'
            sympy.pprint(mat1)
            print ''
            print 'symmetry matrix 2 (after coordinate transformation):'
            sympy.pprint(mat2)
            print ''
            print 'symmetry matrix 3 (after coordinate transformation):'
            sympy.pprint(mat3)


    if debug:
        print ''
        print 'transformation matrix 1 (first index transformation):'
        sympy.pprint(mat1)
        print '' 
        print 'transformation matrix 2 (second index transformation):'
        sympy.pprint(mat2)
        print '' 
        print 'transformation matrix 3 (other indeces transformation):'
        sympy.pprint(mat3)
        print '' 

    trans = tensor(0,tensor_current.dim1,tensor_current.dim2)
    
    for ind1 in trans:
        for ind2 in tensor_current:
            factor = 1
            for i in range(3):
                if i == 0:
                    factor *= mat1[ind1[i],ind2[i]]
                elif i == 1:
                    factor *= mat2[ind1[i],ind2[i]]
                elif i == 2:
                    factor *= mat3[ind1[i],ind2[i]]
            if sym[3] == '-1' and l == 1:
                factor *= -1
            trans[ind1] += factor*tensor_current[ind2]

    return trans

def convert_tensor(ten,T):
    """
    Converts a tensor to a different coordinate system.

    Args:
        tensor: The tensor to be transformed.
        T (matrix): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.
    Returns:
        ten_T: The transformed tensor.
    """


    mat1 = T.inv()
    mat2 = T.T
    ten_T = tensor(0,ten.dim1,ten.dim2)
 
    for ind1 in ten_T:
        for ind2 in ten:
            factor = 1
            for i in range(ten.dim2):
                if i == 0:
                    factor *= mat1[ind1[i],ind2[i]]
                elif i == 1:
                    factor *= mat2[ind1[i],ind2[i]]
                else:
                    factor *= mat2[ind1[i],ind2[i]]
            ten_T[ind1] += factor*ten[ind2]
    
    return ten_T

def sym_type(atom,sym):
    """Returns number of atom to which atom transforms to under sym.

    Args:
        atom (int): Atom number to be transformed.
        sym: The symmetry operation.

    Returns:
        a (int): Transformed atom number.
    """
    a = -1
    for i in sym[4]:
        if i[0] == atom:
                a = i[1]
    if a == -1:
        sys.exit('Could not find the symmetry type. atom='+str(atom)+' symmetry='+str(sym))
    else:
        return a

def sym_part(X):
    #takes  a 3x3 matrix and returns its symmetrical part
    Y  = ( X + X.T() ) / 2
    return Y

def asym_part(X):
    Y = ( X - X.T() ) / 2
    return Y
            
def convert_X(X,T,ren=True,debug=False):
    """
    Converts a matrix of linear response by matrix T from basis A to basis B.

    Args:
        X ([matrix,matrix]): The linear response matrix to be transformed. In basis A. Matrix has to be a tensor class matrix.
            A list containg the even and the odd part of the matrix.
        T (matrix, any type should work): Transition matrix from basis A to basis B.
            This means: Tx_A = x_B, where x_A are coordinates of vector X in basis A.
        ren (Optional): Controls if the transformed tensor should have renamed components using rename.
            Defaults to True.
        debug (Optional): If set to true, debug output is printed. So far only debug output from rename is printed.

    Returns:
        X_T ([matrix,matrix]): Matrix X transformed to basis B. Matrix has to be a tensor class matrix.
    
    """
    X_T = []

    if not isinstance(T,matrix):
        T = mat2ten(T)

    #find the inverse
    Ti = mat2ten(T.mat().inv())

    if ren:
        X_T.append(rename(T*X[0]*Ti,'x',debug))
        X_T.append(rename(T*X[1]*Ti,'x',debug))
    else:
        X_T.append(T*X[0]*Ti)
        X_T.append(T*X[1]*Ti)

    return X_T

def convert_tensor_3op(ten,T):
    """
    Converts a tensor, which describes a 3 operator linear response.
    Args:
        tensor: The tensor to be transformed.
        T (matrix): Coordinate transformation matrix. If it is set, the symmetry operations will be transformed by this matrix.
            Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.
    Returns:
        ten_T: The transformed tensor.
    """


    mat1 = T
    mat2 = T.inv().T
    ten_T = tensor(0,ten.dim1,ten.dim2)
 
    for ind1 in ten_T:
        for ind2 in ten:
            factor = 1
            for i in range(ten.dim2):
                if i == 0:
                    factor *= mat1[ind1[i],ind2[i]]
                elif i == 1:
                    factor *= mat1[ind1[i],ind2[i]]
                else:
                    factor *= mat2[ind1[i],ind2[i]]
            ten_T[ind1] += factor*ten[ind2]
    
    return ten_T

def convert_pos(poss,T,shift=np.array([0,0,0])):
    #converts a positions by a matrix T and shifts them by shift

    poss_T = []
    for pos in poss:
        pos_s = np.array(pos[0:3])
        pos_m = np.array(pos[3:6])

        pos_s_t = np.dot(T,pos_s) + shift
        #shift to first unit cell
        for i in range(3):
            pos_s_t[i] = pos_s_t[i] - math.floor(pos_s_t[i])
            if np.allclose(pos_s_t[i],1):
                pos_s_t[i] = 0
            if np.allclose(pos_s_t[i],0):
                pos_s_t[i] = 0


        pos_m_t = np.dot(T,pos_m)

        pos_t = list(pos_s_t)+list(pos_m_t)
        pos_t.append(pos[6])

        poss_T.append(pos_t)

    return poss_T

def convert_mag(mag_conf,T):

    mag_conf_T = []
    for i in range(len(mag_conf)): 
        mom = sympy.Matrix([[mag_conf[i][3],mag_conf[i][4],mag_conf[i][5]]])
        mom_T = T*mom.T 
        mom_T = mom_T.T
        mom_T = sympy.Matrix([[mom_T[0],mom_T[1],mom_T[2],mag_conf[i][6]]])
        mag_conf_T.append(mom_T)

    return mag_conf_T



def sym2mat(sym,op_type=None):
    """
    Returns the matrix form of the symmetry operation.

    Args:
        sym (As returned by read.py): The symmetry operation to be transformed.
        op_type (Optional): The operator type. Defaults to None.
    
    Returns:
        If op_type is not specified it returns [matrix of space transformation,translation,
            matrix of magnetic moment transformation]
        If the op_type is specified,  then matrix for that operator type is returned.

    Notes:
        I'm not 100% sure that the translation is actually correct!!! It's not needed anywhere so not tested properly.
    """

    if not op_type:

        #space part
        sym_s = np.zeros((3,3))
        for i in range(3):
            trans = convert_op(sym,['x',i])
            for t in trans:
                for l in range(3):
                    if t[0] == l:
                        sym_s[i,l] = t[1]

        #translation
        sym_sv = np.zeros((3))
        for i in range(3):
            trans = convert_op(sym,['translation',i])
            for t in trans:
                sym_sv[i] = sym_sv[i] + t

        #magnetic part
        sym_m = np.zeros(3)
        for i in range(3):
            trans = convert_op(sym,['s',i])
            for t in trans:
                for l in range(3):
                    if t[0] == l:
                        sym_m[i,l] = t[1]

        return [sym_s,sym_sv,sym_m]

    else:

        sym_t = sympy.zeros(3)
        for i in range(3):
            trans = convert_op(sym,[op_type,i])
            for t in trans:
                for l in range(3):
                    if t[0] == l:
                        sym_t[i,l] = t[1]

        return make_rational(sym_t)


def mat2sym(msym):

    #converts the matrix represenation of the symmetry operation back to the original form
    #meaning to the form which is returned by read.py
    #handle some common fractions nicely (1/2,1/3,..), but in general will return a float so the output may
    #not look so nice in general

    sym = []
    sym.append('')

    prec = 1.e-14

    sym_s = msym[0]
    sym_sv = msym[1]
    sym_m = msym[2]

    #space part  
    sym1 = []
    for i in range(3):
        t = ''
        for j in range(3):
            if j == 0:
                symb = 'x'
            if j == 1:
                symb = 'y'
            if j == 2:
                symb = 'z'
            if abs(sym_s[i,j]-1) < prec:
                if t == '':
                    t = t + symb
                else:
                    t = t + '+' + symb
            elif abs(sym_s[i,j]+1) < prec:
                t = t + '-'+symb
            elif abs(sym_s[i,j]) > prec:
                t = t + str(sym_s[i,j])+symb
        sym1.append(t)

    #the translation
    for i in range(3):
        if sym_sv[i] > prec:
            if abs(sym_sv[i]-0.5) < prec:
                sym1[i] = sym1[i] + '+1/2'
            elif abs(sym_sv[i]-1./3.) < prec:
                sym1[i] = sym1[i] + '+1/3'
            elif abs(sym_sv[i]-2./3.) < prec:
                sym1[i] = sym1[i] + '+2/3'
            else:
                sym1[i] = sym1[i] + '+' + str(sym_sv[i])
        elif sym_sv[i] < -prec:
            if abs(sym_sv[i]+0.5) < prec:
                sym1[i] = sym1[i] + '-1/2'
            elif abs(sym_sv[i]+1./3.) < prec:
                sym1[i] = sym1[i] + '-1/3'
            elif abs(sym_sv[i]+2./3.) < prec:
                sym1[i] = sym1[i] + '-2/3'
            else:
                sym1[i] = sym1[i] + str(sym_sv[i])

    #magnetic part  
    sym2 = []
    for i in range(3):
        t = ''
        for j in range(3):
            if j == 0:
                symb = 'mx'
            if j == 1:
                symb = 'my'
            if j == 2:
                symb = 'mz'
            if abs(sym_m[i,j]-1) < prec:
                if t == '':
                    t = t + symb
                else:
                    t = t + '+' + symb
            elif abs(sym_m[i,j]+1) < prec:
                t = t + '-'+symb
            elif abs(sym_m[i,j]) > prec:
                t = t + str(sym_m[i,j])+symb
        sym2.append(t)

    sym.append(sym1)
    sym.append(sym2)

    return sym

def convert_sym(sym,T,shift,debug=False):

    #takes a symmetry operation and transforms it to a different coordinate system, represented by T and shift
    #seems to work fine, but the shifts may not be correct!!!

    prec = 1.e-13

    Ti = np.linalg.inv(T)

    sym_mat = sym2mat(sym)

    trans1 = np.dot(T,np.dot(sym_mat[0],Ti))
    trans2 = np.dot(T,np.dot(sym_mat[2],Ti))

    trans_v1 = np.dot(T,np.array(sym_mat[1])) 
    trans_v2 = np.dot(sym_mat[0],shift)-shift
    trans_v2 = np.dot(T,trans_v2)
    trans_v = trans_v1 + trans_v2

    for i in range(3):
        if math.floor(trans_v[i] + prec) > 0:
            trans_v[i] -= math.floor(trans_v[i])
        if math.floor(trans_v[i]-prec) < -1:
            trans_v[i] -= math.floor(trans_v[i])

    sym_trans = mat2sym([trans1,trans_v,trans2])
    sym_trans[0] = sym[0]
    sym_trans.append(sym[3])
    sym_trans.append(sym[4])

    return sym_trans

def convert_sym_mat(sym,T):
    """
    Converts symmetry operation represented by a matrix to a different coordinate system.

    Args:
        sym: Symmetry operation represented by a matrix. In basis A.
            A sympy matrix.
        T: Matrix repr. the transformation from basis A to basis B.
            That is a matrix such that Tx_A = x_B.
            A sympy matrix.

    Returns:
        sym_B: the symmetry operation in matrix form in basis B. A sympy matrix.
    """

    sym_B = T * sym * T.inv()

    return sym_B


