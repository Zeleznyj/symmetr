import re
import copy
import sys
import math

import sympy
import numpy as np

from tensors import matrix, mat2ten
from read import transform_position

#matrices 3x3 are stored in the following 1d form
#note that python starts numbering from 0, so the following represents a matrix:
# 00 01 02
# 10 11 12
# 20 21 22
#matrix = ['00','01','02','10','11','12','20','21','22']

def convert_index(i,j):
    matrix = ['00','01','02','10','11','12','20','21','22']
    n = i*3+j
    return n


# converts from 1d form index to a 3x3 matrix indeces
def inconvert_index(n):
    matrix = ['00','01','02','10','11','12','20','21','22']
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]


def convert_index_rev(i,j):
    n = 8 - (i*3+j)
    return n


# converts from 1d form index to a 3x3 matrix indeces in a reversed order
def inconvert_index_rev(n):
    matrix = list(reversed(['00','01','02','10','11','12','20','21','22']))
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]


#transforms operator component by a symmetry operation
#op_type = [operator type,component index(0,1 or 2)]
#sym is the symmetry operation
def convert_op(sym,op_type):

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

#transforms matrix of linear response
#matrix_current is the current form of the matrix
def transform_matrix(matrix_current,sym,op1,op2,l):

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


def sym_type(atom,sym):
    #returns the index of atom transformed by sym
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
    #this solves  the system, not all the variabes may be needed, but luckily sympy doesn't complain about that
    solution = sympy.solve_linear_system(Z,b0,b1,b2,b3,b4,b5,b6,b7)
    
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


def convert_X(X,T,ren=True,debug=False):
    #converts a matrix of linear response by matrix T
    X_T = []

    if not isinstance(T,matrix):
        T = mat2ten(T)

    if ren:
        X_T.append(rename(T*X[0]*T.T(),'x',debug))
        X_T.append(rename(T*X[1]*T.T(),'x',debug))
    else:
        X_T.append(T*X[0]*T.T())
        X_T.append(T*X[1]*T.T())

    return X_T


def convert_pos(poss,T,shift):
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


def sym2mat(sym):

    #returns the matrix form of the symmetry operation both for the space and the magnetic part
    #sym is symmetry in the form returned by read.py

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
    sym_m = np.zeros((3,3))
    for i in range(3):
        trans = convert_op(sym,['s',i])
        for t in trans:
            for l in range(3):
                if t[0] == l:
                    sym_m[i,l] = t[1]

    return [sym_s,sym_sv,sym_m]

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


class confs:

    #class to store information about different magnetic configurations
    #self.nconfs is the number of configurations stored
    #self.confs[i] is the configuration (stored as a dict) with number i
    #self.Xs[i] is the tensor of the configuration i

    def __init__(self):

        self.confs = {}
        self.Xs = {}
        self.nconfs = 0

    def add(self,conf,X):

        #adds a new configuration
        #conf should be a dict: {atom number: magnetic moment}
        #X is the tensor of that configuration

        self.confs[self.nconfs] = conf
        self.Xs[self.nconfs] = X
        self.nconfs += 1
        
    def is_in(self,conf):
        #tests if a configuration is in
            
        out = False
        for i in range(self.nconfs):
            out_t = True
            for j in self.confs[i]:
                if not np.allclose(self.confs[i][j],conf[j]):
                    out_t = False
            if out_t:
                out = True

        return out

    def pprint(self,m=-1,latex=False):
        #prints everything in a (somewhat) nice form
        #if m is set it prints only configuration m

        if m == -1:
            ran = range(self.nconfs)
        else:
            ran = [m]
        
        for n in ran:
            if n== 0:
                print 'starting configuration:'
            else:
                print 'configuration %s' %n
            for p in self.confs[n]:
                print 'atom %s, m = %s, %s, %s' %(p,self.confs[n][p][0],self.confs[n][p][1],self.confs[n][p][2])
            print 'even part:'
            self.Xs[n][0].pprint()
            if latex:
                self.Xs[n][0].pprint(latex=True)
            print 'odd part:\n'
            self.Xs[n][1].pprint()
            if latex:
                self.Xs[n][1].pprint(latex=True)
            print ''

    def convert(self,T,shift):
        #converts to a different coordinate system
        #seems to work, but not fully tested, not sure if the shift will work correctly
        #not actually used anywhere, but may be useful in the future

        Tt = mat2ten(T)    

        confs_t = confs()
        for n in range(self.nconfs):

            conf_t = {}
            for j in self.confs[n]:
                pos = self.confs[n][j]
                pos_t = np.dot(T,pos) + np.array(shift)
                conf_t[j] = pos_t

            X = self.Xs[n]
            X_t = []
            X_t.append(Tt*X[0]*Tt.T())
            X_t.append(Tt*X[1]*Tt.T())

            confs_t.add(conf_t,X_t)

        return confs_t


def find_equiv(X,op1,op2,atom,syms,pos,T,shift,debug=False):

    #takes a tensor and a list of nonmagnetic symmetries and find the form of the tensor for all equivalent configurations
    #X the tensor for the starting configuration
    #op1,op2 type of the operators
    #atom, the projection
    #syms are the nonmagnetic symmetry operations
    #pos are the positions and magnetic moments - ie the starting configuration
    #T is the transformation matrix from the nonmagnetic basis to the input basis
    #shift is the shift from the nonmagnetic to the input
    #will only work if the same input basis was used for the nonmagnetic and magnetic!
    
    if debug:
        print 'starting find_equiv'

    #extracts the starting configuration, only the magnetic moments are needed
    start_conf = {}
    for i in range(len(pos)):
        a = round(float(pos[i][3]),5)
        b = round(float(pos[i][4]),5)
        c = round(float(pos[i][5]),5)
        if a!=0 or b !=0 or c != 0:
            start_conf[pos[i][6]] = np.array([a,b,c])

    #creates a conf class, which stores all the configurations and adds the starting one
    C = confs()
    C.add(start_conf,X)

    if debug:
        C.pprint(0)

    for sym in syms: #a loop over all symmetries
        
        #if there is a projection take only the symmetry that keeps the atom invariant
        if atom == -1 or atom == sym_type(atom,sym):

            if debug:
                print 'taking sym: '
                print sym

            #converts the symmetry to the input basis
            sym_mag = convert_sym(sym,T,shift)

            if debug:
                print 'symmetry in the magnetic basis:'
                print sym_mag

            #transforms the starting configuration by the symmetry
            conf_t = {}
            for p in start_conf:
                pos_t = [0,0,0,start_conf[p][0],start_conf[p][1],start_conf[p][2]]
                pos_t = transform_position(pos_t,sym_mag,1.e-13)
                conf_t[sym_type(p,sym)] = np.array([round(pos_t[3],5),round(pos_t[4],5),round(pos_t[5],5)])

            if debug:
                print 'symmetry transforms starting configuration to configuration:'
                for p in conf_t:
                    print 'atom %s, m = %s, %s, %s' %(p,conf_t[p][0],conf_t[p][1],conf_t[p][2])
            
            #if the configuration has not been found before it transforms the tensor by the symmetry
            #operation and adds everything to C
            if not C.is_in(conf_t):
                if debug:
                    print 'configuration has not been found before'
                Xt = []
                for l in range(2):
                    Xt.append(transform_matrix(X[l],sym_mag,op1,op2,l))
                if debug:
                    print 'even part converted to:'
                    Xt[0].pprint()
                    Xt[1].pprint()
                C.add(conf_t,Xt)
    
    return C


#returns a symmetrical form of a response matrix for a given atom and given list of symmetries
#op1 is the type of first operator in the linear response formula, op2 is the second one
#can be set to either 's', which means spin or 'v' which means velocity
#if proj is set to -1 (default) then no projections are taken so all symmetry operations are considered
#if proj is set to atom index then only the symmetry operations that transform this atom into itself are considered
#debug is an optional parameter, if it's true, then the routine outputs lots of information
def symmetr(symmetries,op1,op2,proj=-1,debug=False):

    #this defines starting response matrix
    #we repeat it twice, once for the even part and once for the odd part
    #sympy is a symbolic toolbox
    #the matrices are symbolic matrices
    
    X1 = matrix('s',3)
    X2 = matrix('s',3)
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
            #we do everything separately for the intra and inter band terms
            #most things are the same, the only difference in the physics is that when time-reversal is present, odd part  transformation has
            #minus compared to the even part transformation
            for l in range(2):
                if debug:
                    if l == 0:
                        print ''
                        print 'Intraband term:'
                        print ''
                    if l == 1:
                        print ''
                        print 'Interband term:'
                        print ''

                #this transforms the matrix by the symmetry operation
                X_trans = transform_matrix(X[l],sym,op1,op2,l)
                
                if debug:
                    print ''
                    print 'Current form of the matrix:'
                    print ''
                    sympy.pprint(X[l].mat())
                    print ''
                    print 'Transformed matrix:'
                    print ''
                    sympy.pprint(X_trans.mat())

                #the matrix must be equal to the transformed matrix, this give us a system of 9 linear equations
                #matrix Y is a matrix that represents this system, ie the system X-X_trans = 0
                #we reverse the order of the rows - ie the first row corresponds to x[2,2] and last to x[0,0]
                # it doesn't really matter but the results are more natural this way
                Y = matrix(0,9)

                #we do a loop over all rows of the matrix Y - ie over all linear equations
                for i in range(3):
                    for j in range(3):

                        #convert_index transforms an index in a 3x3 matrix into an index in a 1x9 vector form
                        m = convert_index(i,j)

                        #a loop over all columns of matrix Y
                        for k in range(3):
                            for ll in range(3):

                                #again converts an index from 3x3 matrix form to the 1x9 vector form, but in this case in the reversed order
                                n = convert_index_rev(k,ll)

                                #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                                Y_p = X[l][i,j]-X_trans[i,j]
                                for o in range(3):
                                    for p in range(3):
                                        if o == k and p == ll:
                                            Y_p = Y_p.subs(X[l].x[o,p],1)
                                        else:
                                            Y_p = Y_p.subs(X[l].x[o,p],0)

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
                    i = 8
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
                    for ll in range(j+1,9):
                        tmp = tmp - rref[i,ll]*X[l].x[inconvert_index_rev(ll)[0],inconvert_index_rev(ll)[1]]
                    X[l] = X[l].subs(X[l].x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]],tmp)

                    if debug:
                        print 'substituting ',
                        sympy.pprint(X[l].x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]])
                        print ' for ',
                        sympy.pprint(tmp)
                        print ''


                if debug:
                    print 'Current form of the matrix:'
                    sympy.pprint(X[l].mat())
                    print ''

    if debug:
        print ''
        print 'Symmetrized matrix even part:'
        sympy.pprint(X[0].mat())
        print ''
        print 'Symmetrized matrix odd part:'
        sympy.pprint(X[1].mat())
        print ''
        print '======= End symmetrizing ======='

    if op1 == op2:
        #if operators are the same then the even part has to be symmetrical and odd part antisymetrical

        X[0] = sym_part(X[0])
        X[1] = asym_part(X[1])

    return X


def symmetr_AB(syms,X,op1,op2,atom1,atom2,T=None):
    
    #tries to transform the tensor projected on one atom to a different atom
    #syms are the symmmetry operations
    #X is the input tensor
    #op1 is the first operator
    #op2 is the second operator
    #atom1 is the atom on which X is projected
    #atom2 is the atom on which X is transformed
    #if T is set, it determines the transformation used to transform the symmetries
        #T has to be a numpy matrix

    X_trans = []

    found = False
    for sym in syms:
        #there will usually be more symmetries that transform from atom1 to atom2, we need only one, as they all
        #give the same results
        if sym_type(atom1,sym) == atom2 and not found:
            found = True
            for l in range(2):
                if T != None: #if T is set, the symmetries are transformed by it
                    sym_trans = convert_sym(sym,T,np.zeros(3))
                    X_trans.append(transform_matrix(X[l],sym_trans,op1,op2,l))
                else:
                    X_trans.append(transform_matrix(X[l],sym,op1,op2,l))

    if found:
        return X_trans
    else:
        return None


#renames the matrix so that it has the simplest form possible
#it looks for the relation between components so no information is lost
#debug is an optional parameter, if it's true, then the routine outputs lots of information
def rename(X,name,debug=False):

    #v contains the symbolic variables:
    V = matrix('s',3,name)
    #Y will be the renamed matrix
    Y = matrix(0,3)

    if debug:
        print ''
        print '======= Start renaming ======='
        print 'Input matrix:'
        sympy.pprint(X.mat())

    #a loop over all components of  the input matrix
    for i in range(3):
        for j in range(3):

            #if a components is zero we do nothing with it
            if X[i,j] == 0: 
                Y[i,j] = 0
            else:

                if debug:
                    print ''
                    print 'component: ', i, j, ' = ',
                    sympy.pprint(X[i,j])

                pos = convert_index(i,j)
                if pos == 0: #we always rename the first component of the matrix (unless it's zero)
                    Y[i,j] = V.x[i,j]
                    if debug:
                        print 'renaming to:', ' ',
                        sympy.pprint(V.x[i,j])
                # if we don't have the first component we try if the component is not a linear combination of the previous
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
                        print 'linear equations system: (last columng is the right hand side)'
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
                        print 'solution of the system:'
                        print solution

                    #if the component is not a linear combination of previous components we rename it
                    if not solution:
                        Y[i,j] = V.x[i,j]

                        if debug:
                            print 'renaming to:', ' ',
                            sympy.pprint(V.x[i,j])

                    #if it is we set it equal to the linear combination
                    else:
                        for r in range(len(solution)):
                            Y[i,j] = Y[i,j] + solution[r] * Y[inconvert_index(r)[0],inconvert_index(r)[1]]
                        if debug:
                            print 'renaming to:', ' ',
                            sympy.pprint(Y[i,j])
                    

    if debug:
        print ''
        print 'Input matrix:'
        sympy.pprint(X.mat())
        print ''
        print 'renamed matrix:'
        sympy.pprint(Y.mat())
        print ''
        print '======= end rename ======='
        print '' 

    return Y
