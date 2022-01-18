from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
from builtins import object
import sympy
import re

class Symmetry(object):

    def __init__(self,R,has_T,Rs=None,permutations=None):
        self.R = R
        self.has_T = has_T
        if Rs is None:
            factor = 1
            if has_T:
                factor *= -1
            if self.R.det() == -1:
                factor *= -1
            self.Rs = factor * self.R
        else:
            if abs(abs(Rs.det())-1) > 1e-10:
                raise Warning('Rs must have determinant +1 or -1')
            self.Rs = Rs
        self.permutations = permutations
        self.custom_Rs = {}

    def get_R(self,op=None):
        R_op = None
        if op is None or op == 'x':
            R_op = self.R
        elif op == 's':
            R_op = self.Rs
        elif op == 'v':
            if not self.has_T:
                R_op = self.R
            else:
                R_op = -self.R
        elif op in self.custom_Rs:
            R_op = self.custom_Rs[op]

        if R_op is None:
            print(self.custom_Rs)
            print(op)
            raise Exception('wrong op')
        else:
            return R_op

    def convert(self,T,in_place=False):
        """
        Converts the symmetry operation to a different coordinate system.

        Args:
            T (matrix): Coordinate transformation matrix. If it is set, the symmetry operation will be
                        transformed by this matrix.
                        Symmetry operations are given in basis A. T transforms from A to B, ie Tx_A = x_B.

        Doesn't return anything, the symmetry itself is modified.
        """ 

        R_T = T * self.R * T.inv()
        Rs_T = T * self.Rs * T.inv()
        custom_Rs_T = {}
        for name in self.custom_Rs:
            custom_Rs_T[name] = T * self.custom_Rs[name] * T.inv()

        if in_place:
            self.R = R_T
            self.Rs = Rs_T
            for name in self.custom_Rs:
                self.custom_Rs[name] = self.custom_Rs_T[name]
        else:
            out = Symmetry(R=R_T,has_T=self.has_T,Rs=Rs_T,permutations=self.permutations)
            out.custom_Rs = custom_Rs_T
            return out

    def def_custom_R(self,name,R_c):
        self.custom_Rs[name] = R_c

    def copy(self):
        return Symmetry(self.R,self.has_T,self.Rs,self.permutations)

    def __str__(self):
        out = 'R: '
        out += self.R.__str__()  + '\n'
        out += 'Rs: '
        out += self.Rs.__str__() + '\n'
        out += 'has_T: ' + str(self.has_T) + '\n'
        out += 'permutations: ' + str(self.permutations)

        return out

    def __eq__(self,other):
        if self.R != other.R:
            return False
        if self.has_T != other.has_T:
            return False
        if self.Rs != other.Rs:
            return False
        if self.permutations is None or other.permutations is None:
            if not (self.permutations is None and other.permutations is None):
                return False
            if self.permutations != other.permutations:
                return False
        if self.custom_Rs != other.custom_Rs:
            return False
        return True

    def __mul__(self,other):
        R = self.R * other.R
        Rs = self.Rs * other.Rs
        if not self.has_T and not other.has_T:
            has_T = False
        elif self.has_T and other.has_T:
            has_T = False
        else:
            has_T = True
        if self.permutations is None and other.permutation is None:
            permutations = None
        else:
            if self.permutations is None or other.permutations is None:
                raise Exception('Cannot multiply symmetry operations when only one permutation is defined')
            permutations = {}
            for i in other.permutations:
                permutations[i] = self.permutations[other.permutations[i]]
        if self.custom_Rs.keys() != other.custom_Rs.keys():
            raise Exception('For multiplication the custom_Rs must be the same')
        custom_Rs = {}
        for name in custom_Rs:
            custom_Rs[name] = self.custom_Rs[name]*other.custom_Rs[name]
        out = Symmetry(R,has_T,Rs,permutations)
        out.custom_Rs = custom_Rs
        return out

    def pprint(self):
        sympy.pprint(self.R)
        sympy.pprint(self.Rs)
        for custom_R in self.custom_Rs:
            print(custom_R)
            sympy.pprint(self.custom_Rs[custom_R])
        print('has T:', self.has_T)
        if self.permutations is not None:
            print('permutations: ', self.permutations)

    def inv(self):
        Ri = self.R.inv()
        Rsi = self.Rs.inv()
        has_Ti = self.has_T
        if self.permutations is not None:
            permuationsi = {}
            for p in self.permutations:
                permuationsi[self.permutations[p]] = p
        else:
            permutationsi = None
        out = Symmetry(Ri,has_Ti,Rsi,permuationsi)
        for Rc in self.custom_Rs:
            out.custom_Rs[Rc] = self.custom_Rs[Rc].inv()

        return out

def findsym2sym(sym_findsym):
    R = sym2R(sym_findsym)
    Rs = sym2Rs(sym_findsym)
    if sym_findsym[3] == '-1':
        has_T = True
    elif sym_findsym[3] == '+1':
        has_T = False
    else:
        print(sym_findsym[3])
        raise Exception('Wrong findsym format')
    if len(sym_findsym) > 4:
        permutations = {}
        for (i,j) in sym_findsym[4]:
            permutations[i] = j
    else:
        permutations = None
    return Symmetry(R=R,Rs=Rs,has_T=has_T,permutations=permutations)

def matsym2sym(sym):
    R = sym[0]
    Rs = sym[2]
    if sym[3] == '-1':
        has_T = True
    elif sym[3] == '+1':
        has_T = False
    else:
        print(sym[3])
        raise Exception('Wrong findsym format')
    permutations = {}
    for (i,j) in sym[4]:
        permutations[i] = j
    return Symmetry(R=R,Rs=Rs,has_T=has_T,permutations=permutations)

def sym2R(sym):

    R = sympy.zeros(3)
    for i in range(3):
        trans = convert_op(sym,['x',i])
        for t in trans:
            for l in range(3):
                if t[0] == l:
                    R[i,l] = t[1]

    return R

def sym2Rs(sym):

    sym_s = sympy.zeros(3)
    for i in range(3):
        trans = convert_op(sym,['s',i])
        for t in trans:
            for l in range(3):
                if t[0] == l:
                    sym_s[i,l] = t[1]

    return sym_s

def sym2trans(sym):
    out = []
    for i in range(3):
        out.append(sym[1][i])
        for s in ['x','y','z']:
            out[i] = out[i].replace(s,'0')
        out[i] = sympy.sympify(out[i])
    return sympy.Matrix(out)

def create_T():

    R = sympy.diag(1,1,1)
    return Symmetry(R=R,has_T=True)

def create_I():
    R = sympy.diag(1,1,1)
    return Symmetry(R=R,has_T=False)

def create_P():
    """
    Note this doesn't contain any permutations, so it's not always correct!
    """
    R = -sympy.diag(1,1,1)
    return Symmetry(R=R,has_T=False)

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
        op = [_f for _f in op if _f] #remove empty strings from the list  
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
        op = [_f for _f in op if _f] #remove empty strings from the list  
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
        op = [_f for _f in op if _f] #remove empty strings from the list  
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
        op = [_f for _f in op if _f] #remove empty strings from the list  
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
        op = [_f for _f in op if _f] #remove empty strings from the list  
        out = []
        for j in range(len(op)):
            match = False
            if re.match('-?[0-9]*[./]?[0-9]+',op[j]):
                if re.match('-?[0-9.]+/[0-9.]+',op[j]):
                    op_s = op[j].split('/')
                    op[j] = old_div(float(op_s[0]), float(op_s[1]))
                    match = True
            if match:
                out.append(op[j])

        return out

