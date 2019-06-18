# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Defines a symbolic tensor class.

All the tensors that are symmetrized use this class.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
from  six import string_types

import sympy
from sympy.core.numbers import Zero
import numpy as np
import copy
from .symmetry import create_T
import prettytable
from prettytable import PrettyTable

class tensor(object):
    """
    Creates a symbolic tensor.

    Usage:
        Use by t=tensor('s',dim1,dim2), where 's' means the tensor should be symbolic, dim2 is the number of indeces
             of the tensor, dim1 is the number of values each index can have.
        By default variables will be called x### for a different variable name, Run with tensor('s',dim1,dim2,'name').
        For nonsymbolic zero tensor, run  tensor(0,dim1,dim2).
        Print by print t
        Write out element of the tensor by t[1,2,3] or t[(1,2,3)] for example
        Assign a value to element: t[1,2,3] = 0.
        Sum or substract two tensors: t1+t2.
        Multiply by a number or symbolic variable: t*2.
        Loop over all elements, following will print all elements for example. Note that the loop is over indeces, not elements!
         for i in t:
           print t[i]
        If you need to use the symbolic variables use t.x[index], for example for x011, use t.x[0,1,1] (or t.x[(0,1,1)]).
            You can do for example t[0,0,0] = t.x[1,1,1], then t[0,0,0] will be equal to x111.
            You can also access x011 by t['x011'].
    """

    def __init__(self,kind,dim1,dim2,name='x',ind_types=None):
        """
        Creates a symbolic tensor.

        Args:
            kind : Sets what kind of tensor. 's' for symbolic, 0 for zero tensor
            dim1 (int): How many values can each index have.
            dim2 (int): How many indeces are there in the tensor.
            name (optional[str]): Name of the symbolic variables. Defaults to 'x'.
            ind_types (optinal[tuple]): Specifies whether individual indeces are covariant
                or contravariant. 1 specifies contravariant and -1 covariant. If not specified
                all the indeces are assumed contravariant.
        """

        self.inds = makeinds(dim1,dim2)
        self.dim1 = dim1
        self.dim2 = dim2
        self.name = name
        if ind_types is None:
            self.ind_types = (1,)*self.dim2
        else:
            self.ind_types = ind_types
        self.trans_type = None
        if kind == 's':
            self.v = {}
            for ind in self.inds:
                s_tot = var_name(name,ind,self.dim1)
                self.v[s_tot] = sympy.symbols(s_tot)
        # self.t contains the tensor itself, stored as a dictionary
        self.t = {}
        #self.x contains the symbolic variables in a form that can be easily retrieved
        self.x = {}
        type_found = False
        for ind in self.inds:
            if kind == 0:
                type_found = True
                self.t[ind] = 0
            if kind == 's':
                type_found = True
                n_ind = var_name(name,ind,self.dim1)
                self.t[ind] = self.v[n_ind]
                self.x[ind] = self.v[n_ind]
        if type_found == False:
            sys.exit('wrong tensor type')

        self.G = None
        self.Gi = None

    def __getitem__(self,key):
        if type(key) == tuple:
            return self.t[key]
        elif isinstance(key,int):
            return self.t[(key,)]
        elif isinstance(key,string_types):
            return self.v[key]
        raise Exception('Unexpected key type for key {} with type {}.'.format(key,type(key)))
    def __setitem__(self,key,value):
        if type(key) != tuple:
            raise TypeError
        elif len(key) != self.dim2:
            raise TypeError
        else:
            for i in key:
                if i not in list(range(self.dim1)):
                    raise TypeError
        self.t[key]=value
    def __len__(self):
        return len(self.inds)
    def __add__(self,other):
        out = self.copy0()
        for ind in self.inds:
            out[ind] = self[ind] + other[ind]
        return out
    def __sub__(self,other):
        out = self.copy0()
        for ind in self.inds:
            out[ind] = self[ind] - other[ind]
        return out
    def __mul__(self,num):
        out = self.copy0()
        for ind in self.inds:
            out[ind] = self[ind] * num
        return out
    def __div__(self,num):
        out = self.copy0()
        for ind in self.inds:
            out[ind] = old_div(self[ind], num)
        return out
    def __radd__(self,other):
        return self + other
    def __rsub__(self,other):
        return self - other
    def __rmul__(self,num):
        return self*num
    def __neg__(self):
        return self*-1
    def __str__(self):
        out = ''
        for ind in self:
            out = out + ind.__str__() + ' ' + self[ind].__str__() + '\n'
        return out
    def __iter__(self):
        return iter(self.inds)
    def __eq__(self,other):
        out = True
        if self.dim1 != other.dim1 or self.dim2 != other.dim2:
            out = False
        else:
            for ind in self.inds:
                if sympy.simplify(self[ind]-other[ind]) != 0:
                    out = False
        return out
    def __ne__(self,other):
        if self == other:
            return False
        else:
            return True
    def __type__(self):
        return 'tensor'

    def mat(self,numpy=False):
       #outputs a matrix form of itself either in sympy format (default) or numpy format 
       if self.dim2 != 2:
           raise TypeError
       if not numpy:
           out = sympy.zeros(self.dim1,self.dim1)
       if numpy:
           out = np.zeros((self.dim1,self.dim1))
       for i in range(self.dim1):
           for j in range(self.dim1):
               out[i,j] = self[i,j]
       return out

    def subs(self,old,new='notset'):
        for ind in self:
            if new != 'notset':
                self[ind] = self[ind].subs(old,new)
            else:
                self[ind] = self[ind].subs(old)
        return self
    def reduce(self,comp,value):
        out = tensor(0,self.dim1,self.dim2-1)
        for ind in out:
            ind2 = [0]*self.dim2
            for i in range(self.dim2):
                if i < comp:
                    ind2[i] = ind[i]
                if i == comp:
                    ind2[i] = value
                if i > comp:
                    ind2[i] = ind[i-1]
            ind2 = tuple(ind2)
            out[ind] = self[ind2]
        return out

    def pprint(self,print_format=None,latex=False,no_newline=False,ind_types=False,remove_zeros=False):

        if remove_zeros:
            self.remove_zeros()
        if ind_types:
            print(self.ind_types)
        if self.dim2 == 1:
            vec = [self[i] for i in range(self.dim1)]
            vec = sympy.Matrix([vec])
            if latex:
                print(sympy.latex(vec))
            else:
                sympy.pprint(vec)
        elif self.dim2 == 2:
            if latex:
                if no_newline:
                    print(sympy.latex(self.mat()), end=' ')
                else:
                    print(sympy.latex(self.mat()))
            else:
                sympy.pprint(self.mat())
        else:

            if print_format is None:
                if self.dim2 == 3:
                    print_format = 1
                else:
                    print_format = 2

            if print_format == 0:

                print(self)

            elif print_format == 1:

                if self.dim2 == 3:
                    print('X_ijk' + ' =')
                elif self.dim2 == 4:
                    print('X_ijkl' + ' =')
                else:
                    print('X_ij...pq' + ' =')

                r_inds = makeinds(self.dim1,self.dim2-2)
                for r_ind in r_inds:
                    out_str = 'X_'
                    for i in range(len(r_ind)):
                        if i == 0:
                            X = self.reduce(0,r_ind[i])
                        else:
                            X = X.reduce(0,r_ind[i])
                        #if i > 0:
                        #    out_str += ','
                        out_str += '{0}'.format(r_ind[i])
                    out_str += 'pq ='
                    print(out_str) 
                    if latex:
                        print(sympy.latex(X.mat()))
                    else:
                        X.pprint()

            elif print_format == 2:

                if latex:
                    print('WARNING: latex format not supported for print_format = 2, use print_format = 1')
                    return None

                r_inds = makeinds(self.dim1,self.dim2-2)
                r_inds2 = makeinds(self.dim1,2)
                t = PrettyTable()
                if self.dim2 == 3:
                    ind_names = ['i','j','k']
                elif self.dim2 == 4:
                    ind_names = ['i','j','k','l']
                else:
                    ind_names = ['i','j','p','q']

                t.title = '({0},{1})'.format(ind_names[-2],ind_names[-1])
                t.field_names = ([""]+[" "]+r_inds2)
                n_r_inds = len(r_inds)
                for i,r_ind in enumerate(r_inds):
                    row = [r_ind]
                    for r_ind2 in r_inds2:
                        ind = r_ind + r_ind2
                        row.append(self[ind])
                    if i == int(old_div(n_r_inds,2)):
                        if self.dim2 == 3:
                            label = "({0})".format(ind_names[0])
                        elif self.dim2 == 4:
                            label = "({0},{1})".format(ind_names[0],ind_names[1])
                        else:
                            label = "({0},{1},...)".format(ind_names[0],ind_names[1])

                    else:
                        label = ""
                    t.add_row([label]+row)
                t.hrules = prettytable.ALL
                print(t)
            else:
                raise Exception('print_format not recognized')

    def simplify(self):
        for ind in self:
            self[ind] = sympy.simplify(self[ind])

    def remove_zeros(self,num_digits=14):
        for ind in self:
            symbols = self[ind].free_symbols
            for symbol in symbols:
                if abs(self[ind].coeff(symbol)) < 1e-6:
                    self[ind] = self[ind].subs(symbol,0)

    def evalf(self,acc=15):
        out = 0

    def copy0(self):
        "returns a zero tensor with the same shape and transformation rules"
        out = tensor(0,self.dim1,self.dim2,ind_types=self.ind_types)
        if self.trans_type is not None:
            out.def_trans(ind_trans=self.ind_trans,T_comp=self.T_comp,P_trans=self.P_trans,T_trans=self.T_trans)
        if self.G is not None:
            out.def_metric(self.G)
        return out

    def copy(self):
        "returns a copy of itself"
        out = self.copy0()
        for ind in self:
            out[ind] = self[ind]
        return out

    def is_contravar(self,i):
        """Checkes if the index i is contravariant."""
        if self.ind_types[i] == 1:
            return True
        elif self.ind_types[i] == -1:
            return False
        else:
            raise Exception('Wrong ind type')

    def is_covar(self,i):
        """Checkes if the index i is covariant."""
        if self.ind_types[i] == 1:
            return False
        elif self.ind_types[i] == -1:
            return True
        else:
            raise Exception('Wrong ind type')

    def convert(self,T,in_place=True):
        """
        Return the tensor transormed by coordinate transformation matrix T.

        Args:
            T (matrix): Coordinate transformation matrix. T transforms from A to B, ie Tx_A = x_B.

        Returns:
            ten_T (tensor): the transformed tensor
        """
        
        mat1 = T
        mat2 = T.inv().T

        ten_T = self.copy0()
        for ind1 in ten_T:
            for ind2 in self:
                factor = 1
                for i in range(self.dim2):
                    if self.is_contravar(i):
                        factor *= mat1[ind1[i],ind2[i]]
                    else:
                        factor *= mat2[ind1[i],ind2[i]]
                ten_T[ind1] += factor*self[ind2]

        if in_place:
            for ind in self:
                self[ind] = ten_T[ind]

        if self.G is not None:
            G_T = T * self.G * T.T
            Gi_T = T.T.inv() * self.Gi * T.inv()
            if in_place:
                self.G = G_T 
                self.Gi = Gi_T
            else:
                ten_T.G = G_T
                ten_T.Gi = Gi_T
        
        if not in_place:
            return ten_T

    def def_trans(self,ind_trans=None,T_comp=None,P_trans=None,T_trans=None):

        inp_type_1 =  ind_trans is not None and T_comp is not None
        inp_type_2 =  P_trans is not None and T_trans is not None
        if inp_type_1 and inp_type_2:
            raise Exception('You have to specify either ind_trans and T_comp or P_trans and T_trans')
        if (not inp_type_1) and (not inp_type_2):
            raise Exception('You have to specify either ind_trans and T_comp or P_trans and T_trans')

        if inp_type_1:
            self.trans_type = 1
        if inp_type_2:
            self.trans_type = 2
        self.ind_trans = ind_trans
        self.T_comp = T_comp
        self.P_trans = P_trans
        self.T_trans = T_trans

    def transform(self,sym):

        if self.trans_type is None:
            raise Exception('To transform the tensor, the transformation rules must be defined')

        R_list = []
        for i in range(self.dim2):
            if self.trans_type == 1:
                R = sym.get_R(self.ind_trans[i])
            elif self.trans_type == 2:
                R = sym.R
            if self.is_contravar == 1:
                R_list.append(R)
            else:
                R_list.append(R.inv().T)

        ten_R = self.copy0()

        factor_ini = 1
        if self.trans_type == 1:
            if sym.has_T and self.T_comp == -1:
                factor_ini *= -1
        elif self.trans_type == 2:
            if sym.has_T and self.T_trans == -1:
                factor_ini *= -1
            if R_list[0].det() == -1 and self.P_trans == -1:
                factor_ini *= -1

        for ind1 in ten_R:
            for ind2 in self:
                factor = factor_ini
                for i in range(self.dim2):
                    factor *= R_list[i][ind1[i],ind2[i]]
                ten_R[ind1] += factor*self[ind2]

        return ten_R

    def is_even(self,sym=None):

        if self.trans_type is None:
            raise Exception('To transform the tensor, the transformation rules must be defined')

        if sym is None:
            T = create_T()
        else:
            T = sym
        X_T = self.transform(T)
        if X_T == self:
            return True
        elif X_T == -self:
            return False
        else:
            raise Exception('Wrong transformation under time-reversal')

    def def_metric(self,G):
        """
        Defines the metric for the coordinate system in which the tensor is defined

        G should be the covariant metric.
        """

        self.G = G
        self.Gi =G.inv()

    def reverse_index(self,i,in_place=True):

        if self.ind_types[i] == 1:
            if self.Gi is None:
                raise Exception("Metric tensor not defined")
            Gt = self.Gi
        if self.ind_types[i] == -1:
            if self.G is None:
                raise Exception("Metric tensor not defined")
            Gt = self.G

        out = self.copy0()

        for ind in self:
            out[ind] = 0
            for j in range(self.dim1):
                ind2 = list(ind)
                ind2[i] = j
                ind2 = tuple(ind2)
                out[ind] += Gt[ind[i],j] * self[ind2]

        if in_place:
            for ind in self:
                self[ind] = out[ind]

        ind_types = list(self.ind_types)
        ind_types[i] = -ind_types[i]
        if in_place:
            self.ind_types = tuple(ind_types)
        else:
            out.ind_types = tuple(ind_types)

        if not in_place:
            return out

    def raise_index(self,i):

        if self.ind_types[i] == 1:
            raise Exception("Index is already covariant")
        else:
            return self.reverse_index(i)

    def lower_index(self,i):

        if self.ind_types[i] == -1:
            raise Exception("Index is already contravariant")
        else:
            return self.reverse_index(i)

    def rename_vars(self,name=None,xyz=False):
        if name is None:
            name = self.name
        if xyz and self.dim1 != 3:
            raise Exception('xyz allowed only if dim1 == 3')

        new_vars = {} 

        for ind in self:
            if xyz:
                new_vars[ind] = sympy.symbols(var_name_xyz(name,ind))
            else:
                new_vars[ind] = sympy.symbols(var_name(name,ind,self.dim1))
        
        for ind in self:
            for ind2 in self:
                self[ind] = self[ind].subs(self.x[ind2],new_vars[ind2])

    def round(self,prec):
        def float2int(a):
            if abs(a) > 1e-14:
                if abs(old_div(int(a),a)-1) < 1e-14:
                    res = int(a)
                else:
                    res = a
            else:
                res = 0
            return res
        def round_expr(expr, num_digits):
            return expr.xreplace({n : float2int(round(n, num_digits)) for n in expr.atoms(sympy.Number)})
        for ind in self:
            if self[ind] != 0:
                self[ind] = round_expr(self[ind],prec)

class matrix(tensor):

     def __init__(self,kind,dim1,name='x'):
         tensor.__init__(self,kind,dim1,2,name)

     def __mul__(self,other):
         if isinstance(other, matrix):
                 out = mat2ten(self.mat()*other.mat())
         else:
             out = matrix(0,self.dim1)
             for ind in self.inds:
                 out[ind] = self[ind] * other
         return out

     def __rmul__(self,other):
         if isinstance(other, matrix):
             out = mat2ten(other.mat()*self.mat())
         else:
             out = matrix(0,self.dim1)
             for ind in self.inds:
                 out[ind] = self[ind] * other
         return out

     def __add__(self,other):
         out = matrix(0,self.dim1)
         for ind in self.inds:
             out[ind] = self[ind] + other[ind]
         return out

     def __sub__(self,other):
         out = matrix(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = self[ind] - other[ind]
         return out

     def __div__(self,num):
         out = matrix(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = old_div(self[ind], num)
         return out

     def __radd__(self,other):
         return self + other

     def __rsub__(self,other):
         return self - other

     def T(self):
         return mat2ten(self.mat().T)

    
            

def makeinds(dim1,dim2):

     for d in range(dim2):
         if d == 0:
             num = list(range(dim1))
         else:
             num1 = copy.deepcopy(num)
             num = []
             for n in num1:
                 for i in range(dim1):
                     if d == 1:
                         num.append([i]+[n])
                     else:
                         num.append([i]+n)
     for i in range(len(num)):
         if isinstance(num[i], (int,int)):
             num[i] = (num[i],)
         else:
             num[i] = tuple(num[i])

     return num


def var_name(name,num,dim):

     s_ind = ''
     for i in num:
         if dim < 11:
             s_ind = s_ind + str(i)
         else:
             s_ind = s_ind + '_' + str(i)

     if dim <11:
         s_tot = '%s%s' % (name,s_ind)
     else:
         s_tot = '%s%s' % (name,s_ind)

     return s_tot 

def var_name_xyz(name,num):
    xyz = ['x','y','z']
    s_ind = ''
    for i in num:
        s_ind = s_ind + xyz[i] 

    s_tot = name + '_' + s_ind

    return s_tot
def mat2ten(mat):
     try:
         ncols = mat.cols
         nrows = mat.rows
     except AttributeError:
         ncols = mat.shape[0]
         nrows = mat.shape[1]
     if ncols != nrows:
         raise TypeError
     else:
         out = matrix(0,ncols,2)
         for i in range(ncols):
             for j in range(ncols):
                 out[i,j] = mat[i,j]
         return out
     
