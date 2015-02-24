import sympy
import copy

class tensor:
     #this creates a symbolic tensor
     #use by t=tensor('s',dim1,dim2), where 's' means the tensor should be symbolic, dim2 is the number of indeces
     # of the tensor, dim1 is the number of values each index can have
     #by default variables will be called x### for a different variable name, run with tensor('s',dim1,dim2,'name')
     #for nonsymbolic zero tensor, run  tensor(0,dim1,dim2)
     #it's possible to:
     #print by print t
     #write out element of the tensor by t[1,2,3] or t[(1,2,3)] for example
     #assign a value to element: t[1,2,3] = 0
     #sum or substract two tensors: t1+t2
     #multiply by a number or symbolic variable: t*2
     #loop over all elements, following will print all elements for example
     # for i in t:
     #   print t[i]
     #if you need to use the symbolic variables use t.x[index], for example for x011, use t.x[0,1,1] (or t.x[(0,1,1)])
     #you can do for example t[0,0,0] = t.x[1,1,1], then t[0,0,0] will be equal to x111
     #you can also access x011 by t['x011']


     def __init__(self,kind,dim1,dim2,name='x'):

         self.inds = makeinds(dim1,dim2)
         self.dim1 = dim1
         self.dim2 = dim2
         self.name = name

         if kind == 's':

             self.v = {}
             for ind in self.inds:

                 s_tot = var_name(name,ind,self.dim1)
                 self.v[s_tot] = sympy.symbols(s_tot)

         #self.t contains the tensor itself, stored as a dictionary
         self.t = {}
         #self.x contains the symbolic variables in a form that can be easily retrieved
         self.x = {}
         for ind in self.inds:
             if kind == 0:
                 self.t[ind] = 0
             if kind == 's':
                 n_ind = var_name(name,ind,self.dim1)
                 self.t[ind] = self.v[n_ind]
                 self.x[ind] = self.v[n_ind]


     def __getitem__(self,key):
         if type(key) == tuple:
             return self.t[key]
         if type(key) == str:
             return self.v[key]

     def __setitem__(self,key,value):
         if type(key) != tuple:
             raise TypeError
         elif len(key) != self.dim2:
             raise TypeError
         else:
             for i in key:
                 if i not in range(self.dim1):
                     raise TypeError

         self.t[key]=value

     def __len__(self):
         return len(self.inds)

     def __add__(self,other):
         out = tensor(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = self[ind] + other[ind]
         return out

     def __sub__(self,other):
         out = tensor(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = self[ind] - other[ind]
         return out

     def __mul__(self,num):
         out = tensor(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = self[ind] * num
         return out

     def __div__(self,num):
         out = tensor(0,self.dim1,self.dim2)
         for ind in self.inds:
             out[ind] = self[ind] / num
         return out

     def __radd__(self,other):
         return self + other

     def __rsub__(self,other):
         return self - other

     def __rmul__(self,num):
         return self*num

     def __neg__(self):
         return t*-1

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

     def mat(self):
         if self.dim2 != 2:
             raise TypeError
         out = sympy.zeros(self.dim1,self.dim1)
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



class matrix(tensor):

     def __init__(self,kind,dim1,name='x'):
         tensor.__init__(self,kind,dim1,2,name)

     def __mul__(self,other):
         if isinstance(other, matrix):
                 out = mat2ten(self.mat()*other.mat())
         else:
             out = tensor.__mul__(self,other)
         return out

     def __rmul__(self,other):
         if isinstance(other, matrix):
             out = mat2ten(other.mat()*self.mat())
         else:
             out = tensor.__mul__(self,other)
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
             out[ind] = self[ind] / num
         return out

     def __radd__(self,other):
         return self + other

     def __rsub__(self,other):
         return self - other

     def T(self):
         return mat2ten(self.mat().T)

     def pprint(self):
         sympy.pprint(self.mat())

def makeinds(dim1,dim2):

     for d in range(dim2):
         if d == 0:
             num = range(dim1)
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
     
