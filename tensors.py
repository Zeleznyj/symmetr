import sympy
import copy

class tensor:

  def __init__(self,kind,dim1,dim2,name='x'):

    self.inds = makeinds(dim1,dim2)
    self.dim1 = dim1
    self.dim2 = dim2
    self.name = name

    if kind == 's':

      self.v = {}
      for ind in self.inds:

        s_tot = var_name(name,ind)
        self.v[s_tot] = sympy.symbols(s_tot)

    self.t = {}
    self.x = {}
    for ind in self.inds:
      if kind == 0:
        self.t[ind] = 0
      if kind == 's':
        n_ind = var_name(name,ind)
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


def var_name(name,num):

  s_ind = ''
  for i in num:
    s_ind = s_ind + str(i)

  s_tot = '%s%s' % (name,s_ind)
  return s_tot 
  
