import re
import copy
import mat
import sympy
import sys

#matrices 3x3 are stored in the following 1d form
#note that python starts numbering from 0, so the following represents a matrix:
# 00 01 02
# 10 11 12
# 20 21 22
#matrix = ['00','01','02','10','11','12','20','21','22']

def convert_index(i,j):
  n = i*3+j
  return n

# converts from 1d form index to a 3x3 matrix indeces
def inconvert_index(n):
  matrix = ['00','01','02','10','11','12','20','21','22']
  i=int(matrix[n][0])
  j=int(matrix[n][1])
  return [i,j]


def convert_space(s):
  if re.match('x((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 0
    minus = 1
  if re.match('-x((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 0
    minus = -1
  if re.match('y((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 1
    minus = 1
  if re.match('-y((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 1
    minus = -1
  if re.match('z((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 2
    minus = 1
  if re.match('-z((\+|-)[0-9]+(\.[0-9]+)*/[0-9]+(\.[0-9]+)*)*\Z',s):
    index = 2
    minus = -1

  return [index,minus]

def convert_spin(s):
  if s == 'mx':
    index = 0
    minus = 1
  if s == '-mx':
    index = 0
    minus = -1
  if s == 'my':
    index = 1
    minus = 1
  if s == '-my':
    index = 1
    minus = -1
  if s == 'mz':
    index = 2
    minus = 1
  if s == '-mz':
    index = 2
    minus = -1
  return [index,minus]

def create_matrix_trans(sym,l):
  #creates matrix_trans, which is a matric that contains information about transformation of the matrix under symmetry sym
  #format: matrix_trans[i][j] is a list with three components, first two are indices of the transformed component
  #last on is -1 if the component is equal to minus the transformed component
  #if l==0 it means intraband term, if l==1 it means interband term

  #!!!This should be rewritten so that it can take hexagonal transformations into account

  #space transformation: 
  sym1 = sym[1]
  #-1 if the symmetry contains time reversal, 1 otherwise
  time_reversal = sym[3]
  #spin transformation
  sym2 = sym[2]

  matrix_trans = []

  for i in range(3):
    matrix_trans_temp = []
    for j in range(3):

      #index = mat.inconvert_index(n)
      v_trans = convert_space(sym1[j])
      s_trans = convert_spin(sym2[i])

      #if there is a time-reversal symmetry present, v will will have additional minus compared to space transformation
      if time_reversal == '-1':
        v_trans[1] = -1 * v_trans[1]
      
      comp_trans = [0,0,0]
      comp_trans[0] = s_trans[0]
      comp_trans[1] = v_trans[0]
      if v_trans[1] * s_trans[1] > 0:
        comp_trans[2] = 1
      if v_trans[1] * s_trans[1] < 0:
        comp_trans[2] = -1


      #if time-reversal is present, interband has additional minus
      if l == 1:
        if time_reversal == '-1':
          comp_trans[2] = -1 * comp_trans[2]

      matrix_trans_temp.append(comp_trans)

      #matrix_trans contains the transformation rules that follow from the symmetry
    matrix_trans.append(list(matrix_trans_temp))

  return matrix_trans

def update_trans_current(matrix_current,matrix_trans):
  #this is the matrix transformed under symmetry transformation given by transformation matrix matrix_trans

  trans_current = sympy.zeros(3,3)

  for i in range(3):
    for j in range(3):

      trans_current[i,j] = matrix_trans[i][j][2] * matrix_current[matrix_trans[i][j][0],matrix_trans[i][j][1]]

  return trans_current

def sym_type(atom,sym):
  #returns the index of atom transformed by sym
  a = -1
  for i in sym[4]:
    if i[0] == atom:
        a = i[1]
  if a == -1:
    sys.exit('Could not find the symmetry type.')
  else:
    return a


#returns a symmetrical form of a spin-orbit torque response matrix for a given atom and given list of symmetries
def symmetr(symmetries,atom):

  #this defines starting response matrix
  #we repeat it twice, once for intraband term and once for the interband term
  #sympy is a symbolic toolbox
  #the matrices are symbolic matrices
  Xo = sympy.Matrix(sympy.MatrixSymbol('Xo',3,3))
  Xx = sympy.Matrix(sympy.MatrixSymbol('Xx',3,3))
  X = []
  X.append(Xo)
  X.append(Xx)

  #we do a loop over all symmetry, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
  #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
  for sym in symmetries:

    
    #we only consider symmetries that keep the atom invariant
    if sym_type(atom,sym) == atom:
      #we do everything separately for the intra and inter band terms
      #most things are the same, the only difference in the physics is that when time-reversal is present, interband transformation has
      #minus compared to the intraband transformation
      for l in range(2):
        #for all of the components of the matrix we look at how they will be transformed by the symmetry operation

        matrix_trans = create_matrix_trans(sym,l)

        #the response matrix can have a form that is already constrained from a previous symmetry operation, we take this into consideration
        #and store the information in matrix_trans_current
        #for example under symmetry operation R, chi_00, may transform to chi_11, but we may already know that chi_11 must be equal to
        #-chi_00, then matrix_current[0] will be 00, matrix_current[4]=-00, matrix_trans[0]=11 and matrix_trans_current[0]=-00
        matrix_trans_current = update_trans_current(X[l],matrix_trans)


        #now we go through the response matrix and look at what must hold for it due to symmetry
        #only simple rules are implemented, but these should be sufficient if no hexagonal symmetries are present
        #if there are hexagonal symmetries, the program breaks anyway
        #we look at transformation and modify the response matrix according to the rules, if something changes we look at the transformation
        # again and we repeat until nothing changes

        changed = 1

        while changed == 1:

          changed = 0

          for i in range(3):
            for j in range(3):

              n = convert_index(i,j)

              #if boht the matrix and the transformed matrix are zero we do nothing, otherwise the code would loop forever
              if X[l][i,j] != 0 and matrix_trans_current[i,j] !=0:
                #if the component and the transformed component are opposite we set the component to zero
                if sympy.simplify(X[l][i,j]+matrix_trans_current[i,j]) == 0:
                  X[l][i,j] = 0
                  changed = 1
                #if the component is equal to another component we set it equal to that
                #however if we did that for each component, we wouldn't get any information
                #therefore if two components are equal we keep the one that has lower index in the 1d matrix form and the other one set
                #equal to the first one
                elif convert_index(matrix_trans[i][j][0],matrix_trans[i][j][1]) < n and sympy.simplify(X[l][i,j]-matrix_trans_current[i,j]) != 0:
                  changed = 1
                  X[l][i,j] = matrix_trans_current[i,j]

            #since the response matrix may have changed we also have to modify the matrix matrix_trans_current
          matrix_trans_current = update_trans_current(X[l],matrix_trans)


  return X

def rename(X,name):
  #takes a matrix and renames the components so that there is no redundancy
  #for example if there is a component equal to X_11+ X22 and a component X_11-X_22, we can rename first one to X_11 and the other to X_22
  #this is useful when we transform a matrix to a different basis and will also be useful for hexagonal transformations probably

  Y = sympy.zeros(3,3)
  Xname = sympy.Matrix(sympy.MatrixSymbol(name,3,3))

  for i in range(3):
    for j in range(3):

      if X[i,j] == 0:   #if it's zero we do nothing
        Y[i,j] = 0
      else:
        
        pos = convert_index(i,j)
        if pos == 0:
          Y[i,j] = Xname[i,j]
        else:
          changed = False
          for n in range(pos):
            [l,k] = inconvert_index(n)
            if sympy.simplify(X[l,k]-X[i,j]) == 0 and X[l,k] !=0 :
              Y[i,j] = Y[l,k]
              changed=True
            if sympy.simplify(X[l,k]+X[i,j]) == 0 and X[l,k] !=0:
              Y[i,j] = -Y[l,k]
              changed=True
          if not changed:
            Y[i,j] = Xname[i,j]

  return Y


        


      


    

  
