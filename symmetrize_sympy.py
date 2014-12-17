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


def update_trans_current(matrix_current,matrix_trans):

  trans_current = sympy.zeros(3,3)

  for i in range(3):
    for j in range(3):
      n = mat.convert_index(i,j)

      trans_current[i,j] = matrix_trans[n][2] * matrix_current[matrix_trans[n][0],matrix_trans[n][1]]

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


#returns a symmetrical form of a spin-orbit torque response matrix given list of symmetries
def symmetr(symmetries,atom):

  #this defines starting response matrix
  #we repeat it twice, once for intraband term and once for the interband term
  Xo = sympy.Matrix(sympy.MatrixSymbol('Xo',3,3))
  Xx = sympy.Matrix(sympy.MatrixSymbol('Xx',3,3))
  X = []
  X.append(Xo)
  X.append(Xx)

  #we do a loop over all symmetry, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
  #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
  for sym in symmetries:

    #space transformation: 
    sym1 = sym[1]
    #-1 if the symmetry contains time reversal, 1 otherwise
    time_reversal = sym[3]
    #spin transformation
    sym2 = sym[2]

    if sym_type(atom,sym) == atom:
      #print 'Symmetry operation: ', sym_split
      #we do everything separately for the intra and inter band terms
      #most things are the same, the only difference in the physics is that when time-reversal is present, interband transformation has
      #minus compared to the intraband transformation
      for l in range(2):
        #for all of the components of the matrix we look at how they will be transformed by the symmetry operation

        matrix_trans = []

        for n in range(9):

          index = mat.inconvert_index(n)
          v_trans = convert_space(sym1[index[1]])
          s_trans = convert_spin(sym2[index[0]])

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

          #matrix_trans is the transformed matrix - for example, if component 00 of matrix_trans is 11 it means the 00 component of the response matrix transforms to the 11 component
          #the last number is the sign
          matrix_trans.append(list(comp_trans))

        #the response matrix can have a form that is already constrained from a previous symmetry operation, we take this into consideration
        #and store the information in matrix_trans_current
        #for example under symmetry operation R, chi_00, may transform to chi_11, but we may already know that chi_11 must be equal to
        #-chi_00, then matrix_current[0] will be 00, matrix_current[4]=-00, matrix_trans[0]=11 and matrix_trans_current[0]=-00
        matrix_trans_current = update_trans_current(X[l],matrix_trans)

        #if l == 0:
        #  print 'Transformed matrix intraband term'
        #else:
        #  print 'Transformed matrix interband term'
        #print_matrix(matrix_trans)


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

              n = mat.convert_index(i,j)

              #if the transformed component is the same we do nothing
                #if they only differ in sign, the component must be zero
              if X[l][i,j] != 0 and matrix_trans_current[i,j] !=0:

                if sympy.simplify(X[l][i,j]+matrix_trans_current[i,j]) == 0:
                  X[l][i,j] = 0
                  changed = 1
                elif mat.convert_index(matrix_trans[n][0],matrix_trans[n][1]) < n and sympy.simplify(X[l][i,j]-matrix_trans_current[i,j]) != 0:
                  #if the component is equal to another component we set it equal to that
                  #however if we did that for each component, we wouldn't get any information
                  #therefore if two components are equal we keep the one that ha lower index in the 1d matrix form and the  other one set equal
                  #to the first one
                  #index_current = mat.convert_index(matrix_current[l][n][0],matrix_current[l][n][1])
                  #index_current = n
                  changed = 1
                  X[l][i,j] = matrix_trans_current[i,j]

            #since the response matrix may have changed we also have to modify the matrix matrix_trans_current
          matrix_trans_current = update_trans_current(X[l],matrix_trans)

  #sys.exit('kaslu na to')

  return X
    

  
