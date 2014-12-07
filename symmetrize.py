import re
import copy
import mat

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

  trans_current = []

  for n in range(9):
    if matrix_trans[n][2] == 1 or matrix_trans[n][2] == 0:
      trans_current.append(copy.deepcopy(matrix_current[mat.convert_index(matrix_trans[n][0],matrix_trans[n][1])]))
    elif matrix_trans[n][2] == -1:
      trans_current.append(copy.deepcopy(matrix_current[mat.convert_index(matrix_trans[n][0],matrix_trans[n][1])]))
      for l in range(len(trans_current[n])):
        trans_current[n][l][2] = -1 * trans_current[n][l][2]  

  return trans_current


#returns a symmetrical form of a spin-orbit torque response matrix given list of symmetries
def symmetr(symmetries):

  #this defines starting response matrix
  #we repeat it twice, once for intraband term and once for the interband term
  matrix_current = []
  sub1 = []
  sub2 = []
  for n in range(9):
    comp = [0,0,0]
    comp[0]=mat.inconvert_index(n)[0]
    comp[1]=mat.inconvert_index(n)[1]
    comp[2]=1
    sub1.append([list(comp)])
    sub2.append([list(comp)])

  matrix_current.append(list(sub1))
  matrix_current.append(list(sub2))


  #we do a loop over all symmetry, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
  #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
  for sym in symmetries:

    #this just splits the symmetry information into more practical form
    sym_split = sym.split()
    #this is a type of information, so far A means consider the operation and anything else means skip it
    op_type=sym_split[1]
    sym_split2 = sym_split[2].split(',')
    #space transformation: 
    sym1 = sym_split2[0:3]
    #-1 if the symmetry contains time reversal, 1 otherwise
    time_reversal = sym_split2[3]
    #spin transformation
    sym2 = sym_split[3].split(',')

    if op_type == 'A':
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
        matrix_trans_current = copy.deepcopy(update_trans_current(matrix_current[l],matrix_trans))
        
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

          for n in range(9):
            #if the transformed component is the same we do nothing
              #if they only differ in sign, the component must be zero
            if mat.compare_sign(matrix_current[l][n],matrix_trans_current[n]):
              matrix_current[l][n] = [[0,0,0]]
              changed = 1
            elif mat.convert_index(matrix_trans[n][0],matrix_trans[n][1]) < n and matrix_current[l][n] != matrix_trans_current[n]:
              #if the component is equal to another component we set it equal to that
              #however if we did that for each component, we wouldn't get any information
              #therefore if two components are equal we keep the one that ha lower index in the 1d matrix form and the  other one set equal
              #to the first one
              #index_current = mat.convert_index(matrix_current[l][n][0],matrix_current[l][n][1])
              #index_current = n
              changed = 1
              matrix_current[l][n] = copy.deepcopy(matrix_trans_current[n])

          #since the response matrix may have changed we also have to modify the matrix matrix_trans_current
          matrix_trans_current = copy.deepcopy(update_trans_current(matrix_current[l],matrix_trans))



  return matrix_current

def change_basis(matrix,basis_trans):
  
  #loop over elements of the transformed matrix
  for n in range(9):
    

  
