import sys

data = open('sym.dat')

matrix = ['00','01','02','10','11','12','20','21','22']

def convert_index(i,j):
  n = i*3+j
  return n

def inconvert_index(n):
  i=int(matrix[n][0])
  j=int(matrix[n][1])
  return [i,j]

def convert_space(s):
  if s == 'x':
    index = 0
    minus = 1
  if s == '-x':
    index = 0
    minus = -1
  if s == 'y':
    index = 1
    minus = 1
  if s == '-y':
    index = 1
    minus = -1
  if s == 'z':
    index = 2
    minus = 1
  if s == '-z':
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

def print_matrix(matrix):
  for n in range(9):
    if matrix[n][2] == -1:
      sys.stdout.write('-')
      sys.stdout.write(str(matrix[n][0]))
      sys.stdout.write(str(matrix[n][1])+' ')
    elif matrix[n][2] == 1:
      sys.stdout.write(' ')
      sys.stdout.write(str(matrix[n][0]))
      sys.stdout.write(str(matrix[n][1])+' ')
    elif matrix[n][2] == 0:
      sys.stdout.write('x')
      sys.stdout.write('x')
      sys.stdout.write('x'+' ')
    if n == 2 or n == 5 or n == 8:
      sys.stdout.write('\n')
  print ''

matrix_current = []
sub1 = []
sub2 = []
for n in range(9):
  comp = [0,0,0]
  comp[0]=inconvert_index(n)[0]
  comp[1]=inconvert_index(n)[1]
  comp[2]=1
  sub1.append(list(comp))
  sub2.append(list(comp))

matrix_current.append(list(sub1))
matrix_current.append(list(sub2))

for sym in data:

  sym_split = sym.split()
  op_type=sym_split[1]
  sym_split2 = sym_split[2].split(',')
  #space transformation
  sym1 = sym_split2[0:3]
  time_reversal = sym_split2[3]
  #spin transformation
  sym2 = sym_split[3].split(',')

  if op_type == 'A':
    print 'Symmetry operation: ', sym_split
    #we do everything separately for the intra and inter band terms
    for l in range(2):
      #for all of the components of the matrix we look at how they will be transformed by the symmetry operation
      matrix_trans = []
      for n in range(9):

        index = inconvert_index(n)
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

        #matrix_trans is the transformed matrix - for example, if component 00 of matrix_trans is 11 it means the 00 component of the response
        #matrix transforms to the 11 component
        #the last number is the sign
        if l == 1:
          if time_reversal == '-1':
            comp_trans[2] = -1 * comp_trans[2]
        matrix_trans.append(list(comp_trans))

      matrix_trans_current = []
      for n in range(9):
        if matrix_trans[n][2] == 1 or matrix_trans[n][2] == 0:
          matrix_trans_current.append(list(matrix_current[l][convert_index(matrix_trans[n][0],matrix_trans[n][1])]))
        elif matrix_trans[n][2] == -1:
          matrix_trans_current.append(list(matrix_current[l][convert_index(matrix_trans[n][0],matrix_trans[n][1])]))
          matrix_trans_current[n][2] = -1 * matrix_trans_current[n][2]  

      if l == 0:
        print 'Transformed matrix intraband term'
      else:
        print 'Transformed matrix interband term'
      print_matrix(matrix_trans)
      #now we go through the matrix and look at what must hold for it due to symmetry
      changed = 1
      while changed == 1:
        changed = 0
        for n in range(9):
          #if the transofmed mcomponent is the same we do nothing
          if matrix_current[l][n] != matrix_trans_current[n]:
            #if they only differ in sign, the component must be zero
            if matrix_current[l][n][0] == matrix_trans_current[n][0] and matrix_current[l][n][1] == matrix_trans_current[n][1]:
              matrix_current[l][n][2] = 0
              changed = 1
            else:
              index_current = convert_index(matrix_current[l][n][0],matrix_current[l][n][1])
              index_trans = convert_index(matrix_trans_current[n][0],matrix_trans_current[n][1])
              if index_current > index_trans:
                changed = 1
                matrix_current[l][n] = list(matrix_trans_current[n])

        matrix_trans_current = []
        for i in range(9):
          if matrix_trans[i][2] == 1 or matrix_trans[i][2] == 0:
            matrix_trans_current.append(list(matrix_current[l][convert_index(matrix_trans[i][0],matrix_trans[i][1])]))
          elif matrix_trans[i][2] == -1:
            matrix_trans_current.append(list(matrix_current[l][convert_index(matrix_trans[i][0],matrix_trans[i][1])]))
            matrix_trans_current[i][2] = -1 * matrix_trans_current[i][2] 

    print 'Symmetrized matrix intraband term:'
    print_matrix(matrix_current[0])
    print 'Symmetrized matrix interband term:'
    print_matrix(matrix_current[1])
  










      




