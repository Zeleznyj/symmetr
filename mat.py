import sys
import copy

#converts the matrix index to a matrix stored in the 1d form 
def convert_index(i,j):
  n = i*3+j
  return n

# converts from 1d form index to a 3x3 matrix indeces
def inconvert_index(n):
  matrix = ['00','01','02','10','11','12','20','21','22']
  i=int(matrix[n][0])
  j=int(matrix[n][1])
  return [i,j]

#this transforms a space transormation into more appropriate form
#we ignore translations as they have no effect on commutation relations

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

#used for printing matrices nicely
def print_matrix2(matrix):

  for n in range(9):

    length = len(matrix[n])

    if length == 1 and matrix[n][0][2] == 0:
      sys.stdout.write('x')
      sys.stdout.write('x')
      sys.stdout.write('x')
    else:
      for l in range(length): 
        if matrix[n][l][2] == -1:
          sys.stdout.write('-')
          sys.stdout.write(str(matrix[n][l][0]))
          sys.stdout.write(str(matrix[n][l][1]))
        elif matrix[n][l][2] == 1:
          sys.stdout.write('+')
          sys.stdout.write(str(matrix[n][l][0]))
          sys.stdout.write(str(matrix[n][l][1]))
      

    if n == 2 or n == 5 or n == 8:
      sys.stdout.write('\n')
    else:
      sys.stdout.write(', ')

  print ''


def compare_sign(comp1,comp2):
  
  a = True

  if len(comp1) != len(comp2):
    a = False
  else:
    for l in range(len(comp1)):
      if comp1[l][2] == 0 or comp2[l][2] == 0:
        a = False
      elif comp1[l][0] != comp2[l][0] or comp1[l][1] != comp2[l][1] or comp1[l][2] != -1 * comp2[l][2]:
        a =False

  return a


def compare_equal(comp1,comp2):

  a = True

  if len(comp1) != len(comp2):
    a = False
  else:
    for l in range(len(comp1)):
      if comp1[l][0] != comp2[l][0] or comp1[l][1] != comp2[l][1] or comp1[l][2] != comp2[l][2]:
        a =False

  return a


