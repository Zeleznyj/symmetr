import re
import sys

import symmetrize_sympy
import mat

import sympy
from scipy import linalg
import math

def equal_vectors(vec1,vec2,prec):
  if len(vec1) != len(vec2):
    sysexit('Comparing two vectors of unequal size:(')
  a = True
  for i in range(len(vec1)):
    if abs(vec1[i] - vec2[i]) > prec:
      a = False
  return a

def transform_position(pos,sym):

  #transforms the position
  transd = []
  for i in range(3): #loop over components of position
    op = re.sub('-','+-',sym[1][i])
    transd.append(0)
    op = re.split('\+',op) #op contains a list of all operations we have to do with component
    op = filter(None,op) #remove empty strings from the list  
    for j in range(len(op)):
      if re.match('x',op[j]):
        transd[i] = transd[i] + float(pos[0])
      if re.match('-x',op[j]):
        transd[i] = transd[i] - float(pos[0])
      if re.match('y',op[j]):
        transd[i] = transd[i] + float(pos[1])
      if re.match('-y',op[j]):
        transd[i] = transd[i] - float(pos[1])
      if re.match('z',op[j]):
        transd[i] = transd[i] + float(pos[2])
      if re.match('-z',op[j]):
        transd[i] = transd[i] - float(pos[2])
      if re.match('-?[0-9]*[./]?[0-9]+',op[j]):
        if re.match('-?[0-9.]+/[0-9.]+',op[j]):
          op_s = op[j].split('/')
          op[j] = float(op_s[0]) / float(op_s[1])
        transd[i] = transd[i] + float(op[j])

  for i in range(3):
    if math.floor(transd[i]) != 0:
      transd[i] = transd[i] - math.floor(transd[i])

  #transforms the momentum
  transd_m = []
  for i in range(3): #loop over components of magnetic moment
    op = re.sub('-','+-',sym[2][i])
    transd_m.append(0)
    op = re.split('\+',op) #op contains a list of all operations we have to do with component
    op = filter(None,op) #remove empty strings from the list  
    for j in range(len(op)):
      if re.match('mx',op[j]):
        transd_m[i] = transd_m[i] + float(pos[3])
      if re.match('-mx',op[j]):
        transd_m[i] = transd_m[i] - float(pos[3])
      if re.match('my',op[j]):
        transd_m[i] = transd_m[i] + float(pos[4])
      if re.match('-my',op[j]):
        transd_m[i] = transd_m[i] - float(pos[4])
      if re.match('mz',op[j]):
        transd_m[i] = transd_m[i] + float(pos[5])
      if re.match('-mz',op[j]):
        transd_m[i] = transd_m[i] - float(pos[5])
      if re.match('-?[0-9]*[./]?[0-9]+',op[j]):
        if re.match('-?[0-9.]+/[0-9.]+',op[j]):
          op_s = op[j].split('/')
          op[j] = float(op_s[0]) / float(op_s[1])
        transd_m[i] = transd_m[i] + float(op[j])

  return list(transd+transd_m)



data = open('sym.out')
lines = data.readlines()

#read basis transformation
for i in range(len(lines)):
  if "Vectors a,b,c:" in lines[i]:
    pos_vectors = i

vec1 = lines[pos_vectors + 1]
vec2 = lines[pos_vectors + 2]
vec3 = lines[pos_vectors + 3]

vec_a = vec1.split()
vec_b = vec2.split()
vec_c = vec3.split()

print 'using basis:'
print 'a = ', vec_a[0], '* x + ', vec_a[1], '* y + ', vec_a[2], '* z'
print 'b = ', vec_b[0], '* x + ', vec_b[1], '* y + ', vec_b[2], '* z'
print 'c = ', vec_c[0], '* x + ', vec_c[1], '* y + ', vec_c[2], '* z'
print ''

#read space group number
for i in range(len(lines)):
  if "Space Group " in lines[i]:
    group = re.findall('[0-9]+\.[0-9]+',lines[i])

tables=open('magnetic_table_bns.txt')


#read shifts of wickoff positions
#this is neccessary to obtain all positions
found = False
end = False
for line in tables:
  rx = r"BNS: +" + re.escape(group[0])
  if re.match(rx,line):
    found = True
  if found and not end:
    if re.match('Wyckoff positions:',line) or re.match('Wyckoff positions (BNS):',line):
      end = True
      shifts = re.findall('\([0-9./]+,[0-9./]+,[0-9./]+\)\+',line)

if shifts:
  for l in range(len(shifts)):
    shift_sep = re.findall('[0-9./]+',shifts[l])
    for i in range(len(shift_sep)):
      if re.match('[0-9]+/[0-9]+',shift_sep[i]):
        t = shift_sep[i].split('/')
        shift_sep[i] = float(t[0]) / float(t[1])
      else:
        shift_sep[i] = float(shift_sep[i])
    shifts[l] = shift_sep

if shifts[0] != [0,0,0]:
  sysexit('First shift is not zero')
else:
  shifts.pop(0)

#read atomic positions
#format: list containing six floats first three are the positions the other three magnetic moment
for i in range(len(lines)):
  if "Atomic positions and magnetic moments in terms of a,b,c:" in lines[i]:
    pos_pos = i

cont = 1
positions = []
i = 0
while cont == 1:
  i = i+1
  if re.match('\A +[0-9]+ ',lines[pos_pos + i]):
    positions.append(lines[pos_pos + i])
  elif re.match('-------',lines[pos_pos + i]):
    cont = 0

for i in range(len(positions)):
  positions[i] = positions[i].split()
  positions[i].pop(0)
  for l in range(len(positions[i])):
    positions[i][l] = float(positions[i][l])


#read symmetries 
for i in range(len(lines)):
  if "_space_group_symop.magn_operation_mxmymz" in lines[i]:
    pos_syms = i

cont = 1
syms = []
i = 0
while cont == 1:
  i = i+1
  if re.match('\A[0-9]+ ',lines[pos_syms + i]):
    syms.append(lines[pos_syms + i])
  else:
    cont = 0

for i in range(len(syms)):
  syms[i] = syms[i].split()
  syms[i][1] = syms[i][1].split(',')
  syms[i][2] = syms[i][2].split(',')
  syms[i].append(syms[i][1][3])
  syms[i][1] = syms[i][1][0:3]

    
print syms
for l in range(len(syms)):
  sym_trans = []
  for i in range(len(positions)):
    trans = transform_position(positions[i],syms[l])
    for j in range(len(positions)):
      if trans == positions[j]:
        sym_trans.append((i,j))
      else:
        for k in range(len(shifts)):
          pos_shift = []
          for ind in range(3):
            pos_shift.append(positions[j][ind]+shifts[k][ind])
            pos_shift[ind] = pos_shift[ind] - math.floor(pos_shift[ind])
          for ind in range(3):
            pos_shift.append(positions[j][ind+3])
          print ''
          print positions[i]
          print syms[l]
          print trans
          print pos_shift
          if equal_vectors(trans,pos_shift,0.001):
            print 'rovna se'
            sym_trans.append((i,j))
          else:
            print 'nerovna se'

  if len(sym_trans) != len(positions):
    sysexit('Wrong number of transformed atoms. Something\'s wrong')
  syms[l].append(list(sym_trans))

print ''
for pos in positions:
  print pos

print ''
for sym in syms:
  print sym


matrix = symmetrize_sympy.symmetr(syms,0)

print 'Symmetrized matrix in the abc basis intraband term:'
print matrix[0]
print 'Symmetrized matrix in the abc basis interband term:'
print matrix[1]

#transformation matrix:

T = sympy.zeros(3,3)

T[0,0] = vec_a[0]
T[1,0] = vec_b[0]
T[2,0] = vec_c[0]

T[0,1] = vec_a[1]
T[1,1] = vec_b[1]
T[2,1] = vec_c[1]

T[0,2] = vec_a[2]
T[1,2] = vec_b[2]
T[2,2] = vec_c[2]

Ti = sympy.Matrix(linalg.inv(T))

matrix_T = []
matrix_T.append(Ti*matrix[0]*T)
matrix_T.append(Ti*matrix[1]*T)

print 'Symmetrized matrix in the original basis intraband term:'
print matrix_T[0]
print 'Symmetrized matrix in the original basis interband term:'
print matrix_T[1]
