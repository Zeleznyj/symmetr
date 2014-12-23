import re
import sys

import symmetrize_sympy
import read

import sympy
from scipy import linalg


#opens the output from findsym
data = open('sym.out')
lines = data.readlines()

#reads the input
#vec_a,b,c are needed to know the basis transformation
#syms contain the symmetries in the form that is needed by symmetr
[vec_a,vec_b,vec_c] = read.r_basis(lines)
syms = read.r_sym(lines)


print 'using basis:'
print 'a = ', vec_a[0], '* x + ', vec_a[1], '* y + ', vec_a[2], '* z'
print 'b = ', vec_b[0], '* x + ', vec_b[1], '* y + ', vec_b[2], '* z'
print 'c = ', vec_c[0], '* x + ', vec_c[1], '* y + ', vec_c[2], '* z'
print ''


#this returns the symmetrical form of spin-response tensor for sublattice with index atom
atom=int(sys.argv[1])
matrix = symmetrize_sympy.symmetr(syms,atom)

print 'Symmetrized matrix in the abc basis intraband term:'
print matrix[0]
print 'Symmetrized matrix in the abc basis interband term:'
print matrix[1]

#transformation matrix from the original basis to the abc basis:
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

#transforms the response matrix back to the original basis
matrix_T = []
matrix_T.append(Ti*matrix[0]*T)
matrix_T.append(Ti*matrix[1]*T)

print 'Symmetrized matrix in the original basis intraband term:'
print matrix_T[0]
print 'Symmetrized matrix in the original basis interband term:'
print matrix_T[1]

mTo_rename = symmetrize_sympy.rename(matrix_T[0],'Xo')
mTx_rename = symmetrize_sympy.rename(matrix_T[1],'Xx')

print mTo_rename
print mTx_rename
