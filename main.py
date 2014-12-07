import symmetrize
import mat

data = open('sym.out')
lines = data.readlines()
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

symmetries = open('sym.dat')

matrix = symmetrize.symmetr(symmetries)

print 'Symmetrized matrix in the abc basis intraband term:'
mat.print_matrix2(matrix[0])
print 'Symmetrized matrix in the abc basis interband term:'
mat.print_matrix2(matrix[1])



