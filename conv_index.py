def convert_index(i,j):
    matrix = ['00','01','02','10','11','12','20','21','22']
    n = i*3+j
    return n


# converts from 1d form index to a 3x3 matrix indeces
def inconvert_index(n):
    matrix = ['00','01','02','10','11','12','20','21','22']
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]


def convert_index_rev(i,j):
    n = 8 - (i*3+j)
    return n


# converts from 1d form index to a 3x3 matrix indeces in a reversed order
def inconvert_index_rev(n):
    matrix = list(reversed(['00','01','02','10','11','12','20','21','22']))
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]
