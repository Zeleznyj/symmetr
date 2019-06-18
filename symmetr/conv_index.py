from __future__ import division
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from past.utils import old_div
def convert_index(i,j):
    matrix = ['00','01','02','10','11','12','20','21','22']
    n = i*3+j
    return n

def convert_index_3(i,j,k):
    n = i*9 + j*3+k
    return n


# converts from 1d form index to a 3x3 matrix indeces
def inconvert_index(n):
    matrix = ['00','01','02','10','11','12','20','21','22']
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]

def inconvert_index_3(n):
    i = old_div(n,9)
    j = old_div((n-i*9),3)
    k = n-i*9-j*3
    return i,j,k

def convert_index_rev(i,j):
    n = 8 - (i*3+j)
    return n

# converts from 1d form index to a 3x3 matrix indeces in a reversed order
def inconvert_index_rev(n):
    matrix = list(reversed(['00','01','02','10','11','12','20','21','22']))
    i=int(matrix[n][0])
    j=int(matrix[n][1])
    return [i,j]

def convert_index_rev_3(i,j,k):
    n = 26 - convert_index_3(i,j,k)
    return n

def inconvert_index_rev_3(n):
    n = 26 - n
    return inconvert_index_3(n)

