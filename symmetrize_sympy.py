import re
import copy
import sympy
import sys

#matrices 3x3 are stored in the following 1d form
#note that python starts numbering from 0, so the following represents a matrix:
# 00 01 02
# 10 11 12
# 20 21 22
#matrix = ['00','01','02','10','11','12','20','21','22']

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


#transforms operator component by a symmetry operation
#op_type = [operator type,component index(0,1 or 2)]
#sym is the symmetry operation
def convert_op(sym,op_type):

  if op_type[0] == 'v':

    s = sym[1][op_type[1]]

    op = re.sub('-','+-',s)
    op = re.split('\+',op)
    op = filter(None,op) #remove empty strings from the list  
    out = []
    for j in range(len(op)):
      match = False
      if re.match('^x',op[j]):
        t = (0,1)
        match = True
      if re.match('^-x',op[j]):
        t = (0,-1)
        match = True
      if re.match('^y',op[j]):
        t = (1,1)
        match = True
      if re.match('^-y',op[j]):
        t = (1,-1)
        match = True
      if re.match('^z',op[j]):
        t = (2,1)
        match = True
      if re.match('^-z',op[j]):
        t = (2,-1)
        match = True
      #if there is a time-reversal, v has a minus compared to space transformation
      if match:
        if sym[3] == '-1':
          t = (t[0],-1*t[1])
        out.append(t)

    return out

  if op_type[0] == 's':

    s = sym[2][op_type[1]]

    op = re.sub('-','+-',s)
    op = re.split('\+',op)
    op = filter(None,op) #remove empty strings from the list  
    out = []
    for j in range(len(op)):
      if re.match('^mx',op[j]):
        t = (0,1)
        out.append(t)
      if re.match('^-mx',op[j]):
        t = (0,-1)
        out.append(t)
      if re.match('^my',op[j]):
        t = (1,1)
        out.append(t)
      if re.match('^-my',op[j]):
        t = (1,-1)
        out.append(t)
      if re.match('^mz',op[j]):
        t = (2,1)
        out.append(t)
      if re.match('^-mz',op[j]):
        t = (2,-1)
        out.append(t)

    return out


#transforms matrix of linear response
#matrix_current is the current form of the matrix
def transform_matrix(matrix_current,sym,op1,op2,l):

  trans_current = sympy.zeros(3,3)

  for i in range(3):
    for j in range(3):

      op1_trans = list(convert_op(sym,[op1,i]))
      op2_trans = list(convert_op(sym,[op2,j]))

      for v in op2_trans:
        for s in op1_trans:
          if sym[3] == '-1' and l == 1:
            trans_current[i,j] = trans_current[i,j] - v[1] * s[1] * matrix_current[s[0],v[0]]
          else:
            trans_current[i,j] = trans_current[i,j] + v[1] * s[1] * matrix_current[s[0],v[0]]

  return trans_current

def should_rename(X,X_t):
  #looks at the element and the transformation of the element and decides if the element should be renamed

  if X - X_t == 0:
    return False
  else: 
    inds = re.findall("\[([0-9]), ([0-9])\]",str(X))
    if len(inds) > 1:
      return False
    elif len(inds) == 0:
      sys.exit('error in should_rename')
    else:
      ret = True
      inds_t = re.findall("\[([0-9]), ([0-9])\]",str(X_t))
      for ind in inds_t:
        if convert_index(inds[0][0],inds[0][1]) < convert_index(ind[0],ind[1]):
          ret = False
      return ret


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


#outputs a solution to the linear equation system represented by matrix Z
#last column of Z is the right hand side
# if there is more than one solution, this chooses one solution by setting all arbitraty parameters to 0
def solve_lin(Z):

  #we need to loop over variables, this doesn't seem to be possible with sympy so this is a way around that
  b0 = sympy.Symbol('b0')
  b1 = sympy.Symbol('b1')
  b2 = sympy.Symbol('b2')
  b3 = sympy.Symbol('b3')
  b4 = sympy.Symbol('b4')
  b5 = sympy.Symbol('b5')
  b6 = sympy.Symbol('b6')
  b7 = sympy.Symbol('b7')

  def convert_b(i):
    if i == 0:
      return b0
    if i == 1:
      return b1
    if i == 2:
      return b2
    if i == 3:
      return b3
    if i == 4:
      return b4
    if i == 5:
      return b5
    if i == 6:
      return b6
    if i == 7:
      return b7

  #this solves  the system, not all the variabes may be needed, but luckily sympy doesn't complain about that
  solution = sympy.solve_linear_system(Z,b0,b1,b2,b3,b4,b5,b6,b7)
  
  #number of variables of the linear equation system
  dim = len(Z[1,:])-1
  #the list that will contain the solution
  out = []
  #if there is no solution we return None
  if not solution:
    return None
  else:
    #if there is a solution, we look at every component of the solution and set the arbitrary parameters to 0
    #the arbitrary parameters are not in the solution dict, that's how we find them
    for i in range(dim):
      x = convert_b(i)
      if x in solution:
        out_temp = solution[x]
        #we go over all variables and check if they are arbitrary
        for j in range(dim):
          y = convert_b(j)
          #if they are arbitrary we set them to zero
          if y not in solution:
            out_temp = out_temp.subs(y,0)
      else:
        out_temp = 0
      out.append(out_temp)

  return out


#returns a symmetrical form of a spin-orbit torque response matrix for a given atom and given list of symmetries
#op1 is the type of first operator in the linear response formula, op2 is the second one
#can be set to either 's', which means spin or 'v' which means velocity
#if proj is set to -1 (default) then no projections are taken so all symmetry operations are considered
#if proj is set to atom index then only the symmetry operations that transform this atom into itself are considered
#debug is an optional parameter, if it's true, then the routine outputs lots of information
def symmetr(symmetries,op1,op2,proj=-1,debug=False):

  #this defines starting response matrix
  #we repeat it twice, once for intraband term and once for the interband term
  #sympy is a symbolic toolbox
  #the matrices are symbolic matrices
  x = sympy.MatrixSymbol('x',3,3)
  X1 = sympy.Matrix(x)
  X2 = sympy.Matrix(x)
  X = []
  X.append(X1)
  X.append(X2)

  if debug:
    print ''
    print '======= Starting symmetrizing ======='

  #we do a loop over all symmetry operations, for each symmetry, we find what form the response matrix can have, when the system has this symmetry
  #for next symmetry we take the symmetrized matrix from the previous symmetry as a starting point
  for sym in symmetries:
    
    if debug:
      print 'Symmetry:' 
      print sym
      print ''
      if proj != -1:
        print 'Symmetry transforms the atom ', proj, ' into atom ', sym_type(proj,sym)
        if sym_type(proj,sym) != proj:
          print 'Skipping symmetry'
          print ''

    #if there is a projection set up we only consider symmetries that keep the atom invariant
    if proj == -1 :
      take_sym = True
    elif sym_type(proj,sym) == proj:
      take_sym = True
    else:
      take_sym = False

    if take_sym:
      #we do everything separately for the intra and inter band terms
      #most things are the same, the only difference in the physics is that when time-reversal is present, interband transformation has
      #minus compared to the intraband transformation
      for l in range(2):
        if debug:
          if l == 0:
            print ''
            print 'Intraband term:'
            print ''
          if l == 1:
            print ''
            print 'Interband term:'
            print ''

        #this transforms the matrix by the symmetry operation
        X_trans = transform_matrix(X[l],sym,op1,op2,l)
        
        if debug:
          print ''
          print 'Current form of the matrix:'
          print ''
          sympy.pprint(X[l])
          print ''
          print 'Transformed matrix:'
          print ''
          print X_trans
          sympy.pprint(X_trans)

        #the matrix must be equal to the transformed matrix, this give us a system of 9 linear equations
        #matrix Y is a matrix that represents this system, ie the system X-X_trans = 0
        #we reverse the order of the rows - ie the first row corresponds to x[2,2] and last to x[0,0]
        # it doesn't really matter but the results are more natural this way
        Y = sympy.zeros(9,9)

        #we do a loop over all rows of the matrix Y - ie over all linear equations
        for i in range(3):
          for j in range(3):

            #convert_index transforms an index in a 3x3 matrix into an index in a 1x9 vector form
            m = convert_index(i,j)

            #a loop over all columns of matrix Y
            for k in range(3):
              for ll in range(3):

                #again converts an index from 3x3 matrix form to the 1x9 vector form, but in this case in the reversed order
                n = convert_index_rev(k,ll)

                #now in the equation we substite 1 to the matrix component that correponds to the column and 0 to all others
                Y_p = X[l][i,j]-X_trans[i,j]
                for o in range(3):
                  for p in range(3):
                    if o == k and p == ll:
                      Y_p = Y_p.subs(x[o,p],1)
                    else:
                      Y_p = Y_p.subs(x[o,p],0)

                Y[m,n] = Y_p

        #this transforms the matrix into the Reduced row echelon form
        #piv are the indeces o the pivot columns
        [rref,piv] = Y.rref()        

        if debug:
          print ''
          print 'Matrix representing the linear equation system that has to be satisfied: (right hand side is zero)'
          sympy.pprint(Y)
          print ''
          print 'Reduced row echelon form and indeces of the pivot columns:'
          sympy.pprint([rref,piv])
          print ''

        #a loop over all the pivots: it's the pivots that give interesting information
        for j in list(reversed(piv)):

          
          #find the row of pivot j
          found = False
          i = 8
          while found == False:
            if rref[i,j] == 1:
              found = True
            else:
              i = i-1
          
          if debug:
            print ''
            print 'considering pivot ', i,j

          tmp = 0
          #now we just make use of the linear equation that holds for this pivot
          #keep in mind that the rows are in reversed order
          for ll in range(j+1,9):
            tmp = tmp - rref[i,ll]*x[inconvert_index_rev(ll)[0],inconvert_index_rev(ll)[1]]
          X[l] = X[l].subs(x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]],tmp)

          if debug:
            print 'substituting ',
            sympy.pprint(x[inconvert_index_rev(j)[0],inconvert_index_rev(j)[1]],)
            print ' for ',
            sympy.pprint(tmp)
            print ''


        if debug:
          print 'Current form of the matrix:'
          sympy.pprint(X[l])
          print ''

  if debug:
    print ''
    print 'Symmetrized matrix intraband term:'
    sympy.pprint(X[0])
    print ''
    print 'Symmetrized matrix interband term:'
    sympy.pprint(X[1])
    print ''
    print '======= End symmetrizing ======='

  return X


#renames the matrix so that it has the simplest form possible
#it looks for the relation between components so no information is lost
#debug is an optional parameter, if it's true, then the routine outputs lots of information
def rename(X,debug=False):

  #this cretes the neccessary symbolic variables
  x = sympy.MatrixSymbol('x',3,3)
  #Y will be the renamed matrix
  Y = sympy.zeros(3,3)

  if debug:
    print ''
    print '======= Start renaming ======='
    print 'Input matrix:'
    sympy.pprint(X)

  #a loop over all components of  the input matrix
  for i in range(3):
    for j in range(3):

      #if a components is zero we do nothing with it
      if X[i,j] == 0: 
        Y[i,j] = 0
      else:

        if debug:
          print ''
          print 'component: ', i, j, ' = ',
          sympy.pprint(X[i,j])

        pos = convert_index(i,j)
        if pos == 0: #we always rename the first component of the matrix (unless it's zero)
          Y[i,j] = x[i,j]
          if debug:
            print 'renaming to:', ' ',
            sympy.pprint(x[i,j])
        # if we don't have the first component we try if the component is not a linear combination of the previous
        #components, if it's not we give it a new name, if it is then we set it equal to that combination
        else: 
          #to find if the component is a linear combination of a previous components we set a linear equations system
          #matrix of that system will be Z, last column of Z is the right hand side
          Z = sympy.zeros(9,pos+1)
          # al loop over all components of the Z matrix
          for p in range(pos+1):
            for o in range(9):

              #we set the variable that corresponds to the column to 1 and the rest to 0
              tmp = X[inconvert_index(p)[0],inconvert_index(p)[1]]
              for m in range(9):
                if m == o:
                  tmp = tmp.subs(x[inconvert_index(m)[0],inconvert_index(m)[1]],1)
                else:
                  tmp = tmp.subs(x[inconvert_index(m)[0],inconvert_index(m)[1]],0)

              Z[o,p] = tmp
          
          if debug:
            print 'linear equations system: (last columng is the right hand side)'
            sympy.pprint(Z)

          #this will give a solution to the linear equation system
          # if there is no solution it outputs none
          #there is a sympy routine for this, but it outputs a general solution that may contain parameters
          #solve_lin routines sets all arbitrary parameters to 0, ie if there is infinitely many solutions it outputs just one
          solution = solve_lin(Z)

          if debug:
            print 'solution of the system:'
            print solution

          #if the component is not a linear combination of previous components we rename it
          if not solution:
            Y[i,j] = x[i,j]

            if debug:
              print 'renaming to:', ' ',
              sympy.pprint(x[i,j])

          #if it is we set it equal to the linear combination
          else:
            for r in range(len(solution)):
              Y[i,j] = Y[i,j] + solution[r] * Y[inconvert_index(r)[0],inconvert_index(r)[1]]
            if debug:
              print 'renaming to:', ' ',
              sympy.pprint(Y[i,j])
          

  if debug:
    print ''
    print 'Input matrix:'
    sympy.pprint(X)
    print ''
    print 'renamed matrix:'
    sympy.pprint(Y)
    print ''
    print '======= end rename ======='
    print '' 

  return Y

def rename_old(X,name):
  #takes a matrix and _renames the components so that there is no redundancy
  #for example if there is a component equal to X_11+ X22 and a component X_11-X_22, we can rename first one to X_11 and the other to X_22
  #this is useful when we transform a matrix to a different basis

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
            #if sympy.simplify(X[l,k]-X[i,j]) == 0 and X[l,k] !=0 :
            if (X[l,k]-X[i,j]) == 0 and X[l,k] !=0 :
              Y[i,j] = Y[l,k]
              changed=True
            #if sympy.simplify(X[l,k]+X[i,j]) == 0 and X[l,k] !=0:
            if (X[l,k]+X[i,j]) == 0 and X[l,k] !=0:
              Y[i,j] = -Y[l,k]
              changed=True
          if not changed:
            Y[i,j] = Xname[i,j]

  return Y


#returns a symmetrical form of a spin-orbit torque response matrix for a given atom and given list of symmetries
def symmetr_old(symmetries,atom):

  #this defines starting response matrix
  #we repeat it twice, once for intraband term and once for the interband term
  #sympy is a symbolic toolbox
  #the matrices are symbolic matrices
  X1 = sympy.Matrix(sympy.MatrixSymbol('Xo',3,3))
  X2 = sympy.Matrix(sympy.MatrixSymbol('Xx',3,3))
  X = []
  X.append(X1)
  X.append(X2)

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

        #matrix_trans = create_matrix_trans(sym,l)

        #the response matrix can have a form that is already constrained from a previous symmetry operation, we take this into consideration
        #and store the information in matrix_trans_current
        #for example under symmetry operation R, chi_00, may transform to chi_11, but we may already know that chi_11 must be equal to
        #-chi_00, then matrix_current[0] will be 00, matrix_current[4]=-00, matrix_trans[0]=11 and matrix_trans_current[0]=-00
        matrix_trans_current = transform_matrix(X[l],sym,l)


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
                #if sympy.simplify(X[l][i,j]+matrix_trans_current[i,j]) == 0:
                if (X[l][i,j]+matrix_trans_current[i,j]) == 0:
                  X[l][i,j] = 0
                  changed = 1
                #if the component is equal to another component we set it equal to that
                #however if we did that for each component, we wouldn't get any information
                #therefore if two components are equal we keep the one that has lower index in the 1d matrix form and the other one set
                #equal to the first one
                #elif convert_index(matrix_trans[i][j][0],matrix_trans[i][j][1]) < n and sympy.simplify(X[l][i,j]-matrix_trans_current[i,j]) != 0:
                #  changed = 1
                #  X[l][i,j] = matrix_trans_current[i,j]
                elif should_rename(X[l][i,j],matrix_trans_current[i,j]):
                  changed = 1
                  X[l][i,j] = matrix_trans_current[i,j]

            #since the response matrix may have changed we also have to modify the matrix matrix_trans_current
          matrix_trans_current = transform_matrix(X[l],sym,l)
          #matrix_trans_current = update_trans_current(X[l],matrix_trans)


  return X

                
                




  
  

        


      


    

  
