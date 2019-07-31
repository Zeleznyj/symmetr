# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""Contains functions for working with findsym

Includes functions for running findsym, reading findsym input, and extracting
data from findsym output.
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from past.utils import old_div
import re
import math
import sys
import os
import subprocess

from fractions import Fraction
import sympy

def get_newfilename(path,filename):
    if not os.path.exists(path+'/'+filename):
        return filename
    else:
        new_filename = filename + '1'
        i = 1
        while os.path.exists(path+'/'+new_filename):
            i += 1
            new_filename = filename + '{}'.format(i)
    return new_filename

def read_fs_inp(inp,clean=True):
    with open(inp,'r') as f:

        fin = f.readlines()
        
        #fin_c cleans definition of axes from the findsym input
        #otherwise findsym crashes
        fin_c = []
        found = False
        i = 0
        for i in range(len(fin)):
            if 'axes:' in fin[i]:
                found = True
            if not found:
                fin_c.append(fin[i])

    if clean:
        return fin_c
    else:
        return fin

def run_fs_fin(fin):
    dirname, filename = os.path.split(os.path.abspath(__file__))
    my_env = os.environ.copy()
    isodata_path = dirname + '/findsym/'

    #findsym fails if the ISODATA variable is too long.
    #This is a simple fix, which tries to create a symbolic link in the current working directory
    try:
        cwd = os.getcwd()
        findsym_name = get_newfilename(cwd,'findsym')
        os.symlink(isodata_path,cwd+'/'+findsym_name)
        my_env["ISODATA"] = findsym_name+'/'
        link_created = True
    except:
        print('Cannot create a subdirectory in a current working directory. \
            Findym may fail if the path to findsym is too long.')
        my_env["ISODATA"] = isodata_path
        link_created = False

    try:  
        fs = subprocess.Popen([dirname+'/findsym/findsym'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,env=my_env)
        out = fs.communicate(input=''.join(fin).encode())[0]
        lines = out.decode().split('\n')
        if link_created:
            try:
                os.unlink(findsym_name)
            except:
                print('Couldn not remove symlink {}'.format(cwd+'/'+findsym_name))
    except:
        if link_created:
            try:
                os.unlink(findsym_name)
            except:
                print('Couldn not remove symlink {}'.format(cwd+'/'+findsym_name))
        raise Exception('Error in findsym input') 

    return lines

def make_fsinp_nonmag(fin_c):
    start = False
    fin_cnm = []
    for i in range(len(fin_c)):
        if 'magnetic' in fin_c[i]:
            start = True
        if start:
            fin_cnm.append(re.sub(r'([0-9\.\-]+ +[0-9\.\-]+ +[0-9\.\-]+).+',r'\1 0 0 0',fin_c[i],count=1))
        else:
            fin_cnm.append(fin_c[i])
    return fin_cnm

def run_fs(inp):
    fin_c = read_fs_inp(inp)
    return run_fs_fin(fin_c)

def run_fs_nonmag(inp):
    
    fin_c = read_fs_inp(inp)
    #replaces the magnetic moments by 0
    fin_cnm = make_fsinp_nonmag(fin_c)
    
    return run_fs_fin(fin_cnm)

#tests if two vectors are equal with precision prec
def equal_vectors(vec1,vec2,prec):
     if len(vec1) != len(vec2):
         sys.exit('Comparing two vectors of unequal size:(')
     a = True
     for i in range(len(vec1)):
         if abs(vec1[i] - vec2[i]) > prec:
             a = False
     return a

#takes  a position (which includes a magnetic moment,) and symmetry and trasforms the position under the symmetry operation
#it transforms the magnetic moment too
#transformation of the magnetic moment is in principle not neccessary
def transform_position(pos,sym,prec):

     #transforms the position
     transd = []
     for i in range(3): #loop over components of position
         op = re.sub('-','+-',sym[1][i])
         transd.append(0)
         op = re.split('\+',op) #op contains a list of all operations we have to do with component
         op = [_f for _f in op if _f] #remove empty strings from the list  
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
                     op[j] = old_div(float(op_s[0]), float(op_s[1]))
                 transd[i] = transd[i] + float(op[j])

     #we want all positions to be in the first unit cell, ie they must lie between 0 and 1
     #if they are not we tranform it there
     for i in range(3):
         transd[i] = transd[i] - math.floor(transd[i])
         #due to rounding error it may happen that the vector still lies outside of the unit cell
         #we have to take this into account
         if math.floor(transd[i]+prec) != 0: #if this occurs it means trands[i] lies within prec to 1 so we can se it to 0
             transd[i] = 0

     #transforms the momentum
     transd_m = []
     for i in range(3): #loop over components of magnetic moment
         op = re.sub('-','+-',sym[2][i])
         transd_m.append(0)
         op = re.split('\+',op) #op contains a list of all operations we have to do with component
         op = [_f for _f in op if _f] #remove empty strings from the list  
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
                     op[j] = old_div(float(op_s[0]), float(op_s[1]))
                 transd_m[i] = transd_m[i] + float(op[j])

     return list(transd+transd_m)


def r_basis(lines):
#read basis transformation
     for i in range(len(lines)):
         if "Vectors a,b,c:" in lines[i]:
             pos_vectors = i

     vec1 = lines[pos_vectors + 1]
     vec2 = lines[pos_vectors + 2]
     vec3 = lines[pos_vectors + 3]

     vec_a = [float(x) for x in vec1.split()]
     vec_b = [float(x) for x in vec2.split()]
     vec_c = [float(x) for x in vec3.split()]

     return [vec_a,vec_b,vec_c]


def r_origin(lines):
     out = False
     #finds the origin
     for i in range(len(lines)):
         if "Origin at" in lines[i]:
             out = lines[i].split()[2:5]
     for i in range(3):
         out[i] = float(out[i])
     #if we don't find the origin, we assum it's at 0,0,0    
     if not out:
         out = [0,0,0]
     return out
                 


def  r_pos(lines, fix_m=[]):

    #read atomic positions
    #format: list containing six floats for each position first three are the positions the other three magnetic moment
    #fix_m: m is not given in basis of a,b,c, but in basis of unit vectrs along these axis:
        #m = m1*a/|a|+m2*b/|b|+m3*c/|c|
    #since it is more interesting to have m in terms of a,b,c, this can be changed
    #needs magnitudes of a,b,c as input: fix_m=[a,b,c]
    
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
        at_number = int(positions[i].pop(0))
        for l in range(len(positions[i])):
            positions[i][l] = float(positions[i][l])
        positions[i].append(at_number)

    if fix_m:
        for i in range(len(positions)):
            for l in range(3):
                positions[i][l+3] = old_div(positions[i][l+3],fix_m[l])

    return positions


def r_sym(lines,debug=False,syms_only=False,num_prec=None):
#read symmetries 
#format: number of symmetry, space transformation, time transformation, switch controlling
#if there is time-reversal, list of tuples which show to which position each position transforms
#under this symmetry
     if num_prec is None:
         #tolerance = r_tolerance(lines)
         num_prec = 1.1e-5
     debug = False
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

     if debug:
         print('found symmetries:')
         for sym in syms:
             print(sym)

     if not syms_only:
         #read space group number
         for i in range(len(lines)):
             if "Space Group " in lines[i]:
                 group = re.findall('[0-9]+\.[0-9]+',lines[i])

         positions = r_pos(lines)
         check_pos_prec(positions,num_prec)

         #magnetic group tables
         #needs a file 'tables_wyckoff.txt' which are magnetic tables with all unneccessary information deleted
         dirname, filename = os.path.split(os.path.abspath(__file__))
         tables_loc = str(dirname)+'/findsym/tables_wyckoff.txt'
         tables=open(tables_loc)

         #read shifts of wickoff positions
         #this is neccessary to obtain all positions
         #in body centered structure for example, the output from findsym only contains half of the positions,
         #but symmetries also consider those that trasform one atom into another one shifted by 1/2 1/2 1/2
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
                         shift_sep[i] = old_div(float(t[0]), float(t[1]))
                     else:
                         shift_sep[i] = float(shift_sep[i])
                 shifts[l] = shift_sep

             if shifts[0] != [0,0,0]:
                 sys.exit('First shift is not zero')
             else:
                 shifts.pop(0)

                 if debug:
                     print('found wyckoff shifts:')
                     for shift in shifts:
                         print(shift)
         

         if debug:
             print('finding info about sublattice transformations')
         #this adds the information about transformation of sublattices 
         for l in range(len(syms)):
             if debug:
                 print('')
                 print('taking symmetry', syms[l])
             sym_trans = []
             for i in range(len(positions)):
                 trans = transform_position(positions[i],syms[l],num_prec)
                 if debug:
                     print('taking position:', positions[i])
                     print('transformed to:', trans)
                 found_atoms = []
                 for j in range(len(positions)):
                     if equal_vectors(trans[0:3],positions[j][0:3],num_prec):
                         found_atoms.append(positions[j][6])
                         #sym_trans.append((positions[i][6],positions[j][6]))
                         if debug:
                             print('transformed atom identified as atom ', positions[j][6], ' with position ', positions[j])
                     else:
                         for k in range(len(shifts)):
                             sym_temp =['',['x+'+str(shifts[k][0]),'y+'+str(shifts[k][1]),'z+'+str(shifts[k][2])],['mx','my','mz']]
                             pos_shift = transform_position(positions[j][0:6],sym_temp,num_prec)
                             if equal_vectors(trans[0:3],pos_shift[0:3],num_prec):
                                 #sym_trans.append((positions[i][6],positions[j][6]))
                                 found_atoms.append(positions[j][6])
                                 if debug:
                                     print('transformed atom identified as atom ', positions[j][6], ' with position ', positions[j])
                 if len(found_atoms) != 1:
                    print(found_atoms)
                    raise Exception('Problem with identifying transformation of atom {0}. '
                            'Try changing --pos-prec.'.format(positions[i][6]))
                 else:
                    sym_trans.append((positions[i][6],found_atoms[0]))
                        

             if len(sym_trans) != len(positions):
                 print(syms[l])
                 print(sym_trans)
                 raise Exception('Wrong number of transformed atoms. Try modifying --pos-prec parameter.')
             syms[l].append(list(sym_trans))

     return syms


def r_abc(lines):

     for i in range(len(lines)):
         if "Values of a,b,c,alpha,beta,gamma" in lines[i]:
             pos_abc = i

     abc = lines[pos_abc+1].split()

     return abc

def r_mag_fin(fin):
    start = False
    mags = []
    for i in range(len(fin)):
        if start:
            mag = sympy.zeros(1,3)  
            for j in range(3):
                mag[j] = sympy.sympify(fin[i].split()[j+3])
            mags.append(mag)
        elif (not start) and 'magnetic' in fin[i]:
            start = True
    return mags

def r_tolerance(lines):
    for line in lines:
        if 'Tolerance:' in line:
            num_prec = float(line.split(':')[1])
    if abs(num_prec) < 1e-14:
        num_prec = 1e-6
    return num_prec

def check_pos_prec(positions,num_prec):
    for i,pos in enumerate(positions):
        for j,pos2 in enumerate(positions):
            if i != j:
                if equal_vectors(pos[0:3],pos2[0:3],num_prec):
                    raise Exception(
                            'Atoms {0} and {1} too close together try lower --pos-prec'.format(pos[6],pos2[6]))


