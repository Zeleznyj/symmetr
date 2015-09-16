import re
import copy
import sys
import math

import sympy
import numpy as np

from tensors import matrix, mat2ten
from read import transform_position
from funcs import *

class confs:
    """
    Class to store equivalent magnetic configurations and their response tensors.

    Description:
        self.nconfs is the number of configurations stored
        self.confs[i] is the configuration (stored as a dict) with number i
        self.Xs[i] is the tensor of the configuration i

    Usage:
        __init__ creates an empty object. 
        Add a configuration by add(conf,X).
        Print all by pprint.
    """

    def __init__(self):
        """Creates an empty confs object."""

        self.confs = {}
        self.Xs = {}
        self.nconfs = 0

    def add(self,conf,X):
        """
        Adds a new configuration.
        
        Args:
            conf (dictionary): The magnetic configuration: {atom number: magnetic moment}
            X: The tensor of that configuration
        """

        self.confs[self.nconfs] = conf
        self.Xs[self.nconfs] = X
        self.nconfs += 1
        
    def is_in(self,conf):
        #tests if a configuration is in
            
        out = False
        for i in range(self.nconfs):
            out_t = True
            for j in self.confs[i]:
                if not np.allclose(self.confs[i][j],conf[j]):
                    out_t = False
            if out_t:
                out = True

        return out

    def pprint(self,m=-1,latex=False):
        #prints everything in a (somewhat) nice form
        #if m is set it prints only configuration m

        if m == -1:
            ran = range(self.nconfs)
        else:
            ran = [m]
        
        for n in ran:
            if n== 0:
                print 'starting configuration:'
            else:
                print 'configuration %s' %n
            for p in self.confs[n]:
                print 'atom %s, m = %s, %s, %s' %(p,self.confs[n][p][0],self.confs[n][p][1],self.confs[n][p][2])
            print 'even part:'
            self.Xs[n][0].pprint()
            if latex:
                self.Xs[n][0].pprint(latex=True)
            print 'odd part:\n'
            self.Xs[n][1].pprint()
            if latex:
                self.Xs[n][1].pprint(latex=True)
            print ''

    def convert(self,T,shift):
        #converts to a different coordinate system
        #seems to work, but not fully tested, not sure if the shift will work correctly
        #not actually used anywhere, but may be useful in the future

        Tt = mat2ten(T)    

        confs_t = confs()
        for n in range(self.nconfs):

            conf_t = {}
            for j in self.confs[n]:
                pos = self.confs[n][j]
                pos_t = np.dot(T,pos) + np.array(shift)
                conf_t[j] = pos_t

            X = self.Xs[n]
            X_t = []
            X_t.append(Tt*X[0]*Tt.T())
            X_t.append(Tt*X[1]*Tt.T())

            confs_t.add(conf_t,X_t)

        return confs_t


def find_equiv(X,op1,op2,atom,syms,pos,T,shift,debug=False):
    """
    Takes a tensor and a list of nonmagnetic symmetries and find the form of the tensor for all equivalent configurations.

    !!!Will only work if the same input basis was used for the nonmagnetic and magnetic!!!

    Args:
        X: The tensor for the starting configuration.
        op1: First operator type.
        op2: Second operator type.
        atom: Sets a projection to an atom.
        syms: The nonmagnetic symmetry operations.
        pos: The positions and magnetic moments - ie the starting configuration.
        T: is the transformation matrix from the nonmagnetic basis to the input basis

    Returns:
        C (confs class): Contains all the equivalent configurations and the transformed tensors for each.
    """
    
    if debug:
        print 'starting find_equiv'

    #extracts the starting configuration, only the magnetic moments are needed
    start_conf = {}
    for i in range(len(pos)):
        a = round(float(pos[i][3]),5)
        b = round(float(pos[i][4]),5)
        c = round(float(pos[i][5]),5)
        if a!=0 or b !=0 or c != 0:
            start_conf[pos[i][6]] = np.array([a,b,c])

    #creates a conf class, which stores all the configurations and adds the starting one
    C = confs()
    C.add(start_conf,X)

    if debug:
        C.pprint(0)

    for sym in syms: #a loop over all symmetries
        
        #if there is a projection take only the symmetry that keeps the atom invariant
        if atom == -1 or atom == sym_type(atom,sym):

            #convert the symmetry to the matrix form and transform to the input basis
            sym_mat = sym2mat(sym,op_type='s')
            sym_mat_T = convert_sym_mat(sym_mat,T)

            if debug:
                print 'taking sym (in the nonmagnetic basis): '
                print sym
                print ''
                print 'in matrix form:'
                sympy.pprint(sym_mat)
                print 'in the input basis:'
                sympy.pprint(sym_mat_T)
                print ''

            #transforms the starting configuration by the symmetry
            conf_t = {}
            for p in start_conf:
                mom = sympy.Matrix([[start_conf[p][0],start_conf[p][1],start_conf[p][2]]])
                mom_T = sym_mat_T * mom.T
                mom_T = mom_T.T
                conf_t[sym_type(p,sym)] = np.array([round(mom_T[0],5),round(mom_T[1],5),round(mom_T[2],5)])

            if debug:
                print 'symmetry transforms starting configuration to configuration:'
                for p in conf_t:
                    print 'atom %s, m = %s, %s, %s' %(p,conf_t[p][0],conf_t[p][1],conf_t[p][2])
            
            #if the configuration has not been found before it transforms the tensor by the symmetry
            #operation and adds everything to C
            if not C.is_in(conf_t):
                if debug:
                    print 'configuration has not been found before'
                Xt = []
                for l in range(2):
                    Xt.append(transform_matrix(X[l],sym,op1,op2,l,T=T))
                if debug:
                    print 'even part converted to:'
                    Xt[0].pprint()
                    print 'odd part converted to:'
                    Xt[1].pprint()
                C.add(conf_t,Xt)
    
    return C


