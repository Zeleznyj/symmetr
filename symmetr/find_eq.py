from __future__ import print_function
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from builtins import range
from builtins import object

from collections import OrderedDict
import sympy

class confs(object):
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
                c = self.confs[i][j] - conf[j]
                for l in range(3):
                    if c[l].round(4)!= 0:
                        out_t = False
            if out_t:
                out = True

        return out

    def pprint(self,eo=None,m=-1,latex=False,print_format=None):
        #prints everything in a (somewhat) nice form
        #if m is set it prints only configuration m

        if m == -1:
            ran = list(range(self.nconfs))
        else:
            ran = [m]
        
        for n in ran:
            if n== 0:
                print('starting configuration:')
            else:
                print('configuration %s' %n)
            for p in sorted(self.confs[n].keys()):
                print('atom %s, m = %s, %s, %s' %(p,self.confs[n][p][0].round(4),self.confs[n][p][1].round(4),self.confs[n][p][2].round(4)))
                #print 'atom %s, m = %s, %s, %s' %(p,self.confs[n][p][0],self.confs[n][p][1],self.confs[n][p][2])
            print('First part of the response tensor')
            self.Xs[n][0].pprint(print_format=print_format)
            if latex:
                self.Xs[n][0].pprint(latex=True)
            print('Second part of the response tensor')
            self.Xs[n][1].pprint(print_format=print_format)
            if latex:
                self.Xs[n][1].pprint(latex=True)
            print('')

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

def find_equiv(Xs,mag,syms,atom,debug=False,round_prec=None):
    """
    Takes a tensor and a list of nonmagnetic symmetries and find the form of the tensor for all equivalent configurations.

    !!!Will only work if the same input basis was used for the nonmagnetic and magnetic!!!

    Args:
        X: The tensor for the starting configuration.
        op1: First operator type.
        op2: Second operator type.
        atom: Sets a projection to an atom.
        syms: The nonmagnetic symmetry operations.
        mag: The magnetic moments - ie the starting configuration.
        T: is the transformation matrix from the nonmagnetic basis to the input basis

    Returns:
        C (confs class): Contains all the equivalent configurations and the transformed tensors for each.
    """
    
    if debug:
        print('starting find_equiv')

    n_at = []
    for sym in syms:
        n_at.append(len(sym.permutations))
    n_at = list(set(n_at))
    if len(n_at) != 1:
        raise Exception('Different number of atoms permutations for different symmetries. Something wrong.')
    else:
        n_at = n_at[0]

    #extracts the starting configuration, only the magnetic moments are needed
    start_conf = OrderedDict()
    for i in range(len(mag)):
        a = mag[i][0]
        b = mag[i][1]
        c = mag[i][2]
        if a.round(4) !=0 or b.round(6) !=0 or c.round(6) != 0:
            start_conf[i+1] = sympy.Matrix([[a,b,c]])

    if n_at != len(start_conf):
        message = ('Warning: nonmagnetic unit cell smaller than the magnetic!\n'
                  'Output will show values of only some magnetic moments.'
                  'Furthermore some equivalent configurations will be missing.'
                  'Use with caution, this is not properly tested.')
        print('')
        print(message)
        start_conf2 = OrderedDict()
        for p in start_conf:
            if p in syms[0].permutations:
                start_conf2[p] = start_conf[p]
        start_conf = start_conf2

    #creates a conf class, which stores all the configurations and adds the starting one
    C = confs()
    C.add(start_conf,Xs)

    if debug:
        C.pprint(0)

    for sym in syms: #a loop over all symmetries
        
        #if there is a projection take only the symmetry that keeps the atom invariant
        if atom == -1 or sym.permutations[atom] == atom:

            if debug:
                print('')
                print('taking sym : ')
                print(sym)
                print('')

            #transforms the starting configuration by the symmetry
            conf_t = OrderedDict()
            for p in start_conf:
                mom = sympy.Matrix([[start_conf[p][0],start_conf[p][1],start_conf[p][2]]])
                mom_T = sym.get_R('s') * mom.T
                mom_T = mom_T.T
                conf_t[sym.permutations[p]] = sympy.Matrix([[mom_T[0],mom_T[1],mom_T[2]]])

            if debug:
                print('symmetry transforms starting configuration to configuration:')
                for p in conf_t:
                    print('atom %s, m = %s, %s, %s' %(p,conf_t[p][0],conf_t[p][1],conf_t[p][2]))
            
            #if the configuration has not been found before it transforms the tensor by the symmetry
            #operation and adds everything to C
            if not C.is_in(conf_t):
                if debug:
                    print('configuration has not been found before')
                Xt = []
                Xt.append(Xs[0].transform(sym))
                Xt.append(Xs[1].transform(sym))
                if debug:
                    print('even part converted to:')
                    Xt[0].pprint()
                    print('odd part converted to:')
                    Xt[1].pprint()
                if round_prec is not None:
                    for X in Xt:
                        X.round(round_prec)
                C.add(conf_t,Xt)
    
    return C

