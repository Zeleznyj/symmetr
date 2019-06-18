# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
"""
Module for finding a symmetrical form of classical magnetic Hamiltonians.
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from builtins import object
import sympy as sp

from .symmetrize import symmetr,SymmetrOpt
from .tensors import tensor

#testing
import sys
import os
import subprocess

import sympy

def convert_mag_ham(H,T):
    """Converts the tensor representing the magnetic Hamiltonian to a different basis.

    Note: This is just a tensor transformation, there are several other functions very
    similar to this. The only difference is which indeces transform by T and which by
    T^-1. This just corresponds to covariant and contravariant indeces and this would be
    nice to eventually implement in the code.

    Args:
        H (tensor): the magnetic Hamiltonian in basis A
        T (matrix): . If it is set, the symmetry operations
            will be transformed by this matrix. T transforms from A to B, ie Tx_A = x_B.
    Returns:
        H_T: the magnetic Hamiltonian in basis B
    """

    Ti = T.inv()
    H_T = H.copy0()
    for ind1 in H:
        for ind2 in H_T:
            factor = 1
            for i in range(H.dim2):
                    factor *= Ti[ind1[i],ind2[i]]
            H_T[ind1] += factor*H[ind2]
    return H_T

class params_trans_ham(object):
    """Class for storing parameters of the transformation function."""
    def __init__(self,sites,debug=None,T=None,check_sym=False):
        self.sites = sites
        self.debug = debug
        self.T = T
        self.check_sym = check_sym

def trans_sites(sites,sym):
    sites_t = []
    perms = sym.permutations
    for s in sites:
        sites_t.append(perms[s])
    return sites_t

def trans_mag_ham(H,sym,params):
    """Transforms the magnetic Hamiltonian by a symmetry operation.

    There are two different modes depending on params.check_sym
    If params.chech_sym == True:
        Returns None for any symmetry operation, which transforms sites into a differet
        set of sites. If sym transform sites into the same set of sites then it returns
        transformed H for the original sites. That is:
        if sites = 1,2 and sym(sites) = 2,1 then sym is used and returned H is for
        sites = 1,2.
    If params.check_sym == False:
        Returns transformed H for the transformed sites. That is if sites = 1,2 
        and sym(sites) = 2,1 the H is returned for 2,1 sites. If sym(sites) = 2,3 then
        H is returned for 2,3 sites.

    Args:
        H (class tensor): the magnetic Hamiltonian.
        sym (): symmetry operation in the findsym format  
        params (class params_trans_ham): contains parameters for the transformation
    Returns:
        H_T (class tensor): transformed Hamiltonian
    """

    # this checks whether the symmetry operation leaves the atomic sites invariant
    # if it does then we use it, otherwise we skip it and so we return the unchanged
    # Hamiltonian
    sites = params.sites
    sites_t = trans_sites(sites,sym)
    if params.check_sym:
        if sorted(sites) != sorted(sites_t):
            return None

    R = sym.get_R('s') # select the matrix corresponding the spin transformation
    RiT = R.inv().T
    if params.debug:
        print('symmetry:')
        print(sym)
        print('symmetry in the matrix form')
        sp.pprint(R)

    if params.check_sym:
        # find how does the symmetry operation transforms the atomic sites
        sites_un = set(sites)
        perm = [0]*len(sites)
        pos_s = {}
        pos_st = {}
        for s in sites_un:
            pos_s[s] = [j for j,x in enumerate(sites) if x==s]
            pos_st[s] = [j for j,x in enumerate(sites_t) if x==s]
        for s in sites_un:
            for i,p in enumerate(pos_s[s]):
                perm[p] = pos_st[s][i] 

    if params.debug:
        print('sites:')
        print(sites)
        print('transformed sites:')
        print(sites_t)
        print('the permutation:')
        print(perm)

    # permutes the tensor index according to the permutatin of the atoms
    def perm_ind(ind,perm):
        ind_p = [0]*len(ind)
        for i,ii in enumerate(ind):
            ind_p[perm[i]] = ii
        return tuple(ind_p)

    # transforms the tensor, taking the permutatin of atoms into account
    H_t = H.copy0()
    for ind in H_t:
        for ind2 in H_t:
            factor = 1
            for i in range(len(ind)):
                factor *= RiT[ind[i],ind2[i]]
            if params.check_sym:
                ind_n = perm_ind(ind,perm)
            else:
                ind_n = ind
            H_t[ind_n] += factor*H[ind2]

    if params.debug:
        print('original tensor:')
        print(H)
        print('transformed tensor:')
        print(H_t)

    return H_t

def find_perms(sites):
    """Find permutations of identical atomic sites.

    The Hamiltonian must be invariant under interchanging atomic sites if
    they are the same. This function takes the list of sites and returns
    such permutations.
    Args:
        sites (list): list of atomic sites
    Returns:
        perms : list of permutations
    """
    perms = []
    sites_un = set(sites)
    pos_s = {}
    for s in sites_un:
        pos_s[s] = [j for j,x in enumerate(sites) if x==s]
        if len(pos_s[s]) > 1:
            for i,x in enumerate(pos_s[s]):
                for j in range(1,len(pos_s[s])):
                    if x != pos_s[s][j]:
                        perms.append((x,pos_s[s][j]))
    return perms

def trans_mag_Ham_perms(H,perm,params):
    """Transforms the Hamiltonian by a permutation of atomic sites.
    
    Args:
        H (class tensor): the magnetic Hamiltonian.
        perm (tuple): the permutation
        params: not used right now
    Returns:
        H_T: the transformed Hamiltonian
    
    """
    H_t = H.copy0()
    for ind in H:
        ind_t = list(ind)
        ind_t[perm[0]] = ind[perm[1]]
        ind_t[perm[1]] = ind[perm[0]]
        H_t[ind] = H[tuple(ind_t)]

    return H_t

def equiv(H,sites,syms,T=None,debug=False):
    """Transform H for one set of sites to all sites that are related by symmetry operation

    Args:
        H (class tensor): the magnetic Hamiltonian.
        sites (list): list of atomic sites
        syms (): symmetry operations in the findsym format
        T (sympy matrix): Coordinate transformation matrix used to transform the symmetry operations.
            Must transform syms to the coordinate system in which H is given.
    Returns:
        Asites (dictionary): the keys are all the equivalent sets of sites. Asites[sites] is the
            transformed Hamiltonian for the sites.
    """

    Asites = {}
    Asites[tuple(sorted(sites))] = H

    for sym in syms:
        sites_t = trans_sites(sites,sym)
        sites_th = tuple(sorted(sites_t))
        if sites_th not in Asites:
            params = params_trans_ham(sites,False,None)
            H_t = trans_mag_ham(H,sym,params)
            Asites[sites_th] = H_t

    return Asites

def sym_mag_ham(sites,syms,T=None,s_opt=None):
    """Returns the symmetrized magnetic Hamiltonian

    Args:
        sites (list): The atomic sites for which the Hamiltonian is considered
        syms : the symmetry operations.
        T (optinal[matrix]): Transformation matrix to a different basis. If set,
            it is used to transform the symmetry operations.
    Returns:
        The symmetrized form of the Hamiltonian
    """

    order = len(sites)
    ind_types = (-1,)*len(sites)
    H = tensor('s',3,order,ind_types=ind_types)
    
    params = params_trans_ham(sites,debug=s_opt.debug,T=T,check_sym=True)
    # first symmetrize with normal symmetry operation
    Hs = symmetr(syms,H,trans_mag_ham,params,s_opt)
    # The Hamiltonian must also be symmetric under any permutation of identical atoms
    perms = find_perms(sites)
    Hs = symmetr(perms,Hs,trans_mag_Ham_perms,params,s_opt)

    return Hs

def print_Ham(H,sites,latex=False):
    """Prints the magnetic Hamiltonian in a nice format.
    """
    sites_un = set(sites)
    sp_vars = {}
    Hsymb = 0
    for site in sites_un:
        for i in range(3):
            comps = ['x','y','z']
            varname = 'M_'+str(site)+comps[i]
            sp_vars[(site,i)] = sp.symbols(varname)
    for ind in H:
        Hp = H[ind]
        for i,x in enumerate(ind):
            Hp *= sp_vars[(sites[i],x)]
        Hsymb += Hp

    if not latex:
        sp.pprint(Hsymb)
    else:
        print(sp.latex(Hsymb))
