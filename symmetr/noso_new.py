import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
from symmetr.symmetry import Symmetry
import sympy

def merge_lists(lists):
    if lists is None:
        return None
    merged = []
    for list in lists:
        if list is None:
            merged.append(None)
        else:
            merged += list
    return merged

class RotationType:
    """
    Class that is used to represent the various rotations types.
    typ:
        para: Arbitrary rotation around a given axis
        anti-para: 180 degree rotation around any axis perpendicular to the given axis
        single: a rotation given by axis and angle
        non-mag: completely arbitrary rotation
        identity: Identity matrix, i.e. zero angle rotation
    """

    def __init__(self, typ, axis=None, angle=None):
        if typ not in ['para', 'anti-para', 'single', 'non-mag','identity']:
            raise Exception('Rotation has to be one of "para", "anti-para", "single"!')
        self.typ = typ
        if self.typ not in ['non-mag','identity']:
            if axis is None:
                raise Exception('axis needs to be specified!')
            self.axis = axis / norm(axis)
            if self.typ in ['para', 'anti-para']:
                self.angle = None
                if angle is not None:
                    raise Warning('angle not need for parallel rotation and anti-parallel rotation.')
            else:
                if angle is None:
                    raise Exception('Angle needed for single rotation')
                self.angle = angle

    def __repr__(self):
        out = 'Rotation type: {}\n'.format(self.typ)
        if self.typ == 'non-mag':
            out += 'Arbitrary rotation\n'
        elif self.typ == 'identity':
            out += 'Identity (zero angle rotation)\n'
        elif self.typ in ['para', 'single']:
            out += 'Rotation axis: {}\n'.format(self.axis)
        else:
            out += 'Rotation axis: perpendicular to {}\n'.format(self.axis)

        if self.typ == 'single':
            out += 'Rotation angle: {}'.format(self.angle)
        elif self.typ == 'anti-para':
            out += 'Rotation angle: {}'.format(np.pi)
        elif self.typ == 'para':
            out += 'Rotation angle: arbitrary'
        return out

    def eq(self,other,prec):
        if self.typ != other.typ:
            return False
        if self.typ in ['identity', 'non-mag']:
            return True
        if self.typ in ['para', 'anti-para']:
            if norm(self.axis - other.axis) < prec:
                return True
            if norm(self.axis + other.axis) < prec:
                return True
            return False
        if self.typ == 'single':
            if norm(self.axis - other.axis) < prec:
                if abs(self.angle - other.angle ) < prec:
                    return True
                return False
            if norm(self.axis + other.axis) < prec:
                if abs(self.angle + other.angle - 2 * np.pi) < prec:
                    return True
                return False
            return False
        raise Exception('Wrong rot type')

class NosoSymFinder:
    """
    This class contains the functions that are used for finding the non-relativistic symmetries. They are grouped
    together in a class just so that all have access to the prec and debug parameters.

    prec: General precision parameter used throughout to determine if two floats are the same
    moment_zero: A precision parameter that is used to determine if a magnetic moment is zero. Can be used to ignore
    small moments.
    debug: Debug 0 means no debug output, 1 means standard debug output, 2 means detailed debug.
    """

    def __init__(self,prec=1e-5,moment_zero=1e-3,debug=0):
        self.prec = prec
        self.moment_zero = moment_zero
        self.debug = debug

    def get_i_chain(self,perm,chain):
        pi = perm[chain[-1]]
        if pi == chain[0]:
            return chain
        else:
            chain.append(pi)
            return self.get_i_chain(perm,chain)

    def get_i_chain_nr(self,perm,start):
        chain = []
        chain.append(start)
        finished = False
        current = start
        while not finished:
            next = perm[current]
            if next not in chain:
                chain.append(next)
                current = next
            elif next == chain[0]:
                #chain.append(next)
                finished = True
                chain_ok = True
            else:
                #chain.append(next)
                finished = True
                chain_ok = False
        return chain, chain_ok

    def get_permutation_chains(self,perm):
        """
        Get all permutation chains for given permutation.
        """
        classified = []
        chains = []
        oks = []
        for i in perm:
            if i not in classified:
                chain, ok = self.get_i_chain_nr(perm,i)
                chains.append(tuple(chain))
                oks.append(ok)
                classified += chain
        ok_all = False not in oks
        return tuple(chains), ok_all

    def get_AB_rotations(self,A, B, theta):
        """
        Finds all rotations that transform A to B by a given angle theta.
        """

        An = A / norm(A)
        Bn = B / norm(B)
        para = (norm(A - B) < self.prec)
        antipara = (norm(A + B) < self.prec)

        if para:
            return [RotationType('para', An)]

        if abs(theta) < self.prec:
            return None
        else:
            if antipara:
                if abs(theta - np.pi) > self.prec:
                    return None
                else:
                    return [RotationType('anti-para', An)]
            else:
                if np.dot(An, Bn) < np.cos(theta) - 1e-11:
                    if self.debug > 1:
                        print('angle not compatible')
                        print(np.dot(An, Bn), np.cos(theta))
                    return None

                sqrtt = 2 * (np.dot(An, Bn) - np.cos(theta)) / (1 - np.cos(theta)) / (1 + np.dot(An, Bn))
                if abs(sqrtt) < 1e-12:
                    sqrtt = 0
                rots = []
                for i,sign in enumerate([1,-1]):
                    acos = sign*np.sqrt(sqrtt)
                    if abs(acos) > 1:
                        if self.debug > 1:
                            print(acos)
                        acos = np.sign(acos)
                    alpha = np.arccos(acos)
                    nt = np.cos(alpha) * (A + B) / norm(A + B) + np.sin(alpha) * np.cross(A, B) / norm(np.cross(A, B))

                    r = R.from_rotvec(nt * theta)
                    if norm(r.apply(A) - B) > self.prec:
                        r2 = R.from_rotvec(-nt * theta)
                        if norm(r2.apply(A) - B) < self.prec:
                            if self.debug > 1:
                                print('switching axis')
                            rots.append(RotationType('single', -nt, theta))
                        else:
                            raise Exception('Rotation test failed, somethign is wrong!')
                    else:
                        rots.append(RotationType('single', nt, theta))
                rots = self.remove_repeated(rots)
                return rots


    def get_rotations_overlap(self,r1, r2):
        """
        Finds an overlap between two rotation types.
        """
        if r1 is None or r2 is None:
            return None
        if r1.typ in ['non-mag']:
            return r2
        if r2.typ in ['non-mag']:
            return r1
        if r1.typ == 'identity' or r2.typ == 'identity':
            if r1.typ != 'identity' and r2.typ == 'identity':
                r1i = r2
                r2i = r1
            else:
                r1i = r1
                r2i = r2
            if r2i.typ == 'identity':
                return r1i
            elif r2i.typ == 'para':
                return r1i
            else:
                return None
        if r1.typ != 'single' and r2.typ == 'single':
            r1i = r2
            r2i = r1
        elif (r1.typ != 'single' and r2.typ != 'single') and (r1.typ != 'para' and r2.typ == 'para'):
            r1i = r2
            r2i = r1
        else:
            r1i = r1
            r2i = r2

        if r1i.typ == 'single':

            if r2i.typ == 'single':
                if norm(r1i.axis - r2i.axis) < self.prec:
                    if abs(r1i.angle - r2i.angle) < self.prec:
                        return r1
                if norm(r1i.axis + r2i.axis) < self.prec:
                    if abs(r1i.angle + r2i.angle - 2 * np.pi) < self.prec:
                        return r1i

            elif r2i.typ == 'para':
                if norm(r1i.axis - r2i.axis) < self.prec or norm(r1i.axis + r2i.axis) < self.prec:
                    return r1i

            elif r2i.typ == 'anti-para':
                if abs(r1i.angle - np.pi) < self.prec and abs(np.dot(r1i.axis, r2i.axis)) < self.prec:
                    return r1i

            else:
                raise Exception('Wrong r2i.type')

        elif r1i.typ == 'para':

            if r2i.typ == 'para':
                if norm(r1i.axis - r2i.axis) < self.prec or norm(r1i.axis + r2i.axis) < self.prec:
                    return r1i
                else:
                    return RotationType('identity')

            elif r2i.typ == 'anti-para':
                if abs(np.dot(r1i.axis, r2i.axis)) < self.prec:
                    return RotationType('single', r1i.axis, np.pi)
            else:
                raise Exception('Wrong r2i.type')

        elif r1i.typ == 'anti-para':

            if r2i.typ == 'anti-para':
                if norm(r1i.axis - r2i.axis) < self.prec or norm(r1i.axis + r2i.axis) < self.prec:
                    return r1i
                else:
                    return RotationType('single',np.cross(r1i.axis,r2i.axis),np.pi)
            else:
                raise Exception('Wrong r2i.type')

        else:
            raise Exception('Wrong r1i.type')

        return None

    def get_rotations_overlap_multi(self, r1, r2):
        """
        Finds an overlap between two lists of rotations.
        That is it finds all overlaps between r1i and r2j where r1i are rotations from r1 and r2j rotations from r2.
        """
        if self.debug > 1:
            print('rotation_overlap_multi:')
            print('r1',r1)
            print(type(r2))
            print('r2',r2)
        overlaps = []
        if r1 is None or r2 is None:
            return None
        for ri in r1:
            for rj in r2:
                r = self.get_rotations_overlap(ri,rj)
                if r is not None:
                    overlaps.append(r)
        if len(overlaps) > 0:
            overlaps = self.remove_repeated(overlaps)
            return overlaps
        else:
            return None

    def remove_repeated(self,rots):
        if rots is None:
            return None
        if len(rots) == 0:
            return rots
        rots_u = []
        rots_u.append(rots.pop(0))
        for rot in rots:
            is_in = False
            for rot2 in rots_u:
                if rot.eq(rot2,self.prec):
                    is_in = True
            if not is_in:
                rots_u.append(rot)
        return rots_u

    def merge_rotations(self,rots):
        """
        Takes a list of rotations, where each element is a list of rotations and find a common overlap.
        """
        if self.debug > 1:
            print('starting merge')
            print(rots)

        #This part seems wrong and unnecessary
        #rots_all = merge_lists(rots)
        #if any(x is None for x in rots_all):
        #    if self.debug > 1:
        #        print('finished merge')
        #    return None
        rots_merged = None
        if len(rots) == 1:
            rots_merged = rots[0]
        else:
            if self.debug > 1:
                print('starting get_rotations_overlap_multi')
            ro = self.get_rotations_overlap_multi(rots[0],rots[1])
            if self.debug > 1:
                print('finished get_rotations_overlap_multi')
            if ro is None:
                if self.debug > 1:
                    print('finished merge')
                return None
            rots_merged = ro
            for irot in range(2,len(rots)):
                    if self.debug > 1:
                        print('starting get_rotations_overlap_multi')
                    ro = self.get_rotations_overlap_multi(rots_merged,rots[irot])
                    if self.debug > 1:
                        print('finished get_rotations_overlap_multi')
                    if ro is None:
                        if self.debug > 1:
                            print('finished merge')
                        return None
                    else:
                        rots_merged = ro
        if self.debug > 1:
            print('finished merge')
        return rots_merged


    def get_chain_rotations(self, mchain, hasT=False):
        """
        Finds spin rotations for a given chain of magnetic moments.
        """

        if self.debug > 1:
            print('Magnetic chain: ', mchain)

        rotations = []
        para = False
        if not hasT:
            col_dif = max([norm(m - mchain[0]) for m in mchain])
            if col_dif < self.prec:
                para = True
        elif len(mchain) % 2 == 0:
            mchain_f = [2 * (i % 2 - 1 / 2) * m for i, m in enumerate(mchain)]
            col_dif = max([norm(m - mchain_f[0]) for m in mchain_f])
            if col_dif < self.prec:
                para = True

        if para:
            rotations.append(RotationType('para', mchain[0]))
            return rotations

        p = len(mchain)
        if not hasT or len(mchain) % 2 == 0:
            thetas = [2 * np.pi * i / p for i in range(1, int(np.floor(p / 2)) + 1)]
        else:
            thetas = [np.pi * (2 * i + 1) / p for i in range(0, int(np.ceil(p / 2)))]
        for theta in thetas:

            if self.debug > 1:
                print('theta: ', theta)

            rotations_theta = []
            for im in range(len(mchain)):
                im2 = im + 1
                if im2 > len(mchain) - 1:
                    im2 = 0
                if not hasT:
                    f = 1
                else:
                    f = -1
                rot = self.get_AB_rotations(mchain[im], f * mchain[im2], theta)
                if self.debug > 1:
                    print('A: ', mchain[im], 'B: ', mchain[im2])
                    print(rot)
                rotations_theta.append(rot)

            rotations_theta_merged = self.merge_rotations(rotations_theta)
            if rotations_theta_merged is not None:
                rotations.append(rotations_theta_merged)
        rotations_all = merge_lists(rotations)
        return rotations_all

    def get_mchain(self,chain,moments):
        mchain = [moments[i-1] for i in chain]
        mchain_norms = [norm(m) for m in mchain]
        if max(mchain_norms) < self.moment_zero:
            return 0
        if min(mchain_norms) < max(mchain_norms) - 2 * self.moment_zero:
            #raise Exception('Chain contains moments with different magnitudes!')
            return None
        return [m/mchain_norms[i] for i,m in enumerate(mchain)]

    def get_permutations_rotations(self, chains, mags, hasT=False):
        """
        Finds spin rotations for a give permutation, which is represented by the chains.
        """
        # chains = get_permutation_chains(perm)
        rots_chains = []
        for chain in chains:
            mchain = self.get_mchain(chain, mags)
            if self.debug > 0:
                print('chain: ', chain)
                print('mchain: ', mchain)
            if mchain == 0:
                rots_chains.append([RotationType('non-mag')])
            elif mchain is None:
                return None
            else:
                rots = self.get_chain_rotations(mchain, hasT=hasT)
                if self.debug > 0:
                    print(rots)
                #if len(rots) > 1:
                #    raise Exception("More than one rotation per chain, this is not expected behavior, something is wrong.")
                if len(rots) == 0:
                    return None
                else:
                    rots_chains.append(rots)
        rot_merged = self.merge_rotations(rots_chains)

        return rot_merged

    def find_noso_syms(self, syms, mags, n_discretize=4):
        """
        Determinese the non-relativistic symmetry operations of a magnetic system.

        If the system contains continuous spin-rotation symmetry, we discretize it. This is not a perfect solution
        but should be fine in most cases. This only concerns collinear systems.

        Args:
            syms ([class Symmetry]): List of the non-magnetic symmetry operations of the system. Must be in
                some cartesian basis
            mags ([np.array]): List of the magnetic moments of the system. Must be in some cartesian basis.

        returns:
            noso_syms: List of the non-relativistic symmetries.
        """

        nonmag = False
        collinear = False
        if max([norm(m) for m in mags]) < self.moment_zero:
            nonmag = True
        else:
            non_zero_mags = []
            for m in mags:
                if norm(m) > self.moment_zero:
                    non_zero_mags.append(m)
            if self.debug > 0:
                print(non_zero_mags)
            #dots = [abs(np.dot(non_zero_mags[0], m) / norm(non_zero_mags[0]) / norm(m)) for m in non_zero_mags]
            #if min(dots) > 1 - self.prec:
            #    collinear = True
            non_zero_mags_norm = [m/norm(m) for m in non_zero_mags]
            col_dif = max([min(norm(m + non_zero_mags_norm[0]), norm(m - non_zero_mags_norm[0])) for m in non_zero_mags_norm])
            if col_dif < self.prec:
                collinear = True

        syms_noso = []
        if nonmag:
            raise Exception('Not implemented')
        elif collinear:
            col_axis = non_zero_mags[0] / norm(non_zero_mags[0])
            para_rot = RotationType('para', axis=non_zero_mags[0] / norm(non_zero_mags[0]))
            triv_perm = {i: i for i in syms[0].permutations}
            sym_I = Symmetry(sympy.diag(1, 1, 1), False, permutations=triv_perm)
            if self.debug > 0:
                print('Collinear system, axis: {}'.format(col_axis))
            if n_discretize is None:
                sym_I.Rs = para_rot
                syms_noso.append(sym_I)
            else:
                for theta in np.linspace(0, np.pi, n_discretize, endpoint=True):
                    sym_noso = sym_I.copy()
                    r = R.from_rotvec(para_rot.axis * theta)
                    sym_noso.Rs = sympy.Matrix(r.as_matrix())
                    syms_noso.append(sym_noso)

        if self.debug > 0:
            print('Determining spin rotations')
            print(' ')
        chains_all = []
        all_ok = True
        for sym in syms:
            chains, ok =  self.get_permutation_chains(sym.permutations)
            if ok:
                chains_all.append(chains)
            else:
                all_ok = False
                if self.debug > 0:
                    print('WARNING Symmetry that is not compatible with the magnetic unit cell is present!'
                                 'Skipping this symmetry!')
                    print(sym)
                    print(chains)
        if not all_ok:
            print('WARNING Symmetry that is not compatible with the magnetic unit cell is present!')

        chains_u = list(set(chains_all))
        rots = []
        rotsT = []
        for i, chains in enumerate(chains_u):
            if self.debug > 0:
                print('')
                print('taking permutation chains {}'.format(chains))
            # print(chains)
            if self.debug > 0:
                print('')
                print('spin-rotation without T')
            rot = self.get_permutations_rotations(chains, mags, False)
            rot = self.remove_repeated(rot)
            if self.debug > 0:
                print('Compatible rotation: ',rot)
                print('spin-rotation with T')
            rotT = self.get_permutations_rotations(chains, mags, True)
            rotT = self.remove_repeated(rotT)
            rots.append(rot)
            rotsT.append(rotT)
            if self.debug > 0:
                print('Compatible rotation: ', rotT)

        for i, sym in enumerate(syms):
            chains,ok = self.get_permutation_chains(sym.permutations)
            if not ok:
                continue
            index = chains_u.index(chains)
            if self.debug > 0:
                print('')
                print('Taking symmetry {}'.format(i))
                print(sym)
                print('')
            if not sym.has_T:
                rot = rots[index]
            else:
                rot = rotsT[index]
            if rot is None:
                if self.debug > 0:
                    print('No compatible spin rotation')
                continue
            else:
                if self.debug > 0:
                    print('Compatible spin rotation:')
                    print(rot)
            if collinear:
                if len(rot) > 1:
                    raise Exception('More than one rotation for collinear, this should not happen.')
                rot = rot[0]
                if rot.typ == 'single' or rot.typ == 'non-mag' or rot.typ == 'identity':
                    raise Exception('Single, non-mag or identity rotation present for collinear system, something is wrong')
                else:
                    sym_noso = sym.copy()
                    if sym.has_T:
                        f = sympy.diag(-1,-1,-1)
                    else:
                        f = sympy.diag(1,1,1)
                    if rot.typ == 'para':
                        sym_noso.Rs = f * sympy.diag(1, 1, 1)
                    else:
                        rot_axis = None
                        for rot_axis_c in [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
                                           np.array([1, 1, 0])]:
                            if abs(np.dot(rot.axis, rot_axis_c)) < self.prec:
                                rot_axis = rot_axis_c
                                break
                        if rot_axis is None:
                            rot_axis = np.cross(rot.axis, [0, 0, 1])
                            rot_axis = rot_axis / norm(rot_axis)
                        sym_noso.Rs = f * sympy.Matrix(R.from_rotvec(rot_axis * np.pi).as_matrix())
                    syms_noso.append(sym_noso)
            else:
                if len(rot) > 1:
                    raise Exception('More than one rotation, this is not implemented (and may suggest somethings wrong).')
                rot = rot[0]
                if rot.typ not in ['single','identity']:
                    raise Exception('Non-collinear system and non-single rotation present, someting is wrong!')

                if rot.typ == 'single':
                    r = R.from_rotvec(rot.axis * rot.angle).as_matrix()
                else:
                    r = np.diag([1,1,1])
                sym_noso = sym.copy()
                if not sym.has_T:
                    sym_noso.Rs = sympy.Matrix(r)
                else:
                    sym_noso.Rs = sympy.Matrix(np.matmul(np.diag([-1, -1, -1]), r))
                syms_noso.append(sym_noso)

        return syms_noso

