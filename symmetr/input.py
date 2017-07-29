"""Module for parsing user input."""

import argparse
import sys

class options:
    """Class to store all the input options.

    Is meant to be used only in conjuction with the parse function.
    Initialize by running
        opt = options(args),
    where args is a dictionary of all arguments.
    Check method controls whether there are comflicting options, however this is not complete and it does not
    check whether the options have allowed values!
    You can access the options by:
        opt['option']
    and also change by:
        opt['option'] = new_value
    """
    def __init__(self,args):
        """During init only the mandatory arguments are set, the rest are set to default values."""
        self.args = args
        self.check()
    def __setitem__(self,key,value):
        self.args[key] = value
        self.check()
    def __getitem__(self,key):
        return self.args[key]
    def __str__(self):
        return self.args.__str__()
    def check(self):
        if self['atom2'] != -1:
            if self['atom'] == -1:
                sys.exit('projection2 can be setonly if projection1 is set')

        if self['inp'] and self['group']:
            sys.exit('You cannot specify both the symmetry group and Findsym input file')

        if (not self['inp']) and (not self['group']):
            sys.exit('You have to specify either the symmetry group or the Findsym input file')

        if ( self['atom'] != -1 or self['atom2'] !=-1 ) and self['group']:
            sys.exit('Projections not possible with group name input. Use Findsym input instead.')

        if self['equiv'] and self['group']:
            sys.exit('Equivalent configurations are not possible with group name input. Use findsym input file.')

        if self['op3'] and self['equiv']:
            sys.exit('Equivalent configurations not implemented for three operators.')

        if self['op3'] and (self['exp'] != -1):
            sys.exit('Expansions are not implemented for three operators.')

        if self['op3'] and (self['equiv']):
            sys.exit('Equivalent magnetic configurations are not implemented for three operators.')

        if self['noso'] and self['group']:
            sys.exit('Spin-orbit coupling cannot be ignored when group name is used as an input.')

        if self['noso'] and (self['equiv'] or (self['exp'] != -1) or (self['atom2'] != -1)):
            sys.exit('This is not implemented.')

        if self['transform_result'] and self['transform_syms']:
            sys.exit('You cannot specify both --transform-result and --transform-syms')

        if (self['equiv'] or (self['exp'] != -1)) and self['syms_sel'] != -1:
            sys.exit('You cannot select symmetries when exp or equiv is set since this is not implemented')

        if self['equiv'] and ( self['exp'] != -1 ):
            sys.exit('You cannot use --equiv and --exp together. Equivalent configurations not supported for expansions.')
        if self['op1'] == self['op2'] and self['op3'] == None and self['group']:
            sys.exit('You have to set \'--ignore-op1e1op2\' when using group name input and two same operators \
                        since this is not implemented yet.')
        

def parse(clargs=None):
    """Parses the input line arguments

    If clargs is not set then this parses the input arguments. Otherwise clargs is a string that contains the input arguments
    like on command line.
    Args:
        clargs (optional[string]): The input arguments to be parsed.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('op1', help='Type of the first operator')
    parser.add_argument('op2', help='Type of the second operator')
    parser.add_argument('-p','--projection',help='Sets a projection on an atom.',default=-1,dest='atom',type=int)
    parser.add_argument('-p2','--projection2',help='Sets a projection on a second atom. Tries to find a relation between tensors on the first \
            atom and on the second atom.',default=-1,dest='atom2',type=int)
    parser.add_argument('-f','--findsym',help='Findsym input file',default=None,dest='inp')
    parser.add_argument('-g','--group',help='group name',default=None)
    parser.add_argument('-op3',help='third operator in the linear response formula',default=None)
    parser.add_argument('-b','--basis',help='Sets a coordinate basis: abc for conventional crystallographic basis, i for the one used in input \
    (default). cart for a cartesian basis in which the input basis is define. \
            abc_c for orthogonalized crystalographic basis (not tested much).',default='i')
    parser.add_argument('-e','--equivalent',action='store_const',const=True,default=False,help='finds response matrices for equivalent magnetic configurations.',dest='equiv')
    parser.add_argument('--no-rename',action='store_true')
    parser.add_argument('--debug',help='Controls if debug output is printed. all means all debug output is printed, symmetrize means debug\
            output for symmetrizing, rename for renaming, equiv for finding the equivalent configurations,\
            noso for the symmetry without spin-orbit coupling. op1eqop2 is also possible.',default='')
    parser.add_argument('--latex',action='store_const',const=True,default=False,help='If set, the matrices are printed also in a latex format.')
    parser.add_argument('--exp',default=-1,help=\
            'Prints the tensor, which describes the expansion of the linear response tensor in magnetic moments.\
            Choose which order term in the expansion.\
            Only works for ferromagnets and collinear antiferromagnets. In antiferromagnets only for quantities\
            on a sublattice!!!',type=int)
    parser.add_argument('--print-syms',action='store_const',const=True,default=False,help='Prints all symmetry operations.')
    parser.add_argument('--transform-result',action='store_const',const=True,default=False,help=
            'Keep the symmetry operations in findsym basis and transform the result to the correct basis.')
    parser.add_argument('--transform-syms',action='store_const',const=True,default=False,help=\
            'Transform the symmetry operations to the correct basis instead of transforming the result.')
    parser.add_argument('--syms',default=-1,help='Choose which symmetry operations to take, the rest is ignored. Insert symmetry operation\
             numbers separated by commas with no spaces. They are numbered as they appear in the findsym output file.\
             Also can include ranges. Example: 1-3,7,9-12',dest='syms_sel')
    parser.add_argument('--syms-noso',default=-1,help='Like --syms, but fot noso symmetry operations',dest='syms_sel_noso')
    parser.add_argument('--noso',action='store_const',const=True,default=False,help='Symmetry without spin-orbit coupling.')
    parser.add_argument('--ignore-op1eqop2',action='store_const',const=True,default=False,help='When op1=op2, the even part has to be \
            symmetric and the odd part antisymmetric. If this is selected, this property is ignored.',dest='ig_op1eqop2')
    if clargs != None:
        args = parser.parse_args(clargs.split())
    else:
        args = parser.parse_args()

    args_dict = vars(args)
    args_dict['debug'] = args_dict['debug'].split(',')

    #set the debug variables
    args_dict['debug_sym'] = False
    args_dict['debug_rename'] = False
    args_dict['debug_equiv'] = False
    args_dict['debug_tensor'] = False
    args_dict['debug_time'] = False
    args_dict['debug_symY'] = False
    args_dict['debug_noso'] = False
    args_dict['debug_op1eqop2'] = False

    if 'symmetrize' in args_dict['debug'] or 'all' in args_dict['debug'] or 'symmetrizeY' in args_dict['debug']:
        args_dict['debug_sym'] = True
    if 'rename' in args_dict['debug'] or 'all' in args_dict['debug']:
        args_dict['debug_rename'] = True
    if 'equiv' in args_dict['debug'] or 'all' in args_dict['debug']:
        args_dict['debug_equiv'] = True
    if 'exp' in args_dict['debug'] or 'all' in args_dict['debug']:
        args_dict['debug_tensor'] = True
    if 'time' in args_dict['debug'] or 'all' in args_dict['debug']:
        args_dict['debug_time'] = True
    if 'symmetrizeY' in args_dict['debug']:
        args_dict['debug_symY'] = True
    if 'noso' in args_dict['debug']:
        args_dict['debug_noso'] = True
    if 'op1eqop2' in args_dict['debug']:
        args_dict['debug_op1eqop2'] = True

    opt = options(args_dict)

    if (not opt['transform_result']) and (not opt['transform_syms']):
        if opt['op3']== None:
            opt['transform_result'] = True
        else:
            opt['transform_result'] = False

    return opt
