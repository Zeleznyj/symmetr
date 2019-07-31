from __future__ import print_function
from __future__ import absolute_import
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#!/usr/bin/env python

import sys
from . import funcs_main as fmain
from . import symT
from . import input
from . import fslib
from .version import __version__

if '--version' in sys.argv:
    print('Symmetr version: {}'.format(__version__))
    print('Python version: {}'.format(sys.version))

opt = input.parse()

if opt['print_syms']:
    syms = symT.get_syms(opt)
    print('Symmetry operations:')
    print('Format: Number, space transformation, magnetic moment transformation, time-reversal, transformation of the sublattices')
    for i,sym in enumerate(syms):
        print('Symmetry: ', i)
        print(sym)
    if opt['mode'] == 'res' and opt['noso']:
        syms_noso = symT.get_syms_noso(opt)
        print('')
        print('Noso symmetry operations (in the magnetic basis)')
        for i,sym in enumerate(syms_noso):
            print(i+1,sym)

if opt['print_opt']:
    print(opt)

if opt['print_pos']:
    lines = fslib.run_fs(opt['inp'])
    pos = fslib.r_pos(lines)
    print('Atomic sites as used by findsym:')
    for p in pos:
        print(p[-1],p[0:3])
    print('')

if opt['mode'] == 'res':
    fmain.sym_res(opt,printit=True)
else:
    fmain.sym_mham(opt,printit=True) 
