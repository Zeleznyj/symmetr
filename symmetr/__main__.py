#!/usr/bin/python

import funcs_main as fmain
import input

opt = input.parse()
    
if opt['print_syms']:
    syms = fmain.get_syms(opt)
    print 'Symmetry operations:'
    print 'Format: Number, space transformation, magnetic moment transformation, time-reversal, transformation of the sublattices'
    for sym in syms:
        print sym
    if opt['noso']:
        syms_noso = fmain.get_syms_noso(opt)
        print ''
        print 'Noso symmetry operations (in the magnetic basis)'
        for i,sym in enumerate(syms_noso):
            print i+1,sym

if opt['exp'] == -1:
    fmain.sym_linres(opt,printit=True)
else:
    fmain.sym_exp(opt,printit=True)
