#!/usr/bin/python

import funcs_main as fmain
import symT
import input

opt = input.parse()
    
if opt['print_syms']:
    syms = symT.get_syms(opt)
    print 'Symmetry operations:'
    print 'Format: Number, space transformation, magnetic moment transformation, time-reversal, transformation of the sublattices'
    for sym in syms:
        print sym
    if opt['mode'] == 'res' and opt['noso']:
        syms_noso = fmain.get_syms_noso(opt)
        print ''
        print 'Noso symmetry operations (in the magnetic basis)'
        for i,sym in enumerate(syms_noso):
            print i+1,sym

if opt['print_opt']:
    print opt

if opt['mode'] == 'res':
    fmain.sym_res(opt,printit=True)
else:
    fmain.sym_mham(opt,printit=True) 
