#!/bin/bash

files=()
files+=('IrMn3/out_vv')
files+=('IrMn3/out_vv_equiv')
files+=('IrMn3/out_svv')
files+=('IrMn3/out_mham_23')
files+=('IrMn3/out_mham_22_equiv')
files+=('Mn2Au/out_sv0')
files+=('Mn2Au/out_sv0_1')
files+=('Mn2Au/out_vv')
files+=('Mn2Au/out_sv0_equiv')
files+=('NiMnSb/out_exp_1')
files+=('NiMnSb/out_exp_2')
files+=('groups/out_jE_-43m')
files+=('groups/out_sE_4mm_exp1')

for f in ${files[@]};do
    #cp ${f}_ref ${f}_ref_old
    mv $f ${f}_ref_py3 
done
