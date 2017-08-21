#!/bin/bash

set -e

files=()
files+=('IrMn3/out_vv')
files+=('IrMn3/out_vv_equiv')
files+=('IrMn3/out_mham_23')
files+=('IrMn3/out_mham_22_equiv')
files+=('Mn2Au/out_sv0')
files+=('Mn2Au/out_sv0_1')
files+=('Mn2Au/out_vv')
files+=('Mn2Au/out_sv0_equiv')
files+=('NiMnSb/out_exp_1')
files+=('NiMnSb/out_exp_2')

for f in ${files[@]};do
    cp $f ${f}_ref
done
