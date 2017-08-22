#!/bin/bash

set -e

files=()
files+=('IrMn3/out_vv')
files+=('IrMn3/out_vv_equiv')
files+=('IrMn3/out_svv')
files+=('Mn2Au/out_sv0')
files+=('Mn2Au/out_sv0_1')
files+=('Mn2Au/out_vv')
files+=('Mn2Au/out_sv0_equiv')
files+=('NiMnSb/out_exp_1')
files+=('NiMnSb/out_exp_2')
files+=('NiMnSb/out_exp_1')
files+=('NiMnSb/out_exp_2')
files+=('groups/out_jE_-43m')
files+=('groups/out_sE_4mm_exp1')


for f in ${files[@]};do 
    if diff $f ${f}_ref > /dev/null
    then
        echo $f ok!
    else
        echo $f differs from ${f}_ref
        vimdiff $f ${f}_ref
    fi
done
