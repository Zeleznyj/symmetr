#!/bin/bash

set -e

set -x
cd IrMn3
symmetr3 res j E -f findsym.in  > out_vv
symmetr3 res j E -f findsym.in -e  > out_vv_equiv
symmetr3 res s.j E -f findsym.in  > out_svv
symmetr3 mham -s 2,3 -f findsym.in  > out_mham_23
symmetr3 mham -s 2,2 -f findsym.in -e  > out_mham_22_equiv
cd ..

cd Mn2Au
symmetr3 res s E -f findsym.in -p 1  > out_sv0
symmetr3 res s E -f findsym.in -p 1 -p2 2  > out_sv0_1
symmetr3 res v E -f findsym.in  > out_vv
symmetr3 res s E -f findsym.in -p 1 -e  > out_sv0_equiv
cd ..

cd NiMnSb
symmetr3 res s E --exp 1 -f NiMnSb.in_nonmag  > out_exp_1
symmetr3 res s E --exp 2 -f NiMnSb.in_nonmag  > out_exp_2
cd ..

cd groups
symmetr3 res j E -g P-43m --ignore-same-op-sym > out_jE_-43m
symmetr3 res s E -g P4mm --exp 1 > out_sE_4mm_exp1
cd ..
