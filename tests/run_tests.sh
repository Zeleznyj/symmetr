#!/bin/bash

set -e

set -x
cd IrMn3
../../exec/symmetr res j E -f findsym.in  > out_vv
../../exec/symmetr res j E -f findsym.in -e  > out_vv_equiv
../../exec/symmetr mham -s 2,3 -f findsym.in  > out_mham_23
../../exec/symmetr mham -s 2,2 -f findsym.in -e  > out_mham_22_equiv
cd ..

cd Mn2Au
../../exec/symmetr res s E -f findsym.in -p 1  > out_sv0
../../exec/symmetr res s E -f findsym.in -p 1 -p2 2  > out_sv0_1
../../exec/symmetr res v E -f findsym.in  > out_vv
../../exec/symmetr res s E -f findsym.in -p 1 -e  > out_sv0_equiv
cd ..

cd NiMnSb
../../exec/symmetr res s E --exp 1 -f NiMnSb.in_nonmag  > out_exp_1
../../exec/symmetr res s E --exp 2 -f NiMnSb.in_nonmag  > out_exp_2
cd ..

cd groups
../../exec/symmetr res j E -g P-43m --ignore-op1eqop2 > out_jE_-43m
../../exec/symmetr res s E -g P4mm --exp 1 > out_sE_4mm_exp1
cd ..
