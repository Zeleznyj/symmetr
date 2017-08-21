#!/bin/bash

set -e

set -x
cd IrMn3
python ../../symmetr res j E -f findsym.in > out_vv
python ../../symmetr res j E -f findsym.in -e > out_vv_equiv
python ../../symmetr mham -s 2,3 -f findsym.in > out_mham_23
python ../../symmetr mham -s 2,2 -f findsym.in -e > out_mham_22_equiv
cd ..

cd Mn2Au
python ../../symmetr res s E -f findsym.in -p 1 > out_sv0
python ../../symmetr res s E -f findsym.in -p 1 -p2 2 > out_sv0_1
python ../../symmetr res v E -f findsym.in > out_vv
python ../../symmetr res s E -f findsym.in -p 1 -e > out_sv0_equiv
cd ..

cd NiMnSb
python ../../symmetr res s E --exp 1 -f NiMnSb.in_nonmag > out_exp_1
python ../../symmetr res s E --exp 2 -f NiMnSb.in_nonmag > out_exp_2
cd ..
