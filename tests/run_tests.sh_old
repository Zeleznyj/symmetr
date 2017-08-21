#!/bin/bash

set -e

set -x
cd IrMn3
 ../../main.py v v -f findsym.in > out_vv
../../main.py v v -f findsym.in -e > out_vv_equiv
cd ..

cd Mn2Au
../../main.py s v -f findsym.in -p 1 > out_sv0
../../main.py s v -f findsym.in -p 1 -p2 2 > out_sv0_1
../../main.py v v -f findsym.in > out_vv
../../main.py s v -f findsym.in -p 1 -e > out_sv0_equiv
cd ..

cd NiMnSb
../../main.py s v --exp 1 -f NiMnSb.in_nonmag > out_exp_1
../../main.py s v --exp 2 -f NiMnSb.in_nonmag > out_exp_2
cd ..
