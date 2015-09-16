#!/bin/bash

set -e

set -x
cd IrMn3
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py v v -e nonmag_out > out_vv_equiv
cd ..

cd Mn2Au
findsym < findsym.in | ../../main.py s v -p 1 > out_sv0
findsym < findsym.in | ../../main.py s v -p 1 -p2 2 > out_sv0_1
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py s v -p 1 -e nonmag_out > out_sv0_equiv
cd ..

cd NiMnSb
findsym < NiMnSb.in_nonmag | ../../expansions.py s v 1 > out_exp_1
findsym < NiMnSb.in_nonmag | ../../expansions.py s v 2 > out_exp_2
cd ..
