#!/bin/bash

cd IrMn3
set -x
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py v v -e nonmag_out > out_vv_equiv
set +x
cd ..

cd Mn2Au
set -x
findsym < findsym.in | ../../main.py s v -p 1 > out_sv0
findsym < findsym.in | ../../main.py s v -p 1 -p2 2 > out_sv0_1
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py s v -p 1 -e nonmag_out > out_sv0_equiv
set +x
cd ..


