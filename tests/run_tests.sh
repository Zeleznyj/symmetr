#!/bin/bash

cd IrMn3
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py v v -e nonmag_out > out_vv_equiv
cd ..

cd Mn2Au
findsym < findsym.in | ../../main.py s v -p 0 > out_sv0
findsym < findsym.in | ../../main.py s v -p 0 -p2 1 > out_sv0_1
findsym < findsym.in | ../../main.py v v > out_vv
findsym < findsym.in | ../../main.py s v -p 0 -e nonmag_out > out_sv0_equiv


