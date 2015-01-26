#!/bin/bash

cd IrMn3
findsym < findsym.in | main.py v v > out_vv
cd ..

cd Mn2Au
findsym < findsym.in | main.py s v 0 > out_sv0
findsym < findsym.in | main.py v v > out_vv


