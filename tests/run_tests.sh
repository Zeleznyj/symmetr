#!/bin/bash

set -e

set -x

py3=false
global=false
while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        --py3)
            py3=true
            ;;
        --global)
            global=true
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            exit 1
            ;;
    esac
    shift
done

if [ "$py3" = true ]; then
    symmetr=symmetr-py3
else
    symmetr=symmetr
fi

if [ "$global" = true ]; then
    path=''
else
    path='../../exec/'
fi

command=$path$symmetr

cd IrMn3
$command res j E -f findsym.in  > out_vv
$command res j E -f findsym.in -e  > out_vv_equiv
$command res s.j E -f findsym.in  > out_svv
$command mham -s 2,3 -f findsym.in  > out_mham_23
$command mham -s 2,2 -f findsym.in -e  > out_mham_22_equiv
cd ..

cd Mn2Au
$command res s E -f findsym.in -p 1  > out_sv0
$command res s E -f findsym.in -p 1 -p2 2  > out_sv0_1
$command res v E -f findsym.in  > out_vv
$command res s E -f findsym.in -p 1 -e  > out_sv0_equiv
cd ..

cd MnTe
$command res j E -f findsym.in > out_vv 
$command res s.v E -f findsym.in > out_svv
cd ..

cd NiMnSb
$command res s E --exp 1 -f NiMnSb.in_nonmag  > out_exp_1
$command res s E --exp 2 -f NiMnSb.in_nonmag  > out_exp_2
cd ..

cd groups
$command res j E -g P-43m --ignore-same-op-sym > out_jE_-43m
$command res s E -g P4mm --exp 1 > out_sE_4mm_exp1
cd ..
