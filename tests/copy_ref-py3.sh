#!/bin/bash

files=$(cat file_list)

for f in ${files[@]};do
    #cp ${f}_ref ${f}_ref_old
    mv $f ${f}_ref_py3 
done
