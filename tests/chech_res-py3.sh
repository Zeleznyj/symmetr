#!/bin/bash

set -e

files=$(cat file_list)

for f in ${files[@]};do 
    if diff $f ${f}_ref_py3 > /dev/null
    then
        echo $f ok!
    else
        echo $f differs from ${f}_ref
        vimdiff $f ${f}_ref_py3
    fi
done
