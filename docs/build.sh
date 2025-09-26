#!/bin/bash

rm source/api/*

cd ..
#sphinx-apidoc -o docs/source/api symmetr symmetr/hsnf symmetr/test.py symmetr/noso.py --force --no-toc --separate --templatedir docs/apidocs_templates
sphinx-apidoc -o docs/source/api symmetr \
   symmetr/hsnf \
   symmetr/test.py \
   symmetr/noso.py \
   symmetr/groups.py \
   symmetr/rename.py \
   symmetr/version.py \
   symmetr/conv_index.py \
   --force \
   --separate
#   --templatedir docs/apidocs_templates
cd docs
rm source/api/modules.rst
make html
