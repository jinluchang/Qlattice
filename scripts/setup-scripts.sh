#!/bin/bash

. scripts/res/conf.sh

name=setup-scripts

{

echo "!!!! build $name !!!!"

mkdir -p "$prefix/bin"

cp -rpv scripts/res/organize-colon-list.py "$prefix"/bin/

./scripts/compiler-wrappers.sh

echo "!!!! $name build !!!!"

} 2>&1 | tee $prefix/log.$name.txt
