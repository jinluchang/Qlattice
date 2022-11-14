#!/bin/bash

. scripts/conf.sh

name=dist-update-hash

{

./scripts/update-sources.sh

cat $distfiles/sha256sums.txt

} 2>&1 | tee $prefix/log.$name.txt
