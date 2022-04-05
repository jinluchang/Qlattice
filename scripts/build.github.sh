#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/cuba.sh
./scripts/eigen.sh
./scripts/qlat.sh

} |& tee $prefix/log.build.txt
