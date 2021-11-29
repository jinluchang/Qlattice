#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-packages.sh

} |& tee $prefix/log.build.txt
