#!/bin/bash

set -e

{

./scripts/setenv.default.sh

./scripts/perl.sh
./scripts/openssl.sh
./scripts/libffi.sh
./scripts/openblas.sh
./scripts/python.sh
./scripts/python-meson.sh
./scripts/python-packages.sh

} 2>&1 | tee $prefix/log.build.txt
