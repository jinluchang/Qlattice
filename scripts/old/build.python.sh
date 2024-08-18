#!/usr/bin/env bash

# Need to have gcc, openssl, zlib, gsl, openblas, perl, ffi installed beforehand.

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libmpfr-dev gcc perl libopenblas-dev libffi-dev libsqlite3-dev

set -e

{

./scripts/setenv.default.sh

./scripts/python.sh
./scripts/python-pip.sh
./scripts/re2c.sh
./scripts/ninja.sh
./scripts/ninja-script.sh
./scripts/python-meson.sh
./scripts/python-meson-py.sh
./scripts/python-packages.sh

} 2>&1 | tee $prefix/log.build.txt
