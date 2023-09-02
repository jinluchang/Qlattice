#!/bin/bash

set -e

# Ubuntu packages: libgsl-dev zlib1g-dev libssl-dev libopenmpi-dev ninja-build libsqlite3-dev libgdbm-dev liblzma-dev libbz2-dev uuid-dev tk-dev libgeos-dev libopenblas-dev node-configurable-http-proxy

{

    ./scripts/setenv.default.sh
    ./scripts/qcore.sh

    ./scripts/python-venv-empty.sh
    ./scripts/python-pip.sh
    ./scripts/ninja-script.sh
    ./scripts/python-meson.sh
    ./scripts/python-meson-py.sh
    ./scripts/python-packages.sh

} 2>&1 | tee $prefix/log.build.txt
