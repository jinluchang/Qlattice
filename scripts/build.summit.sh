#!/bin/bash

./scripts/setenv-summit.sh

./scripts/fftw.sh
./scripts/fftwf.sh
./scripts/cuba.sh
./scripts/zlib.sh
./scripts/eigen.sh
./scripts/qlat.sh

pip3 install numpy
pip3 install simpy

./scripts/c-lime.sh
./scripts/mpfr.sh
./scripts/grid-summit.sh
./scripts/gpt.sh
