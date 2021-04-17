#!/bin/bash

. conf.sh

cd $distfiles

wget -c "http://www.fftw.org/fftw-3.3.9.tar.gz"

wget -c "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"

wget -c "https://zlib.net/zlib-1.2.11.tar.gz"

wget -c "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2"

wget -O c-lime.tar.gz -c "https://github.com/usqcd-software/c-lime/tarball/master"

wget -O Grid-origin.tar.gz -c "https://github.com/waterret/Grid/tarball/feature/gpt"

(
rm -rf waterret-Grid-*
tar xaf Grid-origin.tar.gz
cd waterret-Grid-*
pwd
./bootstrap.sh
)

wget -O gpt.tar.gz -c "https://github.com/waterret/gpt/tarball/master"

sha256sum *.tar.* > sha256sums.txt
