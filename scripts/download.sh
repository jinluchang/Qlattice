#!/bin/bash

. conf.sh

mkdir -p $distfiles

cd $distfiles

wget -c "http://www.fftw.org/fftw-3.3.9.tar.gz"

wget -c "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"

wget -c "https://zlib.net/zlib-1.2.11.tar.gz"

wget -c "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2"

wget -O c-lime.tar.gz -c "https://github.com/usqcd-software/c-lime/tarball/master"

wget -c "https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.bz2"

if [ -d Grid ] ; then
    ( cd Grid ; git pull )
else
    git clone https://github.com/waterret/Grid.git
fi

(
cd Grid
pwd
git clean -f
./bootstrap.sh
)

if [ -d gpt ] ; then
    ( cd gpt ; git pull )
else
    git clone https://github.com/waterret/gpt.git
fi

sha256sum *.tar.* > sha256sums.txt

echo >> sha256sums.txt

echo -n "Grid: " >> sha256sums.txt
(
cd Grid
git rev-parse HEAD >> ../sha256sums.txt
)

echo -n "gpt: " >> sha256sums.txt
(
cd gpt
git rev-parse HEAD >> ../sha256sums.txt
)

cat sha256sums.txt
