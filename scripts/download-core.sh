#!/bin/bash

. scripts/conf.sh

mkdir -p $distfiles

cd $distfiles

dget() {
    name="$1"
    url="$2"
    if [ -f "$name" ] ; then
        echo "$name is downloaded"
    else
        wget --no-check-certificate -O "$name" -c "$url"
    fi
}

dget "Cuba-4.2.1.tar.gz"  "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"

dget "eigen-3.3.9.tar.bz2" "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2"

dget "gsl-2.7.1.tar.gz" "https://gnu.askapache.com/gsl/gsl-2.7.1.tar.gz"

dget "zlib-1.2.11.tar.gz" "https://versaweb.dl.sourceforge.net/project/libpng/zlib/1.2.11/zlib-1.2.11.tar.gz"

dget "fftw-3.3.10.tar.gz" "http://www.fftw.org/fftw-3.3.10.tar.gz"

dget "c-lime.tar.gz" "https://github.com/usqcd-software/c-lime/tarball/master"

if [ -d Grid-lehner ] ; then
    ( cd Grid-lehner ; git pull https://github.com/waterret/Grid.git )
else
    git clone https://github.com/waterret/Grid.git Grid-lehner
fi

if [ -e Grid-lehner/configure ] ; then
    echo "Grid-lehner bootstrapped."
else
    ( cd Grid-lehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d gpt ] ; then
    ( cd gpt ; git pull https://github.com/waterret/gpt.git )
else
    git clone https://github.com/waterret/gpt.git
fi

cd $wd
