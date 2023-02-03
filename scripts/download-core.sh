#!/bin/bash

source qcore/conf.sh

mkdir -p "$distfiles"

cd "$distfiles"

dget() {
    name="$1"
    url="$2"
    if [ -f "$name" ] ; then
        echo "$name is downloaded"
    else
        wget --no-check-certificate -O "$name" -c "$url"
    fi
}

aget() {
    url="$1"
    name="${1##*/}"
    dget "$name" "$url"
}

dget "c-lime.tar.gz" "https://github.com/usqcd-software/c-lime/tarball/master"

aget "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2"

aget "http://www.fftw.org/fftw-3.3.10.tar.gz"

aget "https://feynarts.de/cuba/Cuba-4.2.2.tar.gz"

aget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.bz2"

if [ -d Grid ] ; then
    ( cd Grid ; git pull )
else
    git clone https://github.com/paboyle/Grid.git Grid
fi

if [ -e Grid/configure ] ; then
    echo "Grid bootstrapped."
else
    ( cd Grid ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d Grid-clehner ] ; then
    ( cd Grid-clehner ; git pull https://github.com/jinluchang/Grid.git )
else
    git clone https://github.com/jinluchang/Grid.git Grid-clehner
fi

if [ -e Grid-clehner/configure ] ; then
    echo "Grid-clehner bootstrapped."
else
    ( cd Grid-clehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d gpt ] ; then
    ( cd gpt ; git pull https://github.com/jinluchang/gpt.git )
else
    git clone https://github.com/jinluchang/gpt.git
fi

if [ -d Hadrons ] ; then
    ( cd Hadrons ; git pull )
else
    git clone https://github.com/aportelli/Hadrons.git Hadrons
fi

if [ -e Hadrons/configure ] ; then
    echo "Hadrons bootstrapped."
else
    ( cd Hadrons ; git clean -f ; ./bootstrap.sh )
fi

if [ -d Grid-tblum ] ; then
    ( cd Grid-tblum ; git pull )
else
    git clone https://github.com/tblum2/Grid.git Grid-tblum
fi

if [ -e Grid-tblum/configure ] ; then
    echo "Grid-tblum bootstrapped."
else
    ( cd Grid-tblum ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d Hadrons-tblum ] ; then
    ( cd Hadrons-tblum ; git pull )
else
    git clone https://github.com/jinluchang/Hadrons.git Hadrons-tblum
fi

if [ -e Hadrons-tblum/configure ] ; then
    echo "Hadrons-tblum bootstrapped."
else
    ( cd Hadrons-tblum ; git clean -f ; ./bootstrap.sh )
fi

cd "$wd"
