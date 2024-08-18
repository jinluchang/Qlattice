#!/usr/bin/env bash

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

aget "http://usqcd-software.github.io/downloads/c-lime/lime-1.3.2.tar.gz"

aget "http://usqcd-software.github.io/downloads/qmp/qmp-2.5.4.tar.gz"

aget "http://usqcd-software.github.io/downloads/qio/qio-3.0.0.tar.gz"

aget "https://downloads.sourceforge.net/libtirpc/libtirpc-1.3.4.tar.bz2"

aget "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2"

aget "https://fftw.org/pub/fftw/fftw-3.3.10.tar.gz"

# aget "https://feynarts.de/cuba/Cuba-4.2.2.tar.gz"
aget "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/Cuba-4.2.2.tar.gz"

# aget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.bz2"
aget "https://github.com/jinluchang/Qlattice-distfiles/raw/main/distfiles/hdf5-1.14.3.tar.bz2"

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
    ( cd Grid-clehner ; git pull )
else
    git clone https://github.com/lehner/Grid.git Grid-clehner
fi

if [ -e Grid-clehner/configure ] ; then
    echo "Grid-clehner bootstrapped."
else
    ( cd Grid-clehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
fi

if [ -d gpt ] ; then
    ( cd gpt ; git pull )
else
    git clone https://github.com/lehner/gpt.git
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
    git clone https://github.com/tblum2/Hadrons.git Hadrons-tblum
fi

if [ -e Hadrons-tblum/configure ] ; then
    echo "Hadrons-tblum bootstrapped."
else
    ( cd Hadrons-tblum ; git clean -f ; ./bootstrap.sh )
fi

if [ -d CPS ] ; then
    ( cd CPS ; git pull )
else
    git clone https://github.com/RBC-UKQCD/CPS_public.git CPS
fi

cd "$wd"
