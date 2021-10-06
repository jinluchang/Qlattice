#!/bin/bash

. conf.sh

mkdir -p $distfiles

cd $distfiles

if [ -f "fftw-3.3.9.tar.gz" ] ; then
    echo "fftw is downloaded"
else
    wget -c "http://www.fftw.org/fftw-3.3.9.tar.gz"
fi

if [ -f "Cuba-4.2.1.tar.gz" ] ; then
    echo "cuba is downloaded"
else
    wget -c "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"
fi

if [ -f "zlib-1.2.11.tar.gz" ] ; then
    echo "zlib is downloaded"
else
    wget -c "https://zlib.net/zlib-1.2.11.tar.gz"
fi

if [ -f "eigen-3.3.9.tar.bz2" ] ; then
    echo "eigen is downloaded"
else
    wget -c "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2"
fi

if [ -f "c-lime.tar.gz" ] ; then
    echo "c-lime is downloaded"
else
    wget -O c-lime.tar.gz -c "https://github.com/usqcd-software/c-lime/tarball/master"
fi

if [ -f "mpfr-4.1.0.tar.bz2" ] ; then
    echo "mpfr is downloaded"
else
    wget -c "https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.bz2"
fi

if [ -f "hdf5-1.12.1.tar.bz2" ] ; then
    echo "hdf5 is downloaded"
else
    wget -c "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.bz2"
fi

if [ -d Grid ] ; then
    ( cd Grid ; git pull )
else
    git clone https://github.com/paboyle/Grid.git Grid
fi

if [ -e Grid/configure ] ; then
    echo "Grid bootstrapped."
else
    ( cd Grid ; git clean -f ; ./bootstrap.sh )
fi

if [ -d Grid-lehner ] ; then
    ( cd Grid-lehner ; git pull )
else
    git clone https://github.com/lehner/Grid.git Grid-lehner
fi

if [ -e Grid-lehner/configure ] ; then
    echo "Grid-lehner bootstrapped."
else
    ( cd Grid-lehner ; git clean -f ; ./bootstrap.sh )
fi

if [ -d gpt ] ; then
    ( cd gpt ; git pull )
else
    git clone https://github.com/lehner/gpt.git
fi

if [ -d Grid-tblum ] ; then
    ( cd Grid-tblum ; git pull )
else
    git clone https://github.com/tblum2/Grid.git Grid-tblum
fi

if [ -e Grid-tblum/configure ] ; then
    echo "Grid-tblum bootstrapped."
else
    ( cd Grid-tblum ; git clean -f ; ./bootstrap.sh )
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

cd $wd

./scripts/dist-update-hash.sh
