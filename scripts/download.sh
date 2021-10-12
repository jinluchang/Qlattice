#!/bin/bash

. conf.sh

mkdir -p $distfiles

cd $distfiles

dget() {
    name="$1"
    url="$2"
    if [ -f "$name" ] ; then
        echo "$name is downloaded"
    else
        wget -O "$name" -c "$url"
    fi
}

dget "xz-5.2.5.tar.gz" "https://tukaani.org/xz/xz-5.2.5.tar.gz"

dget "tar-1.34.tar.gz" "https://ftp.gnu.org/gnu/tar/tar-1.34.tar.gz"

dget "fftw-3.3.10.tar.gz" "http://www.fftw.org/fftw-3.3.10.tar.gz"

dget "Cuba-4.2.1.tar.gz"  "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"

dget "zlib-1.2.11.tar.gz" "https://zlib.net/zlib-1.2.11.tar.gz"

dget "eigen-3.3.9.tar.bz2" "https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.bz2"

dget "c-lime.tar.gz" "https://github.com/usqcd-software/c-lime/tarball/master"

dget "gmp-6.2.1.tar.bz2" "https://gmplib.org/download/gmp/gmp-6.2.1.tar.bz2"

dget "mpfr-4.1.0.tar.bz2" "https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.bz2"

dget "mpc-1.2.1.tar.gz" "https://ftp.gnu.org/gnu/mpc/mpc-1.2.1.tar.gz"

dget "autoconf-2.71.tar.gz" "http://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz"

dget "automake-1.16.5.tar.gz" "http://ftp.gnu.org/gnu/automake/automake-1.16.5.tar.gz"

dget "cmake-3.21.3.tar.gz" "https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3.tar.gz"

dget "hdf5-1.10.7.tar.bz2" "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.bz2"

dget "gcc-11.2.0.tar.xz" "http://mirrors.concertpass.com/gcc/releases/gcc-11.2.0/gcc-11.2.0.tar.xz"

dget "openmpi-4.1.1.tar.bz2" "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.bz2"

dget "llvm-project-13.0.0.src.tar.xz" "https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/llvm-project-13.0.0.src.tar.xz"

dget "openssl-3.0.0.tar.gz" "https://www.openssl.org/source/openssl-3.0.0.tar.gz"

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
