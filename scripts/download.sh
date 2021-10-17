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

dget "binutils-2.37.tar.xz" "https://ftp.gnu.org/gnu/binutils/binutils-2.37.tar.xz"

dget "openmpi-4.1.1.tar.bz2" "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.1.tar.bz2"

dget "llvm-project-13.0.0.src.tar.xz" "https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/llvm-project-13.0.0.src.tar.xz"

dget "perl-5.34.0.tar.gz" "https://www.cpan.org/src/5.0/perl-5.34.0.tar.gz"

dget "openssl-3.0.0.tar.gz" "https://www.openssl.org/source/openssl-3.0.0.tar.gz"

dget "Python-3.10.0.tar.xz" "https://www.python.org/ftp/python/3.10.0/Python-3.10.0.tar.xz"

dget "libffi-3.4.2.tar.gz" "https://github.com/libffi/libffi/releases/download/v3.4.2/libffi-3.4.2.tar.gz"

dget "re2c-2.2.tar.xz" "https://github.com/skvadrik/re2c/releases/download/2.2/re2c-2.2.tar.xz"

dget "ninja-1.10.2.tar.gz" "https://github.com/ninja-build/ninja/archive/refs/tags/v1.10.2.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "numpy-1.21.2.zip" "https://files.pythonhosted.org/packages/3a/be/650f9c091ef71cb01d735775d554e068752d3ff63d7943b26316dc401749/numpy-1.21.2.zip"
dget "Cython-0.29.24.tar.gz" "https://files.pythonhosted.org/packages/59/e3/78c921adf4423fff68da327cc91b73a16c63f29752efe7beb6b88b6dd79d/Cython-0.29.24.tar.gz"
dget "setuptools-49.1.3.zip" "https://files.pythonhosted.org/packages/d0/4a/22ee76842d8ffc123d4fc48d24a623c1d206b99968fe3960039f1efc2cbc/setuptools-49.1.3.zip"
dget "wheel-0.36.2.tar.gz" "https://files.pythonhosted.org/packages/ed/46/e298a50dde405e1c202e316fa6a3015ff9288423661d7ea5e8f22f589071/wheel-0.36.2.tar.gz"
dget "sympy-1.9.tar.gz" "https://files.pythonhosted.org/packages/26/86/902ee78db1bab1f0410f799869a49bb03b83be8d44c23b224d9db34f21c3/sympy-1.9.tar.gz"
dget "mpmath-1.2.1.tar.gz" "https://files.pythonhosted.org/packages/95/ba/7384cb4db4ed474d4582944053549e02ec25da630810e4a23454bc9fa617/mpmath-1.2.1.tar.gz"
dget "mpmath-1.2.1-py3-none-any.whl" "https://files.pythonhosted.org/packages/d4/cf/3965bddbb4f1a61c49aacae0e78fd1fe36b5dc36c797b31f30cf07dcbbb7/mpmath-1.2.1-py3-none-any.whl"
)


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
