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

dget "xz-5.2.5.tar.gz" "https://tukaani.org/xz/xz-5.2.5.tar.gz"

dget "tar-1.34.tar.gz" "https://ftp.gnu.org/gnu/tar/tar-1.34.tar.gz"

dget "fftw-3.3.10.tar.gz" "http://www.fftw.org/fftw-3.3.10.tar.gz"

dget "Cuba-4.2.1.tar.gz"  "http://www.feynarts.de/cuba/Cuba-4.2.1.tar.gz"

dget "gsl-2.7.1.tar.gz" "https://gnu.askapache.com/gsl/gsl-2.7.1.tar.gz"

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

dget "OpenBLAS-0.3.18.tar.gz" "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.18/OpenBLAS-0.3.18.tar.gz"

dget "lapack-3.10.0.tar.gz" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz"

dget "zlib-1.2.11.tar.gz" "https://versaweb.dl.sourceforge.net/project/libpng/zlib/1.2.11/zlib-1.2.11.tar.gz"

dget "gnuplot-5.4.3.tar.gz" "https://versaweb.dl.sourceforge.net/project/gnuplot/gnuplot/5.4.3/gnuplot-5.4.3.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "pip-22.2.2.tar.gz" "https://files.pythonhosted.org/packages/4b/30/e15b806597e67057e07a5acdc135216ccbf76a5f1681a324533b61066b0b/pip-22.2.2.tar.gz"
dget "ninja-1.10.2.3.tar.gz" "https://files.pythonhosted.org/packages/00/99/5beedbf09e3ec6b617606df42d04c4251959caddbd98397cce21da4c52d1/ninja-1.10.2.3.tar.gz"
dget "meson-0.63.2.tar.gz" "https://files.pythonhosted.org/packages/a7/f0/565f731cd138a516c2dba8439e47c5622493c82f41c4845d287617ef6ec9/meson-0.63.2.tar.gz"
dget "scikit-build-0.15.0.tar.gz" "https://files.pythonhosted.org/packages/9e/e2/2e440c30e93fc5b505ee56169a4396b05e797a1daadb721aba429adbfd51/scikit-build-0.15.0.tar.gz"
dget "build-0.8.0.tar.gz" "https://files.pythonhosted.org/packages/52/fa/931038182be739955cf83179d9b9a6ce9832bc5f9a917a006f765cb53a1f/build-0.8.0.tar.gz"
dget "pyproject-metadata-0.6.1.tar.gz" "https://files.pythonhosted.org/packages/8e/e3/4314ed332414c059aed340fe2097d5ca065cc8d165c08554c83eb96a59c4/pyproject-metadata-0.6.1.tar.gz"
dget "setuptools_scm-7.0.5.tar.gz" "https://files.pythonhosted.org/packages/d0/43/f038b5009f93bcd77b1b8da9e6d424b739ab17aec9726f3a99eba23d53ca/setuptools_scm-7.0.5.tar.gz"
dget "packaging-21.3.tar.gz" "https://files.pythonhosted.org/packages/df/9e/d1a7217f69310c1db8fdf8ab396229f55a699ce34a203691794c5d1cad0c/packaging-21.3.tar.gz"
dget "typing_extensions-4.3.0.tar.gz" "https://files.pythonhosted.org/packages/9e/1d/d128169ff58c501059330f1ad96ed62b79114a2eb30b8238af63a2e27f70/typing_extensions-4.3.0.tar.gz"
dget "pyparsing-3.0.9.tar.gz" "https://files.pythonhosted.org/packages/71/22/207523d16464c40a0310d2d4d8926daffa00ac1f5b1576170a32db749636/pyparsing-3.0.9.tar.gz"
dget "distro-1.7.0.tar.gz" "https://files.pythonhosted.org/packages/b5/7e/ddfbd640ac9a82e60718558a3de7d5988a7d4648385cf00318f60a8b073a/distro-1.7.0.tar.gz"
dget "flit_core-3.7.1.tar.gz" "https://files.pythonhosted.org/packages/15/d1/d8798b83e953fd6f86ca9b50f93eec464a9305b0661469c8234e61095481/flit_core-3.7.1.tar.gz"
dget "tomli-2.0.1.tar.gz" "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
dget "meson_python-0.8.1.tar.gz" "https://files.pythonhosted.org/packages/c3/94/957a04750188722d09ade6ae3731b115366177faee32ace175c3ca59358b/meson_python-0.8.1.tar.gz"
dget "numpy-1.21.4.zip" "https://files.pythonhosted.org/packages/fb/48/b0708ebd7718a8933f0d3937513ef8ef2f4f04529f1f66ca86d873043921/numpy-1.21.4.zip"
dget "Cython-0.29.24.tar.gz" "https://files.pythonhosted.org/packages/59/e3/78c921adf4423fff68da327cc91b73a16c63f29752efe7beb6b88b6dd79d/Cython-0.29.24.tar.gz"
dget "setuptools-49.1.3.zip" "https://files.pythonhosted.org/packages/d0/4a/22ee76842d8ffc123d4fc48d24a623c1d206b99968fe3960039f1efc2cbc/setuptools-49.1.3.zip"
dget "wheel-0.36.2.tar.gz" "https://files.pythonhosted.org/packages/ed/46/e298a50dde405e1c202e316fa6a3015ff9288423661d7ea5e8f22f589071/wheel-0.36.2.tar.gz"
dget "sympy-1.9.tar.gz" "https://files.pythonhosted.org/packages/26/86/902ee78db1bab1f0410f799869a49bb03b83be8d44c23b224d9db34f21c3/sympy-1.9.tar.gz"
dget "mpmath-1.2.1.tar.gz" "https://files.pythonhosted.org/packages/95/ba/7384cb4db4ed474d4582944053549e02ec25da630810e4a23454bc9fa617/mpmath-1.2.1.tar.gz"
dget "mpmath-1.2.1-py3-none-any.whl" "https://files.pythonhosted.org/packages/d4/cf/3965bddbb4f1a61c49aacae0e78fd1fe36b5dc36c797b31f30cf07dcbbb7/mpmath-1.2.1-py3-none-any.whl"
dget "scipy-1.7.3.tar.gz" "https://files.pythonhosted.org/packages/61/67/1a654b96309c991762ee9bc39c363fc618076b155fe52d295211cf2536c7/scipy-1.7.3.tar.gz"
dget "pybind11-2.7.1.tar.gz" "https://files.pythonhosted.org/packages/66/99/fc60e2287bb2309b8db4d0f080770ecc8d37dc64911e37b86698ec4b6a51/pybind11-2.7.1.tar.gz"
dget "pythran-0.9.12.post1.tar.gz" "https://files.pythonhosted.org/packages/c6/e6/986a967dcca91d89e36f4d4a2f69a052030bce01a7cd48a6b7fba1a50189/pythran-0.9.12.post1.tar.gz"
dget "pythran-0.9.12.post1-py3-none-any.whl" "https://files.pythonhosted.org/packages/e1/0a/a1e28f3c532655032a73482cb11ae63f8c38c8f58dc1a713b3791855fa68/pythran-0.9.12.post1-py3-none-any.whl"
dget "ply-3.11.tar.gz" "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
dget "gast-0.5.0.tar.gz" "https://files.pythonhosted.org/packages/a6/fb/7ff6a4ee66673c5964d3cf515ae85ba2076bc64bc2dcbbbd0153718b005f/gast-0.5.0.tar.gz"
dget "beniget-0.4.0.tar.gz" "https://files.pythonhosted.org/packages/36/09/a4a6a967ca5bcfc0bd6162df4ee93017301fa7d9671483c849300bdba0db/beniget-0.4.0.tar.gz"
dget "psutil-5.9.0.tar.gz" "https://files.pythonhosted.org/packages/47/b6/ea8a7728f096a597f0032564e8013b705aa992a0990becd773dcc4d7b4a7/psutil-5.9.0.tar.gz"
)

if [ -d Grid-paboyle ] ; then
    ( cd Grid-paboyle ; git pull )
else
    git clone https://github.com/paboyle/Grid.git Grid-paboyle
fi

if [ -e Grid-paboyle/configure ] ; then
    echo "Grid bootstrapped."
else
    ( cd Grid-paboyle ; git clean -f ; ./bootstrap.sh )
fi

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
