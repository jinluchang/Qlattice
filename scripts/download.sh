#!/bin/bash

. scripts/conf.sh

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

dget "ninja-1.11.1.tar.gz" "https://github.com/ninja-build/ninja/archive/refs/tags/v1.11.1.zip"

dget "lapack-3.10.1.tar.gz" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz"

aget "https://tukaani.org/xz/xz-5.2.7.tar.gz"

aget "https://ftp.gnu.org/gnu/tar/tar-1.34.tar.gz"

aget "http://www.fftw.org/fftw-3.3.10.tar.gz"

aget "https://feynarts.de/cuba/Cuba-4.2.2.tar.gz"

aget "https://gnu.askapache.com/gsl/gsl-2.7.1.tar.gz"

aget "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2"

aget "https://gmplib.org/download/gmp/gmp-6.2.1.tar.bz2"

aget "https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.bz2"

aget "https://ftp.gnu.org/gnu/mpc/mpc-1.2.1.tar.gz"

aget "http://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz"

aget "http://ftp.gnu.org/gnu/automake/automake-1.16.5.tar.gz"

aget "https://github.com/Kitware/CMake/releases/download/v3.24.3/cmake-3.24.3.tar.gz"

aget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.bz2"

aget "http://mirrors.concertpass.com/gcc/releases/gcc-12.2.0/gcc-12.2.0.tar.xz"

aget "https://ftp.gnu.org/gnu/binutils/binutils-2.39.tar.xz"

aget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.bz2"

aget "https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.2/llvm-project-15.0.2.src.tar.xz"

aget "https://www.cpan.org/src/5.0/perl-5.34.0.tar.gz"

aget "https://www.openssl.org/source/openssl-3.0.7.tar.gz"

aget "https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tar.xz"

aget "https://github.com/libffi/libffi/releases/download/v3.4.4/libffi-3.4.4.tar.gz"

aget "https://github.com/skvadrik/re2c/releases/download/2.2/re2c-2.2.tar.xz"

aget "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.21/OpenBLAS-0.3.21.tar.gz"

aget "https://zlib.net/zlib-1.2.13.tar.gz"

aget "https://gigenet.dl.sourceforge.net/project/gnuplot/gnuplot/5.4.5/gnuplot-5.4.5.tar.gz"

aget "https://cytranet.dl.sourceforge.net/project/tclap/tclap-1.2.5.tar.gz"

(
mkdir -p python-packages
cd python-packages
aget "https://files.pythonhosted.org/packages/a3/50/c4d2727b99052780aad92c7297465af5fe6eec2dbae490aa9763273ffdc1/pip-22.3.1.tar.gz"
aget "https://files.pythonhosted.org/packages/f7/69/938374c8ebfeda683863b22e936f5d465ac9f5bf42be238504c018123190/ninja-1.11.1.tar.gz"
aget "https://files.pythonhosted.org/packages/69/78/ed2b5d1d88ab5f86f93f72a4f5d86989dcea939c274917d53acdd534056b/meson-0.64.0.tar.gz"
aget "https://files.pythonhosted.org/packages/9e/e2/2e440c30e93fc5b505ee56169a4396b05e797a1daadb721aba429adbfd51/scikit-build-0.15.0.tar.gz"
aget "https://files.pythonhosted.org/packages/52/fa/931038182be739955cf83179d9b9a6ce9832bc5f9a917a006f765cb53a1f/build-0.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/8e/e3/4314ed332414c059aed340fe2097d5ca065cc8d165c08554c83eb96a59c4/pyproject-metadata-0.6.1.tar.gz"
aget "https://files.pythonhosted.org/packages/d0/43/f038b5009f93bcd77b1b8da9e6d424b739ab17aec9726f3a99eba23d53ca/setuptools_scm-7.0.5.tar.gz"
aget "https://files.pythonhosted.org/packages/df/9e/d1a7217f69310c1db8fdf8ab396229f55a699ce34a203691794c5d1cad0c/packaging-21.3.tar.gz"
aget "https://files.pythonhosted.org/packages/9e/1d/d128169ff58c501059330f1ad96ed62b79114a2eb30b8238af63a2e27f70/typing_extensions-4.3.0.tar.gz"
aget "https://files.pythonhosted.org/packages/71/22/207523d16464c40a0310d2d4d8926daffa00ac1f5b1576170a32db749636/pyparsing-3.0.9.tar.gz"
aget "https://files.pythonhosted.org/packages/b5/7e/ddfbd640ac9a82e60718558a3de7d5988a7d4648385cf00318f60a8b073a/distro-1.7.0.tar.gz"
aget "https://files.pythonhosted.org/packages/15/d1/d8798b83e953fd6f86ca9b50f93eec464a9305b0661469c8234e61095481/flit_core-3.7.1.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
aget "https://files.pythonhosted.org/packages/ee/96/5143db2524a4400db85f9ab5c11a1cc853b6770c316ecc81798877144c15/meson_python-0.10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/6a/fa/5ec0fa9095c9b72cb1c31a8175c4c6745bf5927d1045d7a70df35d54944f/setuptools-59.6.0.tar.gz"
aget "https://files.pythonhosted.org/packages/ed/46/e298a50dde405e1c202e316fa6a3015ff9288423661d7ea5e8f22f589071/wheel-0.36.2.tar.gz"
aget "https://files.pythonhosted.org/packages/26/86/902ee78db1bab1f0410f799869a49bb03b83be8d44c23b224d9db34f21c3/sympy-1.9.tar.gz"
aget "https://files.pythonhosted.org/packages/95/ba/7384cb4db4ed474d4582944053549e02ec25da630810e4a23454bc9fa617/mpmath-1.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/d4/cf/3965bddbb4f1a61c49aacae0e78fd1fe36b5dc36c797b31f30cf07dcbbb7/mpmath-1.2.1-py3-none-any.whl"
aget "https://files.pythonhosted.org/packages/b4/a2/4faa34bf0cdbefd5c706625f1234987795f368eb4e97bde9d6f46860843e/scipy-1.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/66/99/fc60e2287bb2309b8db4d0f080770ecc8d37dc64911e37b86698ec4b6a51/pybind11-2.7.1.tar.gz"
aget "https://files.pythonhosted.org/packages/c6/e6/986a967dcca91d89e36f4d4a2f69a052030bce01a7cd48a6b7fba1a50189/pythran-0.9.12.post1.tar.gz"
aget "https://files.pythonhosted.org/packages/e1/0a/a1e28f3c532655032a73482cb11ae63f8c38c8f58dc1a713b3791855fa68/pythran-0.9.12.post1-py3-none-any.whl"
aget "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
aget "https://files.pythonhosted.org/packages/a6/fb/7ff6a4ee66673c5964d3cf515ae85ba2076bc64bc2dcbbbd0153718b005f/gast-0.5.0.tar.gz"
aget "https://files.pythonhosted.org/packages/36/09/a4a6a967ca5bcfc0bd6162df4ee93017301fa7d9671483c849300bdba0db/beniget-0.4.0.tar.gz"
aget "https://files.pythonhosted.org/packages/47/b6/ea8a7728f096a597f0032564e8013b705aa992a0990becd773dcc4d7b4a7/psutil-5.9.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4d/19/e11fcc88288f68ae48e3aa9cf5a6fd092a88e629cb723465666c44d487a0/pep517-0.13.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4c/76/1e41fbb365ad20b6efab2e61b0f4751518444c953b390f9b2d36cf97eea0/Cython-0.29.32.tar.gz"
aget "https://files.pythonhosted.org/packages/64/8e/9929b64e146d240507edaac2185cd5516f00b133be5b39250d253be25a64/numpy-1.23.4.tar.gz"
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
    ( cd Grid-lehner ; git pull https://github.com/jinluchang/Grid.git )
else
    git clone https://github.com/jinluchang/Grid.git Grid-lehner
fi

if [ -e Grid-lehner/configure ] ; then
    echo "Grid-lehner bootstrapped."
else
    ( cd Grid-lehner ; git clean -f ; ./bootstrap.sh ; ls -l Eigen )
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
