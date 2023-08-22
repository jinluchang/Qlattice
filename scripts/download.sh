#!/bin/bash

source scripts/download-core.sh

cd "$distfiles"

dget "ninja-1.11.1.tar.gz" "https://github.com/ninja-build/ninja/archive/refs/tags/v1.11.1.tar.gz"

dget "lapack-3.10.1.tar.gz" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.1.tar.gz"

aget "https://tukaani.org/xz/xz-5.4.1.tar.gz"

aget "https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz"

aget "https://zlib.net/zlib-1.2.13.tar.gz"

aget "https://github.com/facebook/zstd/releases/download/v1.5.2/zstd-1.5.2.tar.gz"

aget "https://ftp.gnu.org/gnu/tar/tar-1.34.tar.gz"

aget "https://ftp.gnu.org/gnu/binutils/binutils-2.38.tar.xz"

aget "https://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz"

aget "https://ftp.gnu.org/gnu/automake/automake-1.16.5.tar.gz"

aget "https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz"

aget "https://ftp.gnu.org/gnu/mpfr/mpfr-4.1.1.tar.xz"

aget "https://ftp.gnu.org/gnu/mpc/mpc-1.2.1.tar.gz"

aget "https://gcc.gnu.org/pub/gcc/infrastructure/isl-0.24.tar.bz2"

aget "https://ftp.gnu.org/gnu/gcc/gcc-11.3.0/gcc-11.3.0.tar.xz"

aget "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"

aget "https://github.com/Kitware/CMake/releases/download/v3.24.3/cmake-3.24.3.tar.gz"

aget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.bz2"

aget "https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.6/llvm-project-15.0.6.src.tar.xz"

aget "https://www.cpan.org/src/5.0/perl-5.34.0.tar.gz"

aget "https://www.openssl.org/source/openssl-3.0.7.tar.gz"

aget "https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tar.xz"

aget "https://github.com/libffi/libffi/releases/download/v3.4.4/libffi-3.4.4.tar.gz"

aget "https://github.com/skvadrik/re2c/releases/download/2.2/re2c-2.2.tar.xz"

aget "https://sourceforge.net/projects/tclap/files/tclap-1.2.5.tar.gz"

aget "https://sourceforge.net/projects/gnuplot/files/gnuplot/5.4.5/gnuplot-5.4.5.tar.gz"

aget "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz"

aget "https://prdownloads.sourceforge.net/tcl/tcl8.6.13-src.tar.gz"

aget "https://prdownloads.sourceforge.net/tcl/tk8.6.13-src.tar.gz"

aget "https://ftp.gnu.org/pub/gnu/ncurses/ncurses-6.4.tar.gz"

aget "https://github.com/libevent/libevent/releases/download/release-2.1.12-stable/libevent-2.1.12-stable.tar.gz"

aget "https://github.com/tmux/tmux/releases/download/3.3a/tmux-3.3a.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "ninja-1.0.tar.gz" "https://github.com/jinluchang/Ninja-dummy/archive/refs/tags/1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4d/19/e11fcc88288f68ae48e3aa9cf5a6fd092a88e629cb723465666c44d487a0/pep517-0.13.0.tar.gz"
aget "https://files.pythonhosted.org/packages/bb/f9/a526ea001ecadd48d90ed7a61038b5ba732136d76b60357d23dc37521c58/pybind11-2.10.4.tar.gz"
aget "https://files.pythonhosted.org/packages/e0/47/dd32fa426cc72114383ac549964eecb20ecfd886d1e5ccf5340b55b02f57/mpmath-1.3.0.tar.gz"
aget "https://files.pythonhosted.org/packages/d6/0f/96b7309212a926c1448366e9ce69b081ea79d63265bde33f11cc9cfc2c07/psutil-5.9.5.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
aget "https://files.pythonhosted.org/packages/48/a3/0bd844c54ae8141642088b7ae09dd38fec2ec7faa9b7d25bb6a23c1f266f/gast-0.5.3.tar.gz"
aget "https://files.pythonhosted.org/packages/38/af/b0e6a9eba989870fd26e10889446d1bec2e6d5be0a1bae2dc4dcda9ce199/pyproject-metadata-0.7.1.tar.gz"
aget "https://files.pythonhosted.org/packages/14/e7/50cbac38f77eca8efd39516be6651fdb9f3c4c0fab8cf2cf05f612578737/beniget-0.4.1.tar.gz"
aget "https://files.pythonhosted.org/packages/b9/6c/7c6658d258d7971c5eb0d9b69fa9265879ec9a9158031206d47800ae2213/packaging-23.1.tar.gz"
aget "https://files.pythonhosted.org/packages/d3/20/06270dac7316220643c32ae61694e451c98f8caf4c8eab3aa80a2bedf0df/typing_extensions-4.5.0.tar.gz"
aget "https://files.pythonhosted.org/packages/71/22/207523d16464c40a0310d2d4d8926daffa00ac1f5b1576170a32db749636/pyparsing-3.0.9.tar.gz"
aget "https://files.pythonhosted.org/packages/4b/89/eaa3a3587ebf8bed93e45aa79be8c2af77d50790d15b53f6dfc85b57f398/distro-1.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/10/e5/be08751d07b30889af130cec20955c987a74380a10058e6e8856e4010afc/flit_core-3.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
aget "https://files.pythonhosted.org/packages/25/c1/374304b8407d3818f7025457b7366c8e07768377ce12edfe2aa58aa0f64c/pyproject_hooks-1.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/de/1c/fb62f81952f0e74c3fbf411261d1adbdd2d615c89a24b42d0fe44eb4bcf3/build-0.10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/d2/a3/01d72e506afcf919c191de587c3f77ee7355b0b28951eb1fbc2753bc9d77/scikit_build-0.17.3-py3-none-any.whl"
aget "https://files.pythonhosted.org/packages/98/12/2c1e579bb968759fc512391473340d0661b1a8c96a59fb7c65b02eec1321/setuptools_scm-7.1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/a8/e7/1440b0d19054a5616e9e5beeaa22f68485aa9de20d187f04e52880b7ae7a/setuptools-59.2.0.tar.gz"
aget "https://files.pythonhosted.org/packages/fc/ef/0335f7217dd1e8096a9e8383e1d472aa14717878ffe07c4772e68b6e8735/wheel-0.40.0.tar.gz"
aget "https://files.pythonhosted.org/packages/32/d7/854e45d2b03e1a8ee2aa6429dd396d002ce71e5d88b77551b2fb249cb382/versioneer-0.29.tar.gz"
aget "https://files.pythonhosted.org/packages/cb/1c/af7f886e723b2dfbaea9b8a739153f227b386dd856cf956f9fd0ed0a502b/poetry_core-1.7.0.tar.gz"
aget "https://files.pythonhosted.org/packages/c4/e0/e05fee8b5425db6f83237128742e7e5ef26219b687ab8f0d41ed0422125e/pkgconfig-1.5.5.tar.gz"
aget "https://files.pythonhosted.org/packages/c5/40/7cf58e6230f0e76699f011c6d293dd47755997709a303a4e644823f3a753/h5py-3.7.0.tar.gz"
aget "https://files.pythonhosted.org/packages/ba/19/e63fb4e0d20e48bd2167bb7e857abc0e21679e24805ba921a224df8977c0/pip-23.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/f2/21/3326fe66aa091b3653d12affd1c3928e17c9ac386708ec42b4fad87eefac/meson-1.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/ed/77/786ac8dcc8bc39a927527ba68a016bf9bd8f7ffe5c3622597ad16cd218af/meson_python-0.13.2.tar.gz"
aget "https://files.pythonhosted.org/packages/bc/f2/749af7fd0e7703ddca6cea525ab40f26c3ca6cbe6c23658441c6f9705860/mpi4py-3.1.4.tar.gz"
aget "https://files.pythonhosted.org/packages/7f/a2/fd5ced5dd33597ef291861bfadd46820de417b41bcb6ca2fa0b5f6fa8152/Cython-3.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/57/3485a1a3dff51bfd691962768b14310dae452431754bfc091250be50dd29/sympy-1.12.tar.gz"
aget "https://files.pythonhosted.org/packages/8d/d8/b27e8dc3f3a03dcd317d40d9df0ae07ebbd85444585973ceba07716934d0/pythran-0.13.1.tar.gz"
aget "https://files.pythonhosted.org/packages/a0/41/8f53eff8e969dd8576ddfb45e7ed315407d27c7518ae49418be8ed532b07/numpy-1.25.2.tar.gz"
aget "https://files.pythonhosted.org/packages/9c/ef/87a5565907645998d7c62e76b84b0ca9f0b7c25cd433f5617a968051cec3/scipy-1.11.2.tar.gz"
aget "https://files.pythonhosted.org/packages/b1/a7/824332581e258b5aa4f3763ecb2a797e5f9a54269044ba2e50ac19936b32/pandas-2.0.3.tar.gz"
aget "https://files.pythonhosted.org/packages/34/80/a49e6119ede0ea77ccc0e043338a28dc6baf5c59e75f169d46b87f28c4a3/xarray-2023.8.0.tar.gz"
)

cd "$wd"

./scripts/update-sources.sh
