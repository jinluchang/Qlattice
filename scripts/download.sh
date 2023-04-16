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

aget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.bz2"

aget "https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.6/llvm-project-15.0.6.src.tar.xz"

aget "https://www.cpan.org/src/5.0/perl-5.34.0.tar.gz"

aget "https://www.openssl.org/source/openssl-3.0.7.tar.gz"

aget "https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tar.xz"

aget "https://github.com/libffi/libffi/releases/download/v3.4.4/libffi-3.4.4.tar.gz"

aget "https://github.com/skvadrik/re2c/releases/download/2.2/re2c-2.2.tar.xz"

aget "https://sourceforge.net/projects/tclap/files/tclap-1.2.5.tar.gz"

aget "https://sourceforge.net/projects/gnuplot/files/gnuplot/5.4.5/gnuplot-5.4.5.tar.gz"

aget "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.21/OpenBLAS-0.3.21.tar.gz"

aget "https://prdownloads.sourceforge.net/tcl/tcl8.6.13-src.tar.gz"

aget "https://prdownloads.sourceforge.net/tcl/tk8.6.13-src.tar.gz"

aget "https://ftp.gnu.org/pub/gnu/ncurses/ncurses-6.4.tar.gz"

aget "https://github.com/libevent/libevent/releases/download/release-2.1.12-stable/libevent-2.1.12-stable.tar.gz"

aget "https://github.com/tmux/tmux/releases/download/3.3/tmux-3.3.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "ninja-1.0.tar.gz" "https://github.com/jinluchang/Ninja-dummy/archive/refs/tags/1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/da/bf/1bdbe62f5fbde085351693e3a8e387a59f8220932b911b1719fe65efa2d7/pip-23.1.tar.gz"
aget "https://files.pythonhosted.org/packages/4d/19/e11fcc88288f68ae48e3aa9cf5a6fd092a88e629cb723465666c44d487a0/pep517-0.13.0.tar.gz"
aget "https://files.pythonhosted.org/packages/2f/23/aaf147a5bc31c8be286f07d862b3699d7b49e3411fb75087525b5c31ab3e/pybind11-2.10.1.tar.gz"
aget "https://files.pythonhosted.org/packages/95/ba/7384cb4db4ed474d4582944053549e02ec25da630810e4a23454bc9fa617/mpmath-1.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/3d/7d/d05864a69e452f003c0d77e728e155a89a2a26b09e64860ddd70ad64fb26/psutil-5.9.4.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
aget "https://files.pythonhosted.org/packages/48/a3/0bd844c54ae8141642088b7ae09dd38fec2ec7faa9b7d25bb6a23c1f266f/gast-0.5.3.tar.gz"
aget "https://files.pythonhosted.org/packages/8e/e3/4314ed332414c059aed340fe2097d5ca065cc8d165c08554c83eb96a59c4/pyproject-metadata-0.6.1.tar.gz"
aget "https://files.pythonhosted.org/packages/14/e7/50cbac38f77eca8efd39516be6651fdb9f3c4c0fab8cf2cf05f612578737/beniget-0.4.1.tar.gz"
aget "https://files.pythonhosted.org/packages/6b/f7/c240d7654ddd2d2f3f328d8468d4f1f876865f6b9038b146bec0a6737c65/packaging-22.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e3/a7/8f4e456ef0adac43f452efc2d0e4b242ab831297f1bac60ac815d37eb9cf/typing_extensions-4.4.0.tar.gz"
aget "https://files.pythonhosted.org/packages/71/22/207523d16464c40a0310d2d4d8926daffa00ac1f5b1576170a32db749636/pyparsing-3.0.9.tar.gz"
aget "https://files.pythonhosted.org/packages/4b/89/eaa3a3587ebf8bed93e45aa79be8c2af77d50790d15b53f6dfc85b57f398/distro-1.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/10/e5/be08751d07b30889af130cec20955c987a74380a10058e6e8856e4010afc/flit_core-3.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
aget "https://files.pythonhosted.org/packages/83/40/c0ea8d3072441403aa30d91b900051a961e0f6d58f702c0ec9ca812c8737/meson-1.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/09/51/600a029573c775081d5145fd1b09d22a8d9e73671d84792534f491f38fef/meson_python-0.12.0.tar.gz"
aget "https://files.pythonhosted.org/packages/d6/bd/2d13a273d95f7b7d9903c906c486040b0aebb85e008f93a5dd0891f21f1f/scipy-1.10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/5a/36/4667b08bc45131fe655a27b1a112c1730f3244343c53a338f44d730bd6ba/sympy-1.11.1.tar.gz"
aget "https://files.pythonhosted.org/packages/25/c1/374304b8407d3818f7025457b7366c8e07768377ce12edfe2aa58aa0f64c/pyproject_hooks-1.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/de/1c/fb62f81952f0e74c3fbf411261d1adbdd2d615c89a24b42d0fe44eb4bcf3/build-0.10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/00/91/2c7cd1a6b567e8e40c1677c61aba39e5f14ccf06ac78ae4fa6acb84ae140/scikit-build-0.16.6.tar.gz"
aget "https://files.pythonhosted.org/packages/98/12/2c1e579bb968759fc512391473340d0661b1a8c96a59fb7c65b02eec1321/setuptools_scm-7.1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/a8/e7/1440b0d19054a5616e9e5beeaa22f68485aa9de20d187f04e52880b7ae7a/setuptools-59.2.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4e/be/8139f127b4db2f79c8b117c80af56a3078cc4824b5b94250c7f81a70e03b/wheel-0.37.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e4/a9/6704bb5e1d1d778d3a6ee1278a8d8134f0db160e09d52863a24edb58eab5/numpy-1.24.2.tar.gz"
aget "https://files.pythonhosted.org/packages/4a/1b/059a68158bf65c857cfd6b80aed06a8fd35f2582cf548fb96f0b519b0d2b/pythran-0.12.1.tar.gz"
aget "https://files.pythonhosted.org/packages/dc/f6/e8e302f9942cbebede88b1a0c33d0be3a738c3ac37abae87254d58ffc51c/Cython-0.29.33.tar.gz"
aget "https://files.pythonhosted.org/packages/c5/40/7cf58e6230f0e76699f011c6d293dd47755997709a303a4e644823f3a753/h5py-3.7.0.tar.gz"
aget "https://files.pythonhosted.org/packages/bc/f2/749af7fd0e7703ddca6cea525ab40f26c3ca6cbe6c23658441c6f9705860/mpi4py-3.1.4.tar.gz"
)

cd "$wd"
