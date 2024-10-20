#!/usr/bin/env bash

source scripts/download-core.sh

cd "$distfiles"

dget "ninja-1.12.0.tar.gz" "https://github.com/ninja-build/ninja/archive/refs/tags/v1.12.0.tar.gz"

dget "lapack-3.12.0.tar.gz" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.0.tar.gz"

aget "https://github.com/tukaani-project/xz/releases/download/v5.4.6/xz-5.4.6.tar.gz"

aget "https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz"

aget "https://zlib.net/fossils/zlib-1.3.1.tar.gz"

aget "https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz"

aget "https://ftp.gnu.org/gnu/tar/tar-1.35.tar.gz"

aget "https://ftp.gnu.org/gnu/bison/bison-3.8.2.tar.xz"

aget "https://ftp.gnu.org/gnu/binutils/binutils-2.42.tar.xz"

aget "https://ftp.gnu.org/gnu/autoconf/autoconf-2.72.tar.xz"

aget "https://ftp.gnu.org/gnu/automake/automake-1.16.5.tar.gz"

aget "https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz"

aget "https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz"

aget "https://ftp.gnu.org/gnu/mpc/mpc-1.3.1.tar.gz"

aget "https://gcc.gnu.org/pub/gcc/infrastructure/isl-0.24.tar.bz2"

aget "https://ftp.gnu.org/gnu/gcc/gcc-13.2.0/gcc-13.2.0.tar.xz"

aget "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"

aget "https://github.com/Kitware/CMake/releases/download/v3.29.3/cmake-3.29.3.tar.gz"

# aget "https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.3.tar.bz2"
aget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.bz2"

aget "https://github.com/llvm/llvm-project/releases/download/llvmorg-18.1.5/llvm-project-18.1.5.src.tar.xz"

aget "https://www.cpan.org/src/5.0/perl-5.38.2.tar.gz"

aget "https://github.com/openssl/openssl/releases/download/openssl-3.3.0/openssl-3.3.0.tar.gz"

aget "https://www.python.org/ftp/python/3.12.3/Python-3.12.3.tar.xz"

aget "https://github.com/libffi/libffi/releases/download/v3.4.6/libffi-3.4.6.tar.gz"

aget "https://github.com/skvadrik/re2c/releases/download/3.1/re2c-3.1.tar.xz"

aget "https://sourceforge.net/projects/tclap/files/tclap-1.2.5.tar.gz"

# aget "https://sourceforge.net/projects/gnuplot/files/gnuplot/6.0.0/gnuplot-6.0.0.tar.gz"
aget "https://sourceforge.net/projects/gnuplot/files/gnuplot/5.4.8/gnuplot-5.4.8.tar.gz"

aget "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.27/OpenBLAS-0.3.27.tar.gz"

aget "https://sourceforge.net/projects/tcl/files/Tcl/8.6.14/tcl8.6.14-src.tar.gz"

aget "https://sourceforge.net/projects/tcl/files/Tcl/8.6.14/tk8.6.14-src.tar.gz"

aget "https://ftp.gnu.org/pub/gnu/ncurses/ncurses-6.5.tar.gz"

aget "https://github.com/libevent/libevent/releases/download/release-2.1.12-stable/libevent-2.1.12-stable.tar.gz"

aget "https://github.com/tmux/tmux/releases/download/3.4/tmux-3.4.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "ninja-10.0.tar.gz" "https://github.com/jinluchang/Ninja-dummy/archive/refs/tags/v10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/d0/9641dc7b05877874c6418f8034ddefc809495e65caa14d38c7551cd114bb/pip-24.1.1.tar.gz"
aget "https://files.pythonhosted.org/packages/41/8a/0d1bbd33cd3091c913d298746e56f40586fa954788f51b816c6336424675/sympy-1.12.1.tar.gz"
aget "https://files.pythonhosted.org/packages/d5/f7/2fdd9205a2eedee7d9b0abbf15944a1151eb943001dbdc5233b1d1cfc34e/Cython-3.0.10.tar.gz"
aget "https://files.pythonhosted.org/packages/4e/e5/0230da034a2e1b1feb32621d7cd57c59484091d6dccc9e6b855b0d309fc9/scipy-1.14.0.tar.gz"
aget "https://files.pythonhosted.org/packages/05/35/fb1ada118002df3fe91b5c3b28bc0d90f879b881a5d8f68b1f9b79c44bfe/numpy-2.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/73/32/f892675c5009cd4c1895ded3d6153476bf00adb5ad1634d03635620881f5/pythran-0.16.1.tar.gz"
aget "https://files.pythonhosted.org/packages/e0/47/dd32fa426cc72114383ac549964eecb20ecfd886d1e5ccf5340b55b02f57/mpmath-1.3.0.tar.gz"
aget "https://files.pythonhosted.org/packages/1a/3f/b19e9354c358f5acf322dd1f81ed9f0c633ba4bcccfd32e9c3740c43c9e5/meson_python-0.16.0.tar.gz"
aget "https://files.pythonhosted.org/packages/6f/ee/7ead9d4e69c94a61db0e2a24a5c7dbed6a1ae02b7876030971d76a9857cf/meson-1.4.1.tar.gz"
aget "https://files.pythonhosted.org/packages/44/d7/8f5d2be1a5fed3b0b5ccd3e800153c0f4dd84c2a688d25bce0bb0cb1f87f/pep517-0.13.1.tar.gz"
aget "https://files.pythonhosted.org/packages/b3/17/1d146e0127b66f1945251f130afac430985d2f9d75a3c0330355f21d876a/mpi4py-3.1.6.tar.gz"
aget "https://files.pythonhosted.org/packages/8d/e6/2fc95aec377988ff3ca882aa58d4f6ab35ff59a12b1611a9fe3075eb3019/setuptools-70.2.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4f/a4/00a9ac1b555294710d4a68d2ce8dfdf39d72aa4d769a7395d05218d88a42/setuptools_scm-8.1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/05/3b/23cb81e4cc567c1c4500c0f7ca865225d8cc2a06221099ff5826b99d4e4c/pybind11-2.12.0.tar.gz"
aget "https://files.pythonhosted.org/packages/14/e7/50cbac38f77eca8efd39516be6651fdb9f3c4c0fab8cf2cf05f612578737/beniget-0.4.1.tar.gz"
aget "https://files.pythonhosted.org/packages/c7/07/6f63dda440d4abb191b91dc383b472dae3dd9f37e4c1e4a5c3db150531c6/pyproject_hooks-1.1.0.tar.gz"
aget "https://files.pythonhosted.org/packages/cf/cc/428b057f8c229b7c374efe9d6a6a35e693f79e071e25846ab0c55e59d337/pyproject_metadata-0.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/fc/f8/98eea607f65de6527f8a2e8885fc8015d3e6f5775df186e443e0964a11c3/distro-1.9.0.tar.gz"
aget "https://files.pythonhosted.org/packages/ce/9e/2d725d2f7729c6e79ca62aeb926492abbc06e25910dd30139d60a68bcb19/build-1.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/b8/d6/ac9cd92ea2ad502ff7c1ab683806a9deb34711a1e2bd8a59814e8fc27e69/wheel-0.43.0.tar.gz"
aget "https://files.pythonhosted.org/packages/24/43/a0b5837cf30db1561a04187edd262bdefaffcb61222cb441eadef35f9103/scikit_build-0.18.0-py3-none-any.whl"
aget "https://files.pythonhosted.org/packages/ee/b5/b43a27ac7472e1818c4bafd44430e69605baefe1f34440593e0332ec8b4d/packaging-24.0.tar.gz"
aget "https://files.pythonhosted.org/packages/c4/e0/e05fee8b5425db6f83237128742e7e5ef26219b687ab8f0d41ed0422125e/pkgconfig-1.5.5.tar.gz"
aget "https://files.pythonhosted.org/packages/c4/e6/c1ac50fe3eebb38a155155711e6e864e254ce4b6e17fe2429b4c4d5b9e80/flit_core-3.9.0.tar.gz"
aget "https://files.pythonhosted.org/packages/f2/db/20a9f9cae3f3c213a8c406deb4395698459fd96962cea8f2ccb230b1943c/poetry_core-1.9.0.tar.gz"
aget "https://files.pythonhosted.org/packages/18/c7/8c6872f7372eb6a6b2e4708b88419fb46b857f7a2e1892966b851cc79fc9/psutil-6.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e4/41/f26f62ebef1a80148e20951a6e9ef4d0ebbe2090124bc143da26e12a934c/gast-0.5.4.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
aget "https://files.pythonhosted.org/packages/f6/f3/b827b3ab53b4e3d8513914586dcca61c355fa2ce8252dea4da56e67bf8f2/typing_extensions-4.11.0.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
)

cd "$wd"

./scripts/update-sources.sh
