#!/bin/bash

source scripts/download-core.sh

cd "$distfiles"

dget "ninja-1.11.1.tar.gz" "https://github.com/ninja-build/ninja/archive/refs/tags/v1.11.1.tar.gz"

dget "lapack-3.11.0.tar.gz" "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.11.0.tar.gz"

aget "https://tukaani.org/xz/xz-5.4.4.tar.gz"

aget "https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz"

aget "https://zlib.net/fossils/zlib-1.3.tar.gz"

aget "https://github.com/facebook/zstd/releases/download/v1.5.5/zstd-1.5.5.tar.gz"

aget "https://ftp.gnu.org/gnu/tar/tar-1.35.tar.gz"

aget "https://ftp.gnu.org/gnu/bison/bison-3.8.2.tar.xz"

aget "https://ftp.gnu.org/gnu/binutils/binutils-2.41.tar.xz"

aget "https://ftp.gnu.org/gnu/autoconf/autoconf-2.71.tar.gz"

aget "https://ftp.gnu.org/gnu/automake/automake-1.16.5.tar.gz"

aget "https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz"

aget "https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz"

aget "https://ftp.gnu.org/gnu/mpc/mpc-1.3.1.tar.gz"

aget "https://gcc.gnu.org/pub/gcc/infrastructure/isl-0.24.tar.bz2"

aget "https://ftp.gnu.org/gnu/gcc/gcc-12.3.0/gcc-12.3.0.tar.xz"

aget "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"

aget "https://github.com/Kitware/CMake/releases/download/v3.27.6/cmake-3.27.6.tar.gz"

aget "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.bz2"

aget "https://github.com/llvm/llvm-project/releases/download/llvmorg-17.0.1/llvm-project-17.0.1.src.tar.xz"

aget "https://www.cpan.org/src/5.0/perl-5.38.0.tar.gz"

aget "https://github.com/openssl/openssl/releases/download/openssl-3.1.3/openssl-3.1.3.tar.gz"

aget "https://www.python.org/ftp/python/3.11.5/Python-3.11.5.tar.xz"

aget "https://github.com/libffi/libffi/releases/download/v3.4.4/libffi-3.4.4.tar.gz"

aget "https://github.com/skvadrik/re2c/releases/download/3.1/re2c-3.1.tar.xz"

aget "https://sourceforge.net/projects/tclap/files/tclap-1.2.5.tar.gz"

aget "https://sourceforge.net/projects/gnuplot/files/gnuplot/5.4.8/gnuplot-5.4.8.tar.gz"

aget "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.24/OpenBLAS-0.3.24.tar.gz"

aget "https://sourceforge.net/projects/tcl/files/Tcl/8.6.13/tcl8.6.13-src.tar.gz"

aget "https://sourceforge.net/projects/tcl/files/Tcl/8.6.13/tk8.6.13-src.tar.gz"

aget "https://ftp.gnu.org/pub/gnu/ncurses/ncurses-6.4.tar.gz"

aget "https://github.com/libevent/libevent/releases/download/release-2.1.12-stable/libevent-2.1.12-stable.tar.gz"

aget "https://github.com/tmux/tmux/releases/download/3.3a/tmux-3.3a.tar.gz"

(
mkdir -p python-packages
cd python-packages
dget "ninja-10.0.tar.gz" "https://github.com/jinluchang/Ninja-dummy/archive/refs/tags/v10.0.tar.gz"
aget "https://files.pythonhosted.org/packages/14/e7/50cbac38f77eca8efd39516be6651fdb9f3c4c0fab8cf2cf05f612578737/beniget-0.4.1.tar.gz"
aget "https://files.pythonhosted.org/packages/25/c1/374304b8407d3818f7025457b7366c8e07768377ce12edfe2aa58aa0f64c/pyproject_hooks-1.0.0.tar.gz"
aget "https://files.pythonhosted.org/packages/3f/aa/1a5c72615e0ba4dc30cc36de7e8a9a2eca2158922b0677654fced0d3476c/Cython-3.0.4.tar.gz"
aget "https://files.pythonhosted.org/packages/38/af/b0e6a9eba989870fd26e10889446d1bec2e6d5be0a1bae2dc4dcda9ce199/pyproject-metadata-0.7.1.tar.gz"
aget "https://files.pythonhosted.org/packages/3a/cc/903bb18de90b5d6e15379c97175371ac6414795d94b9c2f6468a9c1303aa/pybind11-2.11.1.tar.gz"
aget "https://files.pythonhosted.org/packages/4b/89/eaa3a3587ebf8bed93e45aa79be8c2af77d50790d15b53f6dfc85b57f398/distro-1.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/4d/19/e11fcc88288f68ae48e3aa9cf5a6fd092a88e629cb723465666c44d487a0/pep517-0.13.0.tar.gz"
aget "https://files.pythonhosted.org/packages/2c/ab/a647b8cc3ac1aa07cde06875157696e4522958fb8363474bce21c302d4d8/pythran-0.14.0.tar.gz"
aget "https://files.pythonhosted.org/packages/98/e3/83a89a9d338317f05a68c86a2bbc9af61235bc55a0c6a749d37598fb2af1/build-1.0.3.tar.gz"
aget "https://files.pythonhosted.org/packages/39/7b/9f265b7f074195392e893a5cdc66116c2f7a31fd5f3d9cceff661ec6df82/scipy-1.11.3.tar.gz"
aget "https://files.pythonhosted.org/packages/78/23/f78fd8311e0f710fe1d065d50b92ce0057fe877b8ed7fd41b28ad6865bfc/numpy-1.26.1.tar.gz"
aget "https://files.pythonhosted.org/packages/a4/99/78c4f3bd50619d772168bec6a0f34379b02c19c9cced0ed833ecd021fd0d/wheel-0.41.2.tar.gz"
aget "https://files.pythonhosted.org/packages/ef/cc/93f7213b2ab5ed383f98ce8020e632ef256b406b8569606c3f160ed8e1c9/setuptools-68.2.2.tar.gz"
aget "https://files.pythonhosted.org/packages/af/54/0a75c590c1aa1908c1036de783c0b7136d8fe4beb8ce4c52e4d92b9e9ca7/setuptools-scm-8.0.3.tar.gz"
aget "https://files.pythonhosted.org/packages/b9/6c/7c6658d258d7971c5eb0d9b69fa9265879ec9a9158031206d47800ae2213/packaging-23.1.tar.gz"
aget "https://files.pythonhosted.org/packages/ba/19/e63fb4e0d20e48bd2167bb7e857abc0e21679e24805ba921a224df8977c0/pip-23.2.1.tar.gz"
aget "https://files.pythonhosted.org/packages/bc/f2/749af7fd0e7703ddca6cea525ab40f26c3ca6cbe6c23658441c6f9705860/mpi4py-3.1.4.tar.gz"
aget "https://files.pythonhosted.org/packages/c4/e0/e05fee8b5425db6f83237128742e7e5ef26219b687ab8f0d41ed0422125e/pkgconfig-1.5.5.tar.gz"
aget "https://files.pythonhosted.org/packages/c4/e6/c1ac50fe3eebb38a155155711e6e864e254ce4b6e17fe2429b4c4d5b9e80/flit_core-3.9.0.tar.gz"
aget "https://files.pythonhosted.org/packages/cb/1c/af7f886e723b2dfbaea9b8a739153f227b386dd856cf956f9fd0ed0a502b/poetry_core-1.7.0.tar.gz"
aget "https://files.pythonhosted.org/packages/d6/0f/96b7309212a926c1448366e9ce69b081ea79d63265bde33f11cc9cfc2c07/psutil-5.9.5.tar.gz"
aget "https://files.pythonhosted.org/packages/e0/47/dd32fa426cc72114383ac549964eecb20ecfd886d1e5ccf5340b55b02f57/mpmath-1.3.0.tar.gz"
aget "https://files.pythonhosted.org/packages/e4/41/f26f62ebef1a80148e20951a6e9ef4d0ebbe2090124bc143da26e12a934c/gast-0.5.4.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/57/3485a1a3dff51bfd691962768b14310dae452431754bfc091250be50dd29/sympy-1.12.tar.gz"
aget "https://files.pythonhosted.org/packages/e5/69/882ee5c9d017149285cab114ebeab373308ef0f874fcdac9beb90e0ac4da/ply-3.11.tar.gz"
aget "https://files.pythonhosted.org/packages/38/a7/ddc350902a1b3b960db8d0e501f61468f925f994e0b4e6d696aeb6a75c00/meson_python-0.14.0.tar.gz"
aget "https://files.pythonhosted.org/packages/f6/36/948280852b1fba7525e7579a5b0e26a9d16fafacb8d841322fc3095e7706/meson-1.2.2.tar.gz"
aget "https://files.pythonhosted.org/packages/fa/af/b3ef8fe0bb96bf7308e1f9d196fc069f0c75d9c74cfaad851e418cc704f4/scikit_build-0.17.6-py3-none-any.whl"
aget "https://files.pythonhosted.org/packages/1f/7a/8b94bb016069caa12fc9f587b28080ac33b4fbb8ca369b98bc0a4828543e/typing_extensions-4.8.0.tar.gz"
aget "https://files.pythonhosted.org/packages/c0/3f/d7af728f075fb08564c5949a9c95e44352e23dee646869fa104a3b2060a3/tomli-2.0.1.tar.gz"
)

cd "$wd"

./scripts/update-sources.sh
