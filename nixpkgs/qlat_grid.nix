{ lib
, stdenv
, fetchPypi
, python
, buildPythonPackage
, setuptools
, setuptools-scm
, cython
, meson-python
, pkg-config
, numpy
, mpi4py
, sympy
, qlat_utils
, qlat
, grid-lehner
, hdf5
, openssl
, gmp
, mpfr
, c-lime
, psutil
, zlib
, eigen
, cuba
, mpi
, git
, fftw
, fftwFloat
, gsl
, is-pypi-src ? true
}:

buildPythonPackage rec {

  pname = "qlat_grid";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.69";
  src-pypi = fetchPypi {
    inherit pname;
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-1OGzJEin1HeHHBmt0ZbrY6V9OyDM2RQjMqp66GeuhWc=";
  };

  version-local = "${../VERSION}";
  src-local = ../qlat-grid;

  enableParallelBuilding = true;

  build-system = [
    meson-python
    pkg-config
    cython
    numpy
    qlat_utils
    qlat
  ];

  nativeBuildInputs = [
    git
    grid-lehner
  ];

  buildInputs = [
    mpi
    grid-lehner
    hdf5
    openssl
    gmp
    mpfr
    c-lime
    zlib
    eigen
    cython
    fftw
    fftwFloat
    gsl
    cuba
  ];

  dependencies = [
    mpi
    eigen
    cuba
    cython
    numpy
    psutil
    qlat_utils
    qlat
    grid-lehner
    mpi4py
    sympy
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
