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
}:

buildPythonPackage rec {

  pname = "qlat_grid";
  version = "0.69";

  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    hash = "sha256-1OGzJEin1HeHHBmt0ZbrY6V9OyDM2RQjMqp66GeuhWc=";
  };

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
