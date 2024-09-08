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
, qmp
, qio
, cps
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

  pname = "qlat_cps";
  version = "0.69";

  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    hash = "sha256-QQM7k3y+q9K3DDU9JAlZMUcVJ4me0b/8llPxh7VjZ+Y=";
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
  ];

  buildInputs = [
    mpi
    cps
    qmp
    qio
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
    mpi4py
    sympy
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
