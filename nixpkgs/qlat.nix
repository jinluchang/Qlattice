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

  pname = "qlat";
  version = "0.68";

  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    hash = "sha256-JwPJoMo4EQlNT+6h4PQRLgAR6SYkd5Cxisk1HeSxfic=";
  };

  enableParallelBuilding = true;

  build-system = [
    meson-python
    pkg-config
    cython
    numpy
    qlat_utils
  ];

  nativeBuildInputs = [
    git
  ];

  buildInputs = [
    mpi
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
    mpi4py
    sympy
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
