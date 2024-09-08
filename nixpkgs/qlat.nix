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
, is-pypi-src ? true
}:

buildPythonPackage rec {

  pname = "qlat";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.69";
  src-pypi = fetchPypi {
    inherit pname;
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-l/2kAChKFy/1Nnp6zx1OcUmewQi3C5Wl9RVeY0Q3wIY=";
  };

  version-local = "${../VERSION}";
  src-local = ../qlat;

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
