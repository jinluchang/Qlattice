{ fetchPypi
, buildPythonPackage
, meson-python
, pkg-config
, mpi4py
, sympy
, scipy
, jax
, jaxlib
, qlat_utils
, mpi
, git
, fftw
, fftwFloat
, gsl
, cuba
, is-pypi-src ? true
}:

buildPythonPackage rec {

  pname = "qlat";
  version = if is-pypi-src then version-pypi else version-local;

  pyproject = true;

  src = if is-pypi-src then src-pypi else src-local;

  version-pypi = "0.70";
  src-pypi = fetchPypi {
    inherit pname;
    version = version-pypi;
    extension = "tar.gz";
    hash = "sha256-Jm+FcqUt4A5jdlFGHvKBdqNsUa3zU1fNRnWhfWdzDUs=";
  };

  version-local = "${../VERSION}-current";
  src-local = ../qlat;

  enableParallelBuilding = true;

  build-system = [
    meson-python
    pkg-config
    qlat_utils
  ];

  nativeBuildInputs = [
    git
  ];

  propagatedBuildInputs = [
    mpi
    fftw
    fftwFloat
    gsl
    cuba
  ];

  dependencies = [
    qlat_utils
    mpi4py
    sympy
    scipy
    jax
    jaxlib
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
