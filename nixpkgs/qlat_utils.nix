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
, psutil
, zlib
, eigen
, cuba
, git
}:

buildPythonPackage rec {

  pname = "qlat_utils";
  version = "0.69";

  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    hash = "sha256-gBiJ9ilpzXnkfWrU7wazHNePT4vc1qXouaI9InYcRIg=";
  };

  enableParallelBuilding = true;

  build-system = [
    meson-python
    pkg-config
    cython
    numpy
  ];

  nativeBuildInputs = [
    git
  ];

  buildInputs = [
    numpy
    zlib
    eigen
    cuba
  ];

  dependencies = [
    cuba
    cython
    numpy
    psutil
  ];

  postPatch = ''
    sed -i "s/'-j4'/'-j$NIX_BUILD_CORES'/" pyproject.toml
  '';

}
