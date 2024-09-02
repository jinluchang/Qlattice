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
, git
}:

buildPythonPackage rec {

  pname = "qlat_utils";
  version = "0.68";

  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    hash = "sha256-Dl5cmgWMEj8Xkog7+rPX6H17b6QZQUoiYFbij5b7jLs=";
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
  ];

  dependencies = [
    numpy
    psutil
    cython
  ];

}
