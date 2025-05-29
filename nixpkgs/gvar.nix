{ lib
, buildPythonPackage
, fetchPypi
# build-system
, setuptools
, wheel
, numpy
, cython
, oldest-supported-numpy
# dependencies
, scipy
# tests
}:

buildPythonPackage rec {

  env.NIX_CFLAGS_COMPILE = lib.concatStringsSep " " [
    "-Wno-error=int-conversion"
    "-Wno-error=incompatible-pointer-types"
    "-Wno-error=implicit-function-declaration"
    ];

  pname = "gvar";
  version = "13.1.6";
  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    hash = "sha256-U9DAxHrpY2Z19eFdxn2x7YAn6I3/phgkviuDbmbdmPI=";
  };

  build-system = [
    setuptools
    wheel
    numpy
    cython
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    scipy
  ];

  meta = {
    description = "Utilities for manipulating correlated Gaussian random variables.";
    homepage = "https://github.com/gplepage/gvar";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
