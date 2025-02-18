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
, gvar
, vegas
# tests
, hypothesis
}:

buildPythonPackage rec {

  pname = "lsqfit";
  version = "13.2.2";
  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    hash = "sha256-y1wewo6gwWpjr6ss1UXGU48GS4xaKNaq0IVmFwQfsAo=";
  };

  build-system = [
    setuptools
    wheel
    numpy
    cython
    gvar
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    gvar
    vegas
  ];

  meta = {
    description = "Utilities for nonlinear least-squares fits.";
    homepage = "https://github.com/gplepage/lsqfit";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
