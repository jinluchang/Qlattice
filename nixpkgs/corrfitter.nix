{ lib
, buildPythonPackage
, fetchPypi
# build-system
, setuptools
, setuptools-scm
, lsqfit
, gvar
, numpy
# dependencies
# tests
, hypothesis
}:

buildPythonPackage rec {

  pname = "corrfitter";
  version = "8.2";
  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    hash = "sha256-P5qiqeY3raams0vaAZnMI93rIAV7NWxlfXMB1Se85n4=";
  };

  build-system = [
    setuptools
    setuptools-scm
    lsqfit
    gvar
    numpy
  ];

  dependencies = [
    lsqfit
    gvar
    numpy
  ];

  meta = {
    description = "Utilities for fitting correlators in lattice QCD.";
    homepage = "https://github.com/gplepage/corrfitter";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
