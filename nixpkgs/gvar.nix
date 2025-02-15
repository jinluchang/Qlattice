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

  pname = "gvar";
  version = "13.1";
  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    hash = "sha256-CGKB/NvBLhEF7VA6uSjt4irxaqBfu2z97vVOIOLwTts=";
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
