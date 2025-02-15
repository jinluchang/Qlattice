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
# tests
, hypothesis
}:

buildPythonPackage rec {

  pname = "vegas";
  version = "6.1.2";
  pyproject = true;

  src = fetchPypi {
    inherit pname version;
    hash = "sha256-Q5ZOZYYZgIbvxMBE953i12EDOiNe61NMEJr7240ocIw=";
  };

  build-system = [
    setuptools
    numpy
    cython
    oldest-supported-numpy
  ];

  dependencies = [
    numpy
    gvar
  ];

  meta = {
    description = "Tools for adaptive multidimensional Monte Carlo integration.";
    homepage = "https://github.com/gplepage/vegas";
    license = lib.licenses.gpl3Plus;
    maintainers = with lib.maintainers; [];
  };

}
